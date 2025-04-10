#' Calibrate Bertrand model with GNL demand
#'
#' @param param Vector of demand parameters (alpha,mu)
#' @param price Observed prices
#' @param ownership Ownership matrix
#' @param share Observed market shares
#' @param cost Marginal costs for each product
#' @param weight Vector of length four with weights given to prices, shares,
#' diversions, and costs, respectively. Default is c(1,1,1,1).
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param div_matrix A matrix of observed diversions from product in row j to
#' product in column k.
#' @param mu_constraint_matrix is a (K-by-K') matrix indicating which nesting
#' parameters are constrained to be equal to each other, where K is the
#' number of nests and K' is the number of freely varying nesting
#' parameters. mu_full = mu_constraint_matrix %*% mu_prime. Where mu_full
#' is a vector of length K of the nesting parameter value for each nest,
#' and mu_prime is a vector of length K' of parameters to be calculated.
#' It must be the case that K is greater than K'.
#' @param div_calc_marginal is a logical if function should match to marginal
#' diversions (if TRUE) or second choice diversions (if FALSE). Default
#' to TRUE.
#' @param optimizer Which optimization routine should be used to find
#' equilibrium prices, either BBoptim or multiroot
#' @param optimizer_c Which optimization routine should be used to calibrate
#' costs, either BBoptim or optim
#' @param returnOutcomes logical; should equilibrium objects be returned (mean
#' value parameter, prices, shares, costs) as a list.
#'
#' @returns Difference between model predicted and observed values of
#' prices, shares, and diversions.
#'
#' @details This function calibrates a Bertrand model with generalized nested
#' logit (GNL) demand
#'
#' @examples
#' TO BE ADDED.
#'
#' @export


## useOldWeight is a legacy option in case want to use old weighting in the
## objective function

##################################################################
# Bertrand model first-order conditions for calibration with GNL demand
##################################################################


bertrand_calibrate_gnl <- function(param,ownership,price,shares,cost,
                                   weight = c(1,1,1,1),
                                   nest_allocation,div_matrix,
                                   mu_constraint_matrix = NA,
                                   div_calc_marginal = TRUE,
                                   optimizer="BBoptim",
                                   optimizer_c="optim",
                                   useOldWeight = FALSE,
                                   returnOutcomes = FALSE,
                                   maxitval = 1500,
                                   maxitval_c = 1500){

  # If GNL, define GNL objects
  a_jk <- nest_allocation
  B <- 1*(a_jk > 0)
  K_val <- dim(a_jk)[2]

  alpha <- param[1]


  #### checks on weighting vector ####
  if (useOldWeight == FALSE) {
  if (length(weight) != 4) {
    warning("Weight vector should be of length 4.")
  }
  }


  #### checks on mu_constraint matrix ####

  mu_prime <- param[2:length(param)]
  K_prime <- length(mu_prime)

  ## Intended nesting constraints partly
  ## implied by length of param. Check to make sure consistent, and if not
  ## then throw error that more information needs to be supplied.

  ## give error if K_prime > K_val. This maybe should be a stop()
  if (K_prime > K_val) {warning("K' should not be greater than K")}

  ## if no matrix provided, but K_prime implied by length of param is 1,
  ## then we can assume just one nesting parameter
  if (anyNA(mu_constraint_matrix) & K_prime == 1 ) {
    mu_constraint_matrix <- matrix(1, nrow = length(mu_prime), ncol = 1)
  }

  ## if no matrix provided, but K_prime implied by length of param is K,
  ## then we can assume full flexibility intended
  if (anyNA(mu_constraint_matrix) & K_prime == K_val ) {
    mu_constraint_matrix <- diag(K_val)
  }

  ## if still no mcm determined, give error that more information is need
  if (anyNA(mu_constraint_matrix) ) {
    warning("Please provide more information on nesting parameter calibration
            in mu_constraint_matrix")
  }


  mu <- mu_constraint_matrix %*% mu_prime

  J <- length(price)

  delta0 <- log(shares) - log(1-sum(shares)) - alpha*price  ## assuming alpha<0.

  find_d <- multiroot(f = match_share, start = delta0,
                      price=price,alpha=alpha,nest_allocation=a_jk,mu=mu,
                      shares_obs = shares)

  delta <- find_d$root

  # sometimes, multiroot fails inside of optimization. Set to delta0.
  # Eventually need to understand better when multiroot fails.

  #useOld <- TRUE  # if want to use old version of function, which had no correction
  useOld <- FALSE # if want to use new version, with delta NA correction

  if (useOld == FALSE) {
    if (anyNA(delta)) {
      delta <- delta0
      warning("Multiroot failed to find mean values that matched shares.")
    }
  }

  x0 <- price

  ## NEW: Given observed prices, and current guess of demand parameter values,
  ## back out costs consistent with FOCs.
  #x06 <- price * 0.5
  ## Add check that at least one element of the cost vector is non-missing.
  ## Or, adjust this to allow for starting value when cost vector is completely
  ## missing.

  if (anyNA(cost) == FALSE) {
    x06 <- cost
  }
  if (anyNA(cost) == TRUE) {
    meancost <- mean(cost, na.rm = TRUE)
    x06 <- ifelse(!is.na(cost), cost, meancost)
  }

  if (optimizer_c == "optim") {
    out_cost <- optim(f = bertrand_foc_c, par = x06,
                      price = price, own = ownership, alpha = alpha,
                      delta = delta,
                      nest_allocation=a_jk, mu=mu, sumFOC = TRUE,
                      control = list(maxit = maxitval_c) )

    cost_cal <- out_cost$par

    if (out_cost$convergence != 0) {
      warning(paste0("Cost calibration did not converge with code "),out_cost$convergence)
    }
  }

  if (optimizer_c == "BBoptim") {
    out_cost <- BBoptim(fn = bertrand_foc_c, par = x06,
                      price = price, own = ownership, alpha = alpha,
                      delta = delta,
                      nest_allocation=a_jk, mu=mu, sumFOC = TRUE,
                      control = list(maxit = maxitval_c) )

    cost_cal <- out_cost$par

    if (out_cost$convergence != 0) {
      warning(paste0("Cost calibration did not converge with code "),out_cost$convergence)
    }
  }

  ## BBoptim or multiroot. If there are missing costs, use BBoptim
  if (optimizer == "BBoptim") {
    out1 <- BBoptim(f = bertrand_foc, par = x0,
                    own = ownership, alpha= alpha,
                    delta = delta, cost = cost_cal,
                    nest_allocation = a_jk, mu = mu,
                    sumFOC = TRUE, control = list(trace=FALSE, maxit = maxitval))

    p_model <- out1$par

    if (out1$convergence != 0) {
      warning(paste0("Price equilibrium did not converge with code "),out1$convergence)
    }

  }
  if (optimizer == "multiroot") {
    out1 <- multiroot(f = bertrand_foc, start = x0,
                      own = ownership, alpha= alpha,
                      delta = delta, cost = cost_cal,
                      nest_allocation = a_jk, mu = mu)

    p_model <- out1$root
  }

  share_m <- share_calc(price=price,delta=delta,alpha=alpha,nest_allocation=a_jk,mu=mu)
  diversions_m <- diversion_calc(price=price,alpha=alpha,delta=delta,
                                 nest_allocation=a_jk,mu=mu,
                                 marginal = div_calc_marginal)

  # Still consider whether matchCost = FALSE should be an available option.

  matchCost <- TRUE

  ## this was the objective function with old weight method
  if (useOldWeight == TRUE) {

    pdiff <- price - p_model
    sdiff <- shares - share_m
    div_diff <- sum((as.numeric(div_matrix - diversions_m)^2), na.rm = TRUE)

    cost_diff <- sum((cost - cost_cal)^2, na.rm = TRUE)

    if (matchCost == FALSE) {
      objfxn <- c(pdiff,sdiff,div_diff) %*% weight %*% c(pdiff,sdiff,div_diff)
    }
    if (matchCost == TRUE) {
      objfxn <- c(pdiff,sdiff,div_diff,cost_diff) %*% weight %*%
        c(pdiff,sdiff,div_diff,cost_diff)
    }
  }

  ## objective function with new weight method
  if (useOldWeight == FALSE) {

    pdiff <- ((price - p_model)^2) * weight[1]
    sdiff <- ((shares - share_m)^2) * 1000 * weight[2]
    div_diff <- ((as.numeric(div_matrix - diversions_m)^2)*100 * weight[3])
    cost_diff <- ((cost - cost_cal)^2) * weight[4]

    if (matchCost == TRUE) {
      objfxn <- sum(pdiff) + sum(sdiff) + sum(div_diff, na.rm = TRUE) +
        sum(cost_diff, na.rm = TRUE)
    }
  }


  if (returnOutcomes == FALSE) {
    return(objfxn)
  }
  if (returnOutcomes == TRUE) {
    return(list("FOCs" = c(pdiff,sdiff,div_diff,cost_diff),
                "delta_cal" = delta,
                "p_model" = p_model,
                "share_m" = share_m,
                "cost_cal" = cost_cal,
                "diversions_m" = diversions_m) )
  }

}

