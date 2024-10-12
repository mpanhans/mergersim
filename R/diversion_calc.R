#' Market share calculations
#'
#' Calculates choice probabilities based on logit or generalized nested logit
#' demand parameters
#'
#' @param price Price
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param mu Nesting parameters for each nest
#' @param returnLogsum logical; whether to return the denominator of the choice
#' probabilities (also known as the log-sum term). Defaults to FALSE, in which
#' case the function returns a vector with each product's choice probability.
#'
#' @returns Returns vector of choice probabilities for each good
#'
#' @details This function calculates choice probabilities based on demand
#' parameters for a logit or generalized nested logit demand system
#'
#' @examples
#' TO BE ADDED.
#'
#' @export



##################################################################
# Calculate choice probabilities
##################################################################
#### Add check for whether rows of nest_allocation sum to 1.

diversion_calc <- function(price,alpha,delta,nest_allocation=NA,mu=NA,
                                           marginal = FALSE){

  # Define GNL objects
  a_jk <- nest_allocation
  B <- 1*(a_jk > 0)

  obs_share <- share_calc(price = price, alpha = alpha, delta = delta,
                              nest_allocation = a_jk, mu = mu)

  J <- length(delta)
  div_out <- matrix(0, nrow = J, ncol = J)

  if (marginal == FALSE) {
    for (i in 1:J) {
      delta_cf <- delta
      delta_cf[i] <- -Inf
      cf_sharei <- share_calc(price = price, alpha = alpha, delta = delta_cf,
                              nest_allocation = a_jk, mu = mu)
      div_out[i,] <- (cf_sharei - obs_share)/obs_share[i]
    }

    diag(div_out) <- 0
    return(div_out)
  }

  if (marginal == TRUE) {
    for (i in 1:J) {
      p_cf <- price
      p_cf[i] <- p_cf[i] * 1.01   # increase price by 1%
      cf_sharei <- share_calc(price = p_cf, alpha = alpha, delta = delta,
                              nest_allocation = a_jk, mu = mu)
      div_out[i,] <- (cf_sharei - obs_share)/(obs_share[i] - cf_sharei[i])
    }

    diag(div_out) <- 0
    return(div_out)
  }

}

