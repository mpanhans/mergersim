#' Bertrand model calibration
#'
#' @param param Price coefficient alpha parameter to calibrate
#' @param price Price
#' @param own Ownership matrix
#' @param cost Marginal costs for each product
#' @param weight Weighting matrix
#'
#' @returns Distance between observed values and model predicted values for
#' prices and shares
#'
#' @details This function calculate the first-order conditions from a Bertrand
#' price-setting model of competition. This function is only for standard logit
#' demand. For nested logit or generalized nested logit, see
#' bertrand_calibrate_gnl().
#'
#' @examples
#' TO BE ADDED.
#'
#' @export



##################################################################
# Bertrand model calibration
##################################################################
# Add warning if alpha > 0. It should be < 0.
# Ownership should be accommodated as either vector of names OR an
# ownership matrix.
# add checks for dimensions of inputs
# add default for weight matrix


bertrand_calibrate <- function(param,own,price,shares,cost,weight){

  J <- length(price)
  alpha <- param[1]


  delta0 <- log(shares) - log(1-sum(shares)) - alpha*price

  find_d <- multiroot(f = match_share, start = delta0,
                      price = price, alpha = alpha, shares_obs = shares,
                      nest_allocation = NA, mu = NA)

  delta <- find_d$root

  x0 <- price
  out1 <- BBoptim(f = bertrand_foc, par = x0,
                  own = own, alpha = alpha,
                  delta = delta, cost = cost,
                  sumFOC = TRUE)

  p1 <- out1$par
  share1 <- (exp(delta + alpha*p1))/(1+sum(exp(delta + alpha*p1)))

  pdiff <- price - p1
  sdiff <- shares - share1

  objfxn <- c(pdiff,sdiff) %*% weight %*% c(pdiff,sdiff)
  return(objfxn)
}
