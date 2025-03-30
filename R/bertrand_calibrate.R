#' Bertrand model calibration
#'
#' @param param Price coefficient alpha and mean values delta, parameters to
#' calibrate
#' @param price Price
#' @param own Ownership matrix
#' @param cost Marginal costs for each product
#' @param weight Weighting matrix
#'
#' @returns The first-order conditions
#'
#' @details This function calculate the first-order conditions from a Bertrand
#' price-setting model of competition
#'
#' @examples
#' TO BE ADDED.
#'
#' @export



##################################################################
# Bertrand model calibration
##################################################################
# Add warning if alpha > 0. It should be < 0.
# Would be nice if ownership could be either vector of names OR an
# ownership matrix.
# also add checks for dimensions of inputs
# also add default for weight matrix
# share1 <- assignment needs to be to generalized to account for GNL. As
# written this function can only handle standard logit.
# Note: this function should be re-written to use berry inversion instead of
# searching over delta in param vector.
# This function is only for standard logit demand. For nested logit or
# generalized nested logit, see bertrand_calibrate_gnl().

bertrand_calibrate <- function(param,own,price,shares,cost,weight){

  J <- length(price)
  alpha <- param[1]
  delta <- param[2:(1+J)]

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
