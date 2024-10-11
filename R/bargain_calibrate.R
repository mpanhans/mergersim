#' Nash bargaining calibration
#'
#' @param price Price
#' @param own Ownership matrix
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost Marginal costs for each product
#' @param lambda Bargaining power of the buyer
#' @param includeMUI logical; TO BE ADDED
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
# Nash bargaining calibration
##################################################################
# Add warning if alpha > 0. It should be < 0.
# Would be nice if ownership could be either vector of names OR an
# ownership matrix.
# also add checks for dimensions of inputs

bargain_calibrate <- function(param,own,price,shares,cost,weight,
                               lambda,includeMUI=TRUE){

  J <- length(price)
  alpha <- param[1]
  delta <- param[2:(1+J)]

  x0 <- price
  out <- multiroot(f = bargain_foc, start = x0, own = own,
                   alpha= alpha, delta = delta, cost = cost,
                   lambda = lambda, includeMUI = includeMUI)
  p3 <- out$root
  share3 <- (exp(delta + alpha*p3))/(1+sum(exp(delta + alpha*p3)))

  pdiff <- price - p3
  sdiff <- shares - share3

  objfxn <- c(pdiff,sdiff) %*% weight %*% c(pdiff,sdiff)
  return(objfxn)
}

