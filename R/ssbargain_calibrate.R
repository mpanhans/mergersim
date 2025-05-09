#' Second score bargaining first-order conditions
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
#' @details This function calculate the first-order conditions from a bargaining
#' model that nests the second score auction.
#'
#' @examples
#' TO BE ADDED.
#'
#' @export



##################################################################
# Second score bargaining first-order conditions
##################################################################
# Add warning if alpha > 0. It should be < 0.
# Would be nice if ownership could be either vector of names OR an
# ownership matrix.
# also add checks for dimensions of inputs

ssbargain_calibrate <- function(param,own,price,shares,cost,weight,
                                   lambda,includeMUI=TRUE){

  J <- length(price)
  alpha <- param[1]
  delta <- param[2:(1+J)]

  x0 <- price
  out <- multiroot(f = ssbargain_foc,start = x0, own = own,
                   alpha= alpha, delta = delta, cost = cost,
                   lambda = lambda, includeMUI = includeMUI)

  price_m <- out$root
  share_m <- (exp(delta + alpha*cost))/(1+sum(exp(delta + alpha*cost)))

  pdiff <- price - price_m
  sdiff <- shares - share_m

  objfxn <- c(pdiff,sdiff) %*% weight %*% c(pdiff,sdiff)
  return(objfxn)
}

