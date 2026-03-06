#' Calibrate lambda for Nash bargain in vertical model with simultaneous timing
#'
#' @param lambda Bargaining power of the buyer/retailer
#' @param price_w Upstream or wholesale prices
#' @param own_down Ownership matrix for downstream firms
#' @param own_up Ownership matrix for upstream firms
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost_r Marginal costs for downstream firm for each product
#' @param cost_w Marginal costs for upstream firm for each product
#' @param price_r Retail prices
#'
#' @returns The value of objective function
#'
#' @details This function can be used to calibrate the bargaining parameter
#' in a Nash bargain which is the upstream market of a vertical supply chain.
#' Assumes logit demand.
#'
#' @examples
#' TO BE ADDED.
#'
#' @export



##################################################################
# Nash Product
##################################################################

## A bargaining calibration function for simple logit demand and calibrates
## lambda. Uses bargain_foc_novert_sim().
## Not sure if need to export both this and the _gnl version.

bargain_vert_sim_calibrate <- function(lambda,price_w,own_down,own_up,
                                       alpha,delta,cost_w,cost_r,price_r){

  out <- BBoptim(par = as.numeric(price_w), fn = bargain_foc_novert_sim,
                 own_down = own_down,
                 own_up = own_up, alpha= alpha, delta = delta,
                 c_W = cost_w, c_R = cost_r, lambda = lambda,
                 p_R = price_r)

  p_W1 <- out$par
  pdiff <- price_w - p_W1

  objfxn <- c(pdiff) %*% diag(length(pdiff)) %*% c(pdiff)
  return(objfxn)
}


