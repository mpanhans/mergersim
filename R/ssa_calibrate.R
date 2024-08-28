#' Second score auction calibration
#'
#' @param param The price coefficient alpha to be calibrated
#' @param own Ownership matrix
#' @param price Price
#' @param share Share
#' @param cost Marginal costs for each product
#' @param weight Weighting matrix
#'
#' @returns The first-order conditions
#'
#' @details This function calculate the first-order conditions from a second
#' score auction model of competition
#'
#' @examples
#' TO BE ADDED.
#'
#' @export



##################################################################
# Second score auction calibration
##################################################################
# Add warning if alpha > 0. It should be < 0.
# Would be nice if ownership could be either vector of names OR an
# ownership matrix.
# also add checks for dimensions of inputs
# also add default for weight matrix
# This needs to be generalized to GNL

ssa_calibrate <- function(param,own,price,share,cost,weight){
  alpha <- param[1]

  p_hat <- cost + log(1 - own%*%share) / (alpha*own%*%share)

  FOC <- p_hat - price
  FOC <- na.omit(FOC)

  objfxn <- c(FOC) %*% weight %*% c(FOC)

  return(objfxn)
}
