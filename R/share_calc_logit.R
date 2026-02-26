#' Market share calculations
#'
#' Calculates choice probabilities based on logit or generalized nested logit
#' demand parameters
#'
#' @param price Price
#' @param alpha Price coefficient
#' @param delta Mean values
#'
#' @returns Returns vector of choice probabilities for each good
#'
#' @details This function calculates choice probabilities based on demand
#' parameters for a logit demand system
#'
#' @noRd



##################################################################
# Calculate choice probabilities
##################################################################

## simplest version of share calculation for logit demand with outside option
## not exported nor documented, this version is for convenience only.

share_calc_logit <- function(price,delta,alpha){
  out <- (exp(delta + alpha*price))/(1+sum(exp(delta + alpha*price)))
  return(out)
}
