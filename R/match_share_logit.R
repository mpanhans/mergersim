#' Market share calculations
#'
#' Calculates choice probabilities based on logit or generalized nested logit
#' demand parameters
#'
#' @param price Price
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param shares Observed shares to match
#'
#' @returns Returns vector of difference between predicted shares and observed
#' shares
#'
#' @details This function calculates the difference between model predicted
#' choice probabilities and observed market shares.
#'
#' @examples
#' TO BE ADDED.
#'

#' @export


##################################################################
# Match shares
##################################################################


## match_share but pure logit. Combine w/ other fxn?

match_share_logit <- function(delta,shares,alpha,price){
  output <- shares - (exp(delta + alpha*price))/(1+sum(exp(delta + alpha*price)))
  return(output)
}

