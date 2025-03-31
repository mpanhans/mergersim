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


match_share <- function(delta,price,alpha,nest_allocation,mu,shares_obs){
  out <- share_calc(price=price, delta=delta, alpha=alpha,
                    nest_allocation=nest_allocation, mu=mu) -
    shares_obs
  return(out)
}

