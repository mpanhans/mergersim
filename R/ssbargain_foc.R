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

ssbargain_foc <- function(price,own,alpha,delta,cost,lambda,includeMUI=TRUE){

  wshares <- (exp(delta + alpha*cost))/(1+sum(exp(delta + alpha*cost)))
  m <- price - cost

  own_excl <- own
  diag(own_excl) <- 0

  Emax <- log(1+sum(exp(delta + alpha*cost)))
  Emax_tilde <- log(1+sum(exp(delta + alpha*cost))-own%*%as.matrix(exp(delta + alpha*cost)))
  cond_total_surplus <- (delta+alpha*cost + -digamma(1) + log(1/wshares))
  diff_surplus <- (1/(own%*%as.matrix(wshares)))*(Emax) - (1/(own%*%as.matrix(wshares)))*(Emax_tilde)
  cond_second_surplus <- cond_total_surplus - diff_surplus

  pi_w <- m
  pi_w_tilde <- 0

  if (includeMUI == TRUE) {
    MUI <- (-1/alpha)
    Dpi_r_Dp <- -1
  } else {
    MUI <- 1
    Dpi_r_Dp <- alpha
  }

  pi_r <- MUI*(delta + alpha*price + -digamma(1) + log(1/wshares) )
  pi_r_tilde <- MUI*cond_second_surplus

  Dpi_w_Dp <- 1

  foc <- lambda*(Dpi_r_Dp)*(pi_w-pi_w_tilde) +
    (1-lambda)*(pi_r - pi_r_tilde)*Dpi_w_Dp

  return(foc)
}

