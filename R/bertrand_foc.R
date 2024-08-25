#' Bertrand first-order conditions
#'
#' @param price Price
#' @param own Ownership matrix
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost Marginal costs for each product
#' @param a_jk For generalized nested logit demand
#' @param B For generalized nested logit demand
#' @param mu Nesting parameters for each nest
#' @param sumFOC logical; whether to return the sum of squares of
#' the first-order conditions. Defaults to FALSE, in which case it returns each
#' product first-order condition as a vector.
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
# Bertrand model first-order conditions
##################################################################
# Add warning if alpha > 0. It should be < 0.
# Would be nice if ownership could be either vector of names OR an
# ownership matrix.
# also add checks for dimensions of inputs

bertrand_foc <- function(price, own, alpha, delta, cost,
                         a_jk=NA, B=NA, mu=NA, sumFOC = FALSE){

  # If no GNL parameters, treat as standard logit. One nest. mu=1.
  J <- length(price)
  if (any(is.na(B))) {
    K <- 1
    B <- matrix(1, ncol = 1, nrow = J)
    a_jk <- B
    mu <- rep(1,K)
  }

  # first create ownership matrices
  #own_fun_down <- function(x) {as.numeric(x == own_down)}
  #own_R <- t(sapply(own_down, own_fun_down) )
  # Since ownership is given by matrix directly:
  own_R <- own

  # then calculate foc
  shares <- share_calc(price = price, alpha = alpha, delta = delta,
                           a_jk = a_jk, B=B, mu = mu)
  m <- price - cost
  dd <- jacobian(share_calc, x = price, delta = delta, alpha = alpha,
                 a_jk=a_jk,B=B,mu=mu)
  omega <- (own_R * t(dd))

  # if all costs available:
  # foc <- omega %*% m + shares

  # To deal with possible missing/NA costs
  exl_na <- which(!is.na(m))
  margin_tilde <- m[exl_na]
  omega_tilde <- omega[exl_na,exl_na]
  shares_tilde <- shares[exl_na]

  foc <- omega_tilde %*% margin_tilde + shares_tilde

  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }

}

