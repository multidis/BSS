#' Caclulate the coefficients a_1 in the expression for K_1
#' @param n an integer
#' @keywords internal
a1Coefficients <- function(n) {
  if (n %% 2 == 0) {

    if (n ==2){
      return (1 / sqrt(2*pi))
    } else {
      return( 2 / sqrt(2*pi) * phangorn::dfactorial(n-3) / factorial(n))
    }
  } else {
    return(0)
  }
}

#' Caclulate the coefficients a_3 in the expression for K_3
#' @param n an integer
#' @keywords internal
a3Coefficients <- function(n) {
  if (n %% 2 == 0) {

    if (n == 2){
      return (6 / sqrt(2*pi))}
    if (n == 4){
      return (1 / 2 / sqrt(2*pi))
    } else {
      return( 12 / sqrt(2*pi) * phangorn::dfactorial(n-5) / factorial(n))
    }
  } else {
    return(0)
  }
}

#' Caclulate the correlation of a fractional Gaussian - needed in the calculation of K
#' @param j an integer
#' @param alpha a float, the smoothness parameter of the BSS process
#' @keywords internal
rhoFractionGaussian <- function(j, alpha) 1/2 * ( (j+1)^(2*alpha + 1) - 2*j^(2*alpha + 1) + (j-1)^(2*alpha + 1) )

#' Caclulate K_1 for a BSS process
#' @param alpha a float, the smoothness parameter of the BSS process
#' @keywords internal
calculateK1 <- function(alpha) {

  rho_vals <- rhoFractionGaussian(1:1e6, alpha)
  summation <- 0

  for (n in 2*(1:50)){
    summation <- summation + factorial(n) * a1Coefficients(n)^2 *(1 + 2*sum(rho_vals^n))
  }

  sqrt(summation)
}

#' Caclulate K_2 for a BSS process
#' @param alpha a float, the smoothness parameter of the BSS process
#' @keywords internal
calculateK2 <- function(alpha) {

  sqrt(2*(1 + 2*sum(rhoFractionGaussian(1:1e6, alpha)^2)))

}

#' Caclulate K_3 for a BSS process
#' @param alpha a float, the smoothness parameter of the BSS process
#' @keywords internal
calculateK3 <- function(alpha) {

  rho_vals <- rhoFractionGaussian(1:1e6, alpha)
  summation <- 0

  for (n in 2*(1:50)){

    summation <- summation + factorial(n) * a3Coefficients(n)^2 *(1 + 2*sum(rho_vals^n))
  }

  sqrt(summation)
}

#' Caclulate K_4 for a BSS process
#' @param alpha a float, the smoothness parameter of the BSS process
#' @keywords internal
calculateK4 <- function(alpha) {

  sqrt(2*6^2 * (1 + 2*sum(rhoFractionGaussian(1:1e6, alpha)^2)) + 24 * (1 + 2*sum(rhoFractionGaussian(1:1e6, alpha)^4)))

}

#' Caclulate K_p for a BSS process for a given value of p
#' @param p an integer - the power to use for K_p
#' @param alpha a float, the smoothness parameter of the BSS process
#' @keywords internal
calculateK <- function(p, alpha) {

  helper_functions <- list(calculateK1, calculateK2, calculateK3, calculateK4)

  helper_functions[[p]](alpha)

}

#' Estimate K_p for a BSS process, for a given power p
#' @param p an integer - the power to use for K_p
#' @param Y a vector of observations of the BSS process
#' @keywords internal
estimateK <- function(Y, p) {

  alpha <- bssAlphaFit(Y, p = 2)

  calculateK(p, alpha)

}










