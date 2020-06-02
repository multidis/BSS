

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


rhoFractionGaussian <- function(j, alpha) 1/2 * ( (j+1)^(2*alpha + 1) - 2*j^(2*alpha + 1) + (j-1)^(2*alpha + 1) )


calculateK1 <- function(alpha) {
  
  rho_vals <- rhoFractionGaussian(1:1e6, alpha)
  summation <- 0
  
  for (n in 2*(1:50)){
    summation <- summation + factorial(n) * a1Coefficients(n)^2 *(1 + 2*sum(rho_vals^n))
  }
  
  sqrt(summation)
}


calculateK2 <- function(alpha) {
  
  sqrt(2*(1 + 2*sum(rhoFractionGaussian(1:1e6, alpha)^2)))
  
}


calculateK3 <- function(alpha) {
  
  rho_vals <- rhoFractionGaussian(1:1e6, alpha)
  summation <- 0
  
  for (n in 2*(1:50)){
    
    summation <- summation + factorial(n) * a3Coefficients(n)^2 *(1 + 2*sum(rho_vals^n))
  }
  
  sqrt(summation)
}


calculateK4 <- function(alpha) {
  
  sqrt(2*6^2 * (1 + 2*sum(rhoFractionGaussian(1:1e6, alpha)^2)) + 24 * (1 + 2*sum(rhoFractionGaussian(1:1e6, alpha)^4)))
  
}


calculateK <- function(p, alpha) {
  
  helper_functions <- list(calculateK1, calculateK2, calculateK3, calculateK4)
  
  helper_functions[[p]](alpha)
  
}


estimateK <- function(Y, p) {
  
  alpha <- bssAlphaFit(Y, p = 2)
  
  calculateK(p, alpha)

}



estimateAccumulatedVolatilityCI <- function(Y, n, p, method = "nonparametric", kernel = "gamma", confidence_level) {

  p_val = 0.5 + 0.5 * confidence_level
  
  z_a = qnorm(p_val)

  K_p <- estimateK(Y, p)
  
  mean <- estimateAccumulatedVolatility(Y, n, p, method = "nonparametric", kernel = "gamma")
  
  var_term <- estimateAccumulatedVolatility(Y, n, 2*p, method = "nonparametric", kernel = "gamma")
  
  var <- K_p * sqrt(var_term)
  
  cbind(mean - var, mean, mean + var)

}








