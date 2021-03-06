#' Realised power variation
#'
#' \code{realisedPowerVariation} calculates the realised power variation for a BSS process, which is then
#' used as a component of estimating the accumulated volatility process.
#'
#' @param Y a vector of observations of a BSS process.
#' @param p the power variation to be calculated.
#'
#' @return The function returns the sum of the pth powers of the absolute first differences of the BSS process.
#' @export
#'
realisedPowerVariation <- function(Y, p) {

  sum(abs(diff(Y))^p)

}




#' Autocorrelation function for the gamma kernel
#'
#' \code{gammaKernelCorrelation} calculates the value of the gamma kernel autocorrelation function
#' directly using the analytic expression.
#'
#' @param alpha the smoothness parameter, alpha, for the gamma kernel.
#' @param lambda the exponent parameter, lambda, for the gamma kernel.
#' @param h the lag to calculate the autocorrelation at.
#'
#' @return The function returns the autocorrelation for the gamma kernel with
#' parameters \code{alpha} and \code{lambda} at lag \code{h}.
#' @export
#'
gammaKernelCorrelation <- function(alpha, lambda, h) {

  ifelse(h == 0, 1,
         2^(-alpha + 1/2) / gamma( alpha + 1/2) * (lambda * h)^(alpha + 1/2) * besselK(lambda*h, nu = alpha + 1/2))

}


#' Autocorrelation function for the power law kernel
#'
#' \code{powerKernelCorrelation} calculates the value of the power law kernel autocorrelation function
#' directly using numerical integration for the numerator (the covariance term) and the analytic expression
#' for the denominator (variance term).
#'
#' @param alpha the smoothness parameter, alpha, for the power law kernel.
#' @param beta the exponent parameter, beta, for the power law kernel.
#' @param h the lag to calculate the autocorrelation at.
#'
#' @return The function returns the autocorrelation for the power law kernel with
#' parameters \code{alpha} and \code{beta} at lag \code{h}.
#' @export
#'
powerKernelCorrelation <- function(alpha, beta, h) {

  cov_integrand <- function(x) x^alpha * (1 + x)^(beta - alpha) * (x + h)^alpha * (1 + x + h)^(beta - alpha)

  max(integrate(cov_integrand, 0, Inf)$val / beta(1 + 2*alpha, -1 - 2*beta), 1)

}



#' Scale factor for the gamma kernel
#'
#' \code{gammaKernelTau} evaluates the scale factor tau_n for the gamma kernel using the
#' exact expression derived from the covariance function.
#'
#' @param n a positive integer indicating the number of observations per unit of time.
#' @param alpha the smoothness parameter, alpha, for the gamma kernel.
#' @param lambda the exponent parameter, lambda, for the gamma kernel.
#'
#' @return The function returns the scale factor (tau_n) for the gamma kernel with
#' parameters \code{alpha} and \code{lambda}, observed at frequency \code{n} per unit of time.
#' @export
#'
gammaKernelTau <- function(n, alpha, lambda) {

  lambda^(-alpha - 1)*sqrt( gamma(alpha + 1)/gamma(1/2) * (gamma(alpha + 1/2) - 2^(-alpha + 1/2) * (lambda/n)^(alpha + 1/2) * besselK(lambda/n , alpha + 1/2))   )

}




#' Scale factor for the power law kernel
#'
#' \code{powerKernelTau} evaluates the scale factor tau_n for the power law kernel using
#' numerical integration for the covariance term, and exact evaluation for the variance term.
#'
#' @param n a positive integer indicating the number of observations per unit of time.
#' @param alpha the smoothness parameter, alpha, for the power law kernel.
#' @param beta the exponent parameter, beta, for the power law kernel.
#'
#' @return The function returns the scale factor (tau_n) for the power law kernel with
#' parameters \code{alpha} and \code{beta}, observed at frequency \code{n} per unit of time.
#' @export
#'
powerKernelTau <- function(n, alpha, beta) {

  cov_integrand <- function(x) x^alpha * (1 + x)^(beta - alpha) * (x + 1/n)^alpha * (1 + x + 1/n)^(beta - alpha)

  cov_delta <- integrate(cov_integrand, 0, Inf)$val

  sqrt(2*beta(1+2*alpha, -1 -2*beta) - 2*cov_delta)

}




#' Asymptotic scale factor for the gamma kernel
#'
#' @param n a positive integer indicating the number of observations per unit of time.
#' @param alpha the smoothness parameter, alpha, for the gamma kernel.
#'
#' @return The function returns an approximation for the scale factor (tau_n) for the gamma kernel with
#' smoothness parameter \code{alpha}, observed at frequency \code{n} per unit of time, using the asymptotic
#' expression for the scale factor.
#' @export
#'
gammaKernelTauAsymptotic <- function(n, alpha) {

    sqrt(2^(-4*alpha - 1) * gamma(2*alpha + 1) * gamma(1/2 - alpha) / gamma(alpha + 3/2) /n^(2*alpha + 1))

}


#' Non-parametric estimate of the scale factor
#'
#' @param Y a vector of observations of a BSS process.
#'
#' @return The function returns the non-parametric estimate for the scale factor.
#' Note that this will be scaled by the expectation of the square of the volatitity.
#' @export
#'
tauNonParametricEstimate <- function(Y) {

  sqrt(mean(diff(Y)^2))

}
