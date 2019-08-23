#' Fitting gamma kernel Brownian semistationary processes
#'
#' \code{gammaKernelBSSFit} uses a method of moments to fit the parameters of a gamma kernel Brownian semistationary process
#' to a vector of observations. A least squares estimate of the parameters is obtained
#' by minimising the mean square error between the true gamma kernel autocorrelation function and the
#' empirical ACF of the data, using lags 0,...,H. The number of lags \code{num_lags} used can be adjusted.
#' The volatility process does not need to be specified.
#'
#' @param Y a vector of observations of a BSS process at frequency \code{n}.
#' @param n positive integer indicating the number of observations per unit of time.
#' @param num_lags the number of lags to be used in the regression. The default is to use the first 10 lags.
#'
#' @return The function returns a list containing the parameters \code{alpha} and \code{lambda}, and also the mean square
#' error \code{mse} of the least squares fit. This can be used to compare model fit when trying different kernels.
#' @export
#'
gammaKernelBSSFit <- function(Y, n, num_lags = 10) {

  gammaKernelCorr <- function(theta, h) {

    alpha <- theta[1]
    lambda <- theta[2]

    2^(-alpha + 1/2) / gamma( alpha + 1/2) * (lambda * h)^(alpha + 1/2) * besselK(lambda*h, nu = alpha + 1/2)

  }

  ## loss function which calculates mean of squares of differences between true and observed acf

  Loss <- function(theta, Y) {

    rho_hat <- acf(Y, lag.max = num_lags, demean = FALSE, plot = FALSE)$acf

    h <- (0:num_lags)/n

    true_vals <- sapply(h, gammaKernelCorr, theta = theta)

    true_vals[1] <- 1

    mean( (true_vals - rho_hat)^2)

  }

  optimimum_value <- optim(Loss, par = c(0,1), Y = Y)

  theta <- optimimum_value$par

  mse <- optimimum_value$value

  list(alpha = theta[1], lambda = theta[2], mse = mse)

}




#' Fitting power law kernel Brownian semistationary processes
#'
#' \code{powerKernelBSSFit} uses a method of moments to fit the parameters of a power law kernel Brownian semistationary process
#' to a vector of observations. A least squares estimate of the parameters is obtained
#' by minimising the mean square error between the true power law kernel autocorrelation function (found by numerical intergration)
#' and the empirical ACF of the data, using lags 0,...,H. The number of lags \code{num_lags} used can be adjusted.
#' The volatility process does not need to be specified.
#'
#' @param Y a vector of observations of a BSS process at frequency \code{n}.
#' @param n positive integer indicating the number of observations per unit of time.
#' @param num_lags the number of lags to be used in the regression. The default is to use the first 10 lags.
#'
#' @return The function returns a list containing the parameters \code{alpha} and \code{beta}, and also the mean square
#' error \code{mse} of the least squares fit. This can be used to compare model fit when trying different kernels.
#' @export
#'
powerKernelBSSFit <- function(Y, n, num_lags = 10) {

  powerKernelCorr <- function(theta, h) {

    alpha <- theta[1]

    beta <- theta[2]

    cov_integrand <- function(x) x^alpha * (1 + x)^(beta - alpha) * (x + h)^alpha * (1 + x + h)^(beta - alpha)

    (integrate(cov_integrand, 0, 1)$val + integrate(cov_integrand, 1, Inf)$val)/ beta(1 + 2*alpha, -1 - 2*beta)

  }

  ## loss function which calculates mean of squares of differences between true and observed acf

  Loss <- function(theta, Y) {

    rho_hat <- acf(Y, lag.max = num_lags, demean = FALSE, plot = FALSE)$acf

    h <- (0:num_lags)/n

    true_vals <- sapply(h, powerKernelCorr, theta = theta)

    true_vals[1] <- 1

    mean( (true_vals - rho_hat)^2)

  }

  optimimum_value <- optim(Loss, par = c(-0.1,-2), Y = Y, method = "L-BFGS-B", upper = c(0.5, -0.5))

  theta <- optimimum_value$par

  mse <- optimimum_value$value

  list(alpha = theta[1], beta = theta[2], mse = mse)

}






#' Estimating the smoothness parameter of a Brownian semistationary process
#'
#' \code{bssAlphaFit} uses the 'Change of Frequency' method to estimate the smoothness parameter, \code{alpha},
#' of a BSS process. The COF method needs only minimal assumptions on the parametric form of the kernel,
#' therefore the estimate can be used in any kernel.
#'
#' @param Y a vector of observations of a BSS process at any frequency.
#' @param p the power to be used in the change of frequency method. The default value is p = 2.
#'
#' @return The function returns a single value - an estimate for the smoothness parameter alpha.
#' @export
#'
bssAlphaFit <- function(Y, p = 2) {

  V <- function(Y, p, v) {
    n <- length(Y)
    filter1 <- 0:((n-1)/v)*v + 1
    filter2 <- 1:(n/v)*v

    Y1 <- Y[filter1]
    Y2 <- Y[filter2]
    sum((diff(diff(Y1)))^p) + sum((diff(diff(Y2)))^p)

  }

  COF <- function(X, p) {

    V(Y, p, 2) / sum((diff(diff(Y)))^p)

  }

  h_p <- function(x, p) log2(x) / p - 1/2

  h_p(COF(Y, p), p)

}




