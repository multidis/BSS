#' Hybrid scheme covariance matrix
#'
#' Generates the covariance matrix used in simulating Brownian semistationary processes by
#' the hybrid scheme.
#'
#' @param kappa number of terms needed for the lower sum in the hybrid scheme.
#' @param n number of observations per unit of time, n = 1/delta.
#' @param alpha smoothness parameter used in the BSS simulation.
#'
#' @return Returns the covariance matrix for the lower sum in the hybrid scheme calculations.
#' The dimensions of the covariance matrix will be (kappa + 1) by (kappa + 1).
#'
#' @export
hybridSchemeCovarianceMatrix <- function(kappa, n, alpha) {
  # create empty matrix
  Sigma <- matrix(0, nrow = kappa + 1, ncol = kappa + 1)
  # fill in top corner = Var(W_i)
  Sigma[1,1] <- 1/n
  # loop over other columns
  for (j in 2:(kappa + 1)) {
    # fill in according to given expressions
    Sigma[1,j] <- ((j-1)^(alpha+1) - (j-2)^(alpha+1))/(alpha+1)/n^(alpha+1)

    Sigma[j,j] <- ((j-1)^(2*alpha+1) - (j-2)^(2*alpha+1))/(2*alpha+1)/n^(2*alpha+1)
    # fill in remaining rows
    if (j < kappa + 1) {
      for (k in (j+1):(kappa+1)) {
        Sigma[j,k] <- 1/(alpha + 1)/n^(2*alpha + 1) *
          ((j - 1)^(alpha + 1) * (k - 1)^alpha *
             hypergeo::hypergeo(-alpha, 1, alpha + 2, (j - 1)/(k - 1) ) -
             (j - 2)^(alpha + 1) * (k - 2)^alpha *
             hypergeo::hypergeo(-alpha, 1, alpha + 2, (j - 2)/(k - 2) ))
      }
    }
  }
  # loop has given an upper triangular (possibly complex) matrix
  # Imaginary part = 0 but display shows complex values (r + 0i) so take real part
  # fill in lower triangle so that S[i,j] = S[j,i]
  Re(Sigma + t(upper.tri(Sigma) * Sigma))
}





exponentiatedOrnsteinUhlenbeck <- function(N, n, T, theta, beta) {
  # initialise values
  v <- numeric(N + n*T + 1)
  # to start in stationary distribution:
  v[1] <- rnorm(1, 0, sqrt(1/(2*theta)))
  # otherwise start from the mean:
  # v[1] <- 0
  for (i in 2:(N + n*T + 1)) {
    v[i] <- v[i-1] - theta * v[i-1] * 1/n + sqrt(1/n) * rnorm(1, 0, 1)
  }
  exp(beta * v)
}






#' Simulation of gamma kernel Brownian semistationary processes
#'
#' \code{gammaKernelBSS} uses the Hybrid scheme to simulate a Brownian semistationary process from the
#' gamma kernel. It simulates a path where the volatility process is independent of the driving Brownian motion of the
#' BSS process.
#'
#' @param N positive integer determining the number of terms in the Riemman sum element of the
#' hybrid scheme calculation. Should be of order at least \code{n}.
#' @param n positive integer indicating the number of observations per unit of time. It represents the fineness or frequency of observations.
#' @param T the time interval to simulate the BSS process over.
#' @param kappa positive integer giving the number of terms to use in the 'lower' sum of the hybrid scheme. Default set to 3.
#' @param alpha the smoothness parameter of the BSS process to simulate.
#' @param lambda the exponent parameter of the BSS process to simulate.
#' @param sigma the volatility process used in the BSS simulation. This should be a vector of length \code{N + n*T + 1}
#' representing the sample path of sigma from -N to nT. By default this is set to by a vector of 1s so that the
#' Gaussian core is simulated.
#'
#' @return The function returns a list of three objects, \code{core} gives the Gaussian core of the process
#' between 0 and T, at intervals of 1/n. \code{bss} gives the BSS sample path on the between 0 and T, at intervals of 1/n,
#' and \code{vol} gives the volatilty process over the same time period.
#' @export
#'
gammaKernelBSS <- function(N, n, T, kappa = 3, alpha, lambda, sigma = rep(1, N + n*T + 1)) {

  ## initialise kernel and discretization parameters:

  # create empty vectors for the 'lower' part of the hybrid scheme sums
  X_lower <- numeric(n*T + 1)
  Y_lower <- numeric(n*T + 1)

  # define gamma kernel
  g <- function(x) x^alpha * exp(-lambda*x)

  # split indices into lower and upper sums
  k_lower <- 1:kappa
  k_upper <- (kappa + 1):N

  # function to calculate the optimal discretization
  b_star <- function(k) ((k^(alpha + 1) - (k-1)^(alpha + 1))/(alpha + 1))^(1/alpha)

  # vector of the L_g(k/n)
  L_g <- exp(-lambda*k_lower/n)

  # vector of g(b*/n) for hybrid scheme
  g_b_star <- c(rep(0, kappa), g(b_star(k_upper)/n))

  ## generate the Brownian increments according to hybrid scheme

  # create the required covariance matrix
  Sigma_W <- hybridSchemeCovarianceMatrix(kappa, n, alpha)

  # sample N + n*T random variables from this multivariate Gaussian
  W <- MASS::mvrnorm(N + n*T, mu = rep(0,kappa + 1), Sigma = Sigma_W)

  ## create the sample hybrid scheme and Riemann sum sample paths
  # split into cases as when kappa = 1, we are dealing with a scalar not a matrix in the first sum
  if (kappa == 1) {
    # loop over each time i/n with i = 0, ..., n*T
    for (i in 1:(n*T + 1)) {
      # calculate X[i] = X(i-1/n) from hybrid scheme
      # sum first kappa terms and remaining N - kappa terms separately

      # generate the Gaussian core element
      X_lower[i] <- sum( L_g  * W[(i + N - kappa):(i+N-1), 2])

      # generate the BSS sample path element
      Y_lower[i] <- sum( L_g  * sigma[(i + N - kappa):(i+N-1)] * W[(i + N - kappa):(i+N-1), 2])
    }

    # add the 'upper' term using the convolution

    # Gaussian core, convolve on with Brownian increments
    X <- X_lower + convolve( g_b_star, rev(W[,1]), type = 'open')[N:(N+n*T)]

    # BSS sample path, convolve with volatility process * Brownian increments
    Y <- Y_lower + convolve( g_b_star, rev(head(sigma,-1) * W[,1]), type = 'open')[N:(N+n*T)]
  } else { # if kappa > 1
    # loop over each time i/n with i = 0, ..., n*T
    for (i in 1:(n*T + 1)) {
      # sum first kappa terms and remaining N - kappa terms separately

      # for the Gaussian core
      X_lower[i] <- sum( L_g  * diag(W[(i + N - kappa):(i+N-1), 2:(kappa + 1)][kappa:1, 1:kappa]))

      # for the BSS sample path
      Y_lower[i] <- sum( L_g  * sigma[(i + N - kappa):(i+N-1)] * diag(W[(i + N - kappa):(i+N-1), 2:(kappa + 1)][kappa:1, 1:kappa]))
    }

    # Gaussian core, convolve on with Brownian increments
    X <- X_lower + convolve( g_b_star, rev(W[,1]), type = 'open')[N:(N+n*T)]

    # BSS sample path, convolve with volatility process * Brownian increments
    Y <- Y_lower + convolve( g_b_star, rev(head(sigma,-1) * W[,1]), type = 'open')[N:(N+n*T)]
  }
  # return Gaussian core, BSS sample path and volatility process for [0,T]
  return list(core = X, bss = Y, vol = tail(sigma, n*T + 1))
}










#' Simulation of power law kernel Brownian semistationary processes
#'
#' \code{powerKernelBSS} uses the Hybrid scheme to simulate a Brownian semistationary process from the
#' power law kernel. It simulates a path where the volatility process is independent of the driving Brownian motion of the
#' BSS process.
#'
#' @param N positive integer determining the number of terms in the Riemman sum element of the
#' hybrid scheme calculation. Should be of order at least \code{n}.
#' @param n positive integer indicating the number of observations per unit of time. It represents the fineness or frequency of observations.
#' @param T the time interval to simulate the BSS process over.
#' @param kappa positive integer giving the number of terms to use in the 'lower' sum of the hybrid scheme. Default set to 3.
#' @param alpha the smoothness parameter of the BSS process to simulate.
#' @param beta the exponent parameter of the BSS process to simulate.
#' @param sigma the volatility process used in the BSS simulation. This should be a vector of length \code{N + n*T + 1}
#' representing the sample path of sigma from -N to nT. By default this is set to by a vector of 1s so that the
#' Gaussian core is simulated.
#'
#' @return The function returns a list of three objects, \code{core} gives the Gaussian core of the process
#' between 0 and T, at intervals of 1/n. \code{bss} gives the BSS sample path on the between 0 and T, at intervals of 1/n,
#' and \code{vol} gives the volatilty process over the same time period.
#' @export
#'
powerKernelBSS <- function(N, n, T, kappa, alpha, beta, sigma = rep(1, N + n*T + 1)) {

  ## initialise kernel and discretization parameters:

  # create empty vectors for the 'lower' part of the hybrid scheme sums
  X_lower <- numeric(n*T + 1)
  Y_lower <- numeric(n*T + 1)

  # define gamma kernel
  g <- function(x) x^alpha * (1 + x)^(beta - alpha)

  # split indices into lower and upper sums
  k_lower <- 1:kappa
  k_upper <- (kappa + 1):N

  # function to calculate the optimal discretization
  b_star <- function(k) ((k^(alpha + 1) - (k-1)^(alpha + 1))/(alpha + 1))^(1/alpha)

  # vector of the L_g(k/n)
  L_g <- (1 + k_lower/n)^(beta - alpha)

  # vector of g(b*/n) for hybrid scheme
  g_b_star <- c(rep(0, kappa), g(b_star(k_upper)/n))

  ## generate the Brownian increments according to hybrid scheme

  # create the required covariance matrix
  Sigma_W <- hybridSchemeCovarianceMatrix(kappa, n, alpha)

  # sample N + n*T random variables from this multivariate Gaussian
  W <- MASS::mvrnorm(N + n*T, mu = rep(0,kappa + 1), Sigma = Sigma_W)

  ## create the sample hybrid scheme and Riemann sum sample paths
  # split into cases as when kappa = 1, we are dealing with a scalar not a matrix in the first sum
  if (kappa == 1) {
    # loop over each time i/n with i = 0, ..., n*T
    for (i in 1:(n*T + 1)) {
      # calculate X[i] = X(i-1/n) from hybrid scheme
      # sum first kappa terms and remaining N - kappa terms separately

      # generate the Gaussian core element
      X_lower[i] <- sum( L_g  * W[(i + N - kappa):(i+N-1), 2])

      # generate the BSS sample path element
      Y_lower[i] <- sum( L_g  * sigma[(i + N - kappa):(i+N-1)] * W[(i + N - kappa):(i+N-1), 2])
    }

    # add the 'upper' term using the convolution

    # Gaussian core, convolve on with Brownian increments
    X <- X_lower + convolve( g_b_star, rev(W[,1]), type = 'open')[N:(N+n*T)]

    # BSS sample path, convolve with volatility process * Brownian increments
    Y <- Y_lower + convolve( g_b_star, rev(head(sigma,-1) * W[,1]), type = 'open')[N:(N+n*T)]
  } else { # if kappa > 1
    # loop over each time i/n with i = 0, ..., n*T
    for (i in 1:(n*T + 1)) {
      # sum first kappa terms and remaining N - kappa terms separately

      # for the Gaussian core
      X_lower[i] <- sum( L_g  * diag(W[(i + N - kappa):(i+N-1), 2:(kappa + 1)][kappa:1, 1:kappa]))

      # for the BSS sample path
      Y_lower[i] <- sum( L_g  * sigma[(i + N - kappa):(i+N-1)] * diag(W[(i + N - kappa):(i+N-1), 2:(kappa + 1)][kappa:1, 1:kappa]))
    }

    # Gaussian core, convolve on with Brownian increments
    X <- X_lower + convolve( g_b_star, rev(W[,1]), type = 'open')[N:(N+n*T)]

    # BSS sample path, convolve with volatility process * Brownian increments
    Y <- Y_lower + convolve( g_b_star, rev(head(sigma,-1) * W[,1]), type = 'open')[N:(N+n*T)]
  }
  # return Gaussian core, BSS sample path and volatility process for [0,T]
  return list(core = X, bss = Y, vol = tail(sigma, n*T + 1))
}











