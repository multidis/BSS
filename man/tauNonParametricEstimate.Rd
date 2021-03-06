% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bss_functions.R
\name{tauNonParametricEstimate}
\alias{tauNonParametricEstimate}
\title{Non-parametric estimate of the scale factor}
\usage{
tauNonParametricEstimate(Y)
}
\arguments{
\item{Y}{a vector of observations of a BSS process.}
}
\value{
The function returns the non-parametric estimate for the scale factor.
Note that this will be scaled by the expectation of the square of the volatitity.
}
\description{
Non-parametric estimate of the scale factor
}
