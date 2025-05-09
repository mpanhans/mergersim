% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssa_calibrate.R
\name{ssa_calibrate}
\alias{ssa_calibrate}
\title{Second score auction calibration}
\usage{
ssa_calibrate(param, own, price, share, cost, weight)
}
\arguments{
\item{param}{The price coefficient alpha to be calibrated}

\item{own}{Ownership matrix}

\item{price}{Price}

\item{share}{Share}

\item{cost}{Marginal costs for each product}

\item{weight}{Weighting matrix}
}
\value{
The first-order conditions
}
\description{
Second score auction calibration
}
\details{
This function calculate the first-order conditions from a second
score auction model of competition
}
\examples{
TO BE ADDED.

}
