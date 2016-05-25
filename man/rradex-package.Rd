\name{rradex-package}
\alias{rradex-package}
\alias{rradex}
\docType{package}
\title{
 RADEX radiative transfer code.
}
\description{
  An interface to the RADEX radiative transfer code.
}
\details{
\tabular{ll}{
Package: \tab rradex\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2016-05-25\cr
License: \tab GPL (>= 2)\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Michael Olberg

Maintainer: Michael Olberg <michael.olberg@chalmers.se>
}
\references{
  Van der Tak, F.F.S. et al. 2007, A&A 468, 627
}
\keyword{ radiative transfer }
\examples{
  library(rradex)
  l = list(partner=c("H2","e","He"), density=c(1.0e5,1.0e4,3.0e3))
\dontrun{
  result <- radex("hco+.dat", freqs=c(50,500), collisions=l)
  print(result)
}
}
