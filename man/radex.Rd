% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radex.R
\name{radex}
\alias{radex}
\title{Run the RADEX radtiative transfer code for given input parameters
and return RADEX output as a data.frame.}
\usage{
radex(molfile, freqs, Tkin = 30, Tbg = 2.73, Ncol = 1e+13, dv = 1,
  collisions = NA, verbose = FALSE)
}
\arguments{
\item{molfile}{filename of molecular data file, including full path}

\item{freqs}{output frequency range in GHz, e.g c(50,500)}

\item{Tkin}{kinetic temperature, default 30 K}

\item{Tbg}{background temperature, default 2.73 K}

\item{Ncol}{column density of molecule under study in cm^-2, default 1.0E13}

\item{dv}{line width in km/s, default 1.0}

\item{collisions}{a list with names "partner" and "density", where "partner"
is an array of strings designating a collision partner (possible
values "H2", "o-H2", "p-H2", "e", "He", "H", "H+") and density is the
corresponding density of that partner per cm^3.}

\item{verbose}{if TRUE will produce some diagnostic messages on the console.}
}
\value{
a data frame with results
}
\description{
Run the RADEX radtiative transfer code for given input parameters
and return RADEX output as a data.frame.
}
\keyword{radiative}
\keyword{transfer}
