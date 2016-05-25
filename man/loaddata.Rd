\name{loaddata}
\alias{loaddata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Load a molecular data file}
\description{
  Load a molecular data file from the Leiden Atomic and Molecular
  Database (LAMDA) database for easy inspection of levels and transitions.
}
\usage{
loaddata(molfile)
}
\arguments{
  \item{molfile}{
    name of the molecular data input file.
  }
}

\value{ Returns a list with components "molecule" (molecule name),
  "amass" (molecular mass in a.m.u.), "levels" (energy levels),
  "transitions" (radiative transitions) and "collisions" (collision rates).}
\author{Michael Olberg <michael.olberg@chalmers.se>}
\seealso{ radex }
\examples{
\dontrun{
result <- loaddata("data/hco+.dat")
print(result)
}
}
\keyword{ radiative }
\keyword{ transfer }
