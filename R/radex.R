#' Run the RADEX radtiative transfer code for given input parameters
#' and return RADEX output as a data.frame.
#'
#' @keywords radiative transfer
#' @param molfile filename of molecular data file, including full path
#' @param freqs output frequency range in GHz, e.g c(50,500)
#' @param Tkin kinetic temperature, default 30 K
#' @param Tbg background temperature, default 2.73 K
#' @param Ncol column density of molecule under study in cm^-2, default 1.0E13
#' @param dv line width in km/s, default 1.0
#' @param collisions a list with names "partner" and "density", where "partner"
#'        is an array of strings designating a collision partner (possible
#'        values "H2", "o-H2", "p-H2", "e", "He", "H", "H+") and density is the
#'        corresponding density of that partner per cm^3.
#' @param verbose if TRUE will produce some diagnostic messages on the console.
#' @return a data frame with results
radex <- function(molfile, freqs, Tkin=30.0, Tbg=2.73, Ncol=1.0e13, dv=1.0,
                  collisions=NA, verbose=FALSE) {
    molfile <- paste(molfile, " ", sep="")
    .Fortran(rdxinput, mol=molfile, f=freqs, tk=Tkin, bg=Tbg,
             ncol=Ncol, dv=dv)
    if (is.list(collisions)) {
        if (all(c("partner","density") %in% names(collisions))) {
            for (i in seq(length(collisions[[1]]))) {
                partner <- collisions$partner[i]
                density <- collisions$density[i]
                .Fortran(rdxdensity, partner, density)
            }
        }
    }
    .Fortran(readdata)

    niter <- .Fortran(rdxrun, niter=integer(1))$niter
    nline <- .Fortran(rdxmoldata, level=integer(1), line=integer(1),
                      coll=integer(1), part=integer(1))$line
    if (verbose) {
        cat("* Radex version        : 30 Nov 2011\n")
        cat("* Geometry             : Uniform sphere\n")
        cat(paste("* Molecular data file  :", molfile, "\n"))
        cat(paste("* T(kin)            [K]:", Tkin, "\n"))
        cat("* Density        [cm-3]:")
        for (i in seq(length(collisions[[1]]))) {
            partner <- collisions$partner[i]
            density <- collisions$density[i]
            if (density > 0.0) {
                cat(paste(" ", partner, ":", density, sep=""))
            }
        }
        cat("\n")
        cat(paste("* T(background)     [K]:", Tbg, "\n"))
        cat(paste("* Column density [cm-2]:", Ncol, "\n"))
        cat(paste("* Line width     [km/s]:", dv, "\n"))
        cat(paste("Calculation finished in", niter, "iterations\n"))
    }
    A <- .Fortran(rdxoutput, up=integer(nline), low=integer(nline),
                  nrg=numeric(nline), ghz=numeric(nline),
                  tx=numeric(nline), tau=numeric(nline), ant=numeric(nline),
                  pu=numeric(nline), pl=numeric(nline),
                  flx=numeric(nline), erg=numeric(nline))
    i <- which(A$ghz >= freqs[1] & A$ghz < freqs[2])
    result <- data.frame(up=A$up[i], low=A$low[i],
                         E.up=A$nrg[i], GHz=A$ghz[i], T.ex=A$tx[i],
                         tau=A$tau[i], T.R=A$ant[i],  pop.up=A$pu[i],
                         pop.low=A$pl[i], Kkms=A$flx[i], ergs=A$erg[i])
    result
}
