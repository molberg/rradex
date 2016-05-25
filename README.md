# rradex

This is an R-interface to
[RADEX](http://home.strw.leidenuniv.nl/~moldata/radex.html), a
statistical equilibrium radiative transfer code, made available for
public use as part of the Leiden Atomic and Molecular Database
(LAMDA). RADEX is a one-dimensional non-LTE radiative transfer code,
that uses the escape probability formulation assuming an isothermal and
homogeneous medium without large-scale velocity fields. RADEX is
comparable to the LVG method and provides a useful tool in rapidly
analyzing a large set of observational data providing constraints on
physical conditions, such as density and kinetic temperature.

The formalism adopted in RADEX is summarized in Van der Tak, F.F.S. et
al.  2007, A&A 468, 627.

Here is an example of what a short session may look like. Note, that the
various columns of the results table will be immediately available as
vectors in R for further calculations, plotting, etc:

``` r
    library(rradex)
    # 
    # set up collision data
    l = list(partner=c("H2","e","He"), density=c(1.0e5,1.0e4,3.0e3))
    # perform calculation, use defaults for Tkin, Tbg, Ncol
    result <- radex("hco+.dat", freqs=c(50,500), collisions=l)
    print(result)
```
