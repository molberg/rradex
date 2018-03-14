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
al. 2007, A&A 468, 627. The formalism is also explained in the
[RADEX manual](http://www.strw.leidenuniv.nl/~moldata/radex_manual.pdf).

Please note, that the RADEX code has not been put under a specific
license by the original authors. For their downloadable version (on
which this package is built) they state:

> Everyone is free to use the program, provided that publications make a reference to our paper.

So the same will apply to this R interface. However, all code added by
me is GPL-licensed.

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

First make sure that you have all the dependencies to build the
package. After having cloned the repository execute the following
command in its parent directory:

```
    R CMD check rradex
```
This may produce errors due to missing algebra libraries which `rradex` is using
to solve the linear equations. E.g. on Linux typically a

```
    sudo apt-get install libblas-dev liblapack-dev
```
or similar might be needed. This step will produce a subdirectory `rradex.Rcheck`
which will contain a manual for the `rradex` package in pdf format, once the above
step succeeds.

Then to install, simply clone this repository and then in its parent directory

```
   R CMD INSTALL rradex
```
