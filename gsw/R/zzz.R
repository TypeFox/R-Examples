.onLoad <- function(libname, pkgname)
{
    ## To get the data() to work and have the package test the build tests,
    ## we must include
    ##     LazyData: no
    ##     Imports:utils
    ## in DESCRIPTION and
    ##     importFrom(utils, data)
    ## in NAMESPACE.
    saar <- NULL
    data("saar", package=pkgname, envir=environment())
    .C("set_up_gsw_data",
       as.integer(saar$gsw_nx),
       as.integer(saar$gsw_ny),
       as.integer(saar$gsw_nz),
       as.double(saar$longs_ref),
       as.double(saar$lats_ref),
       as.double(saar$p_ref),
       as.double(saar$ndepth_ref),
       as.double(saar$saar_ref),
       as.double(saar$delta_sa_ref))
}

.onUnload <- function(libpath)
{
    .C("clear_gsw_data")
}
