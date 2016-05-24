## Please see ../README_developer.md for the scheme used in adding functions
## here. Generally the functions will be added to Part 4.


## PART 1: document the package and the 'gsw' dataset

#' R implementation of the Thermodynamic Equation Of Seawater - 2010 (TEOS-10)
#'
#' @description
#' Provides an R interface to the TEOS-10 / GSW (Gibbs Sea Water) library,
#' partly for use by the \code{oce} package (see \url{http://dankelley.github.io/oce})
#' and partly for general use. It is assumed that users are familiar with
#' the science and methodology of GSW, and that the package vignette
#' (obtained by typing \code{vignette("gsw")} in an R window) provides 
#' enough orientation to get users started with the \code{gsw} functions.
#'
#' @details
#' \code{gsw} was developed using open-source methodologies, on
#' the GitHub site (\url{https://github.com/TEOS-10/GSW-R}), which
#' is part of a set of sites dedicated to GSW formulations in 
#' various languages.
#'
#' The \code{gsw} system is to link R functions with the C version of
#' the TEOS-10 library.  The R function names are chosen to match
#' those of the Matlab version of GSW, and the function arguments
#' also match with one exception: in \code{gsw}, longitude
#' and latitude are indicated with their full names, whereas in
#' Matlab they are indicated with \code{long} and \code{lat};
#' since R permits abbreviated function arguments, the shortened
#' names can be used in \code{gsw} as well.
#'
#' The documentation for the \code{gsw} functions focuses mainly
#' on the arguments and return values, relying on links to the
#' TEOS-10 webpages for details.
#' 
#' See \url{http://www.teos-10.org/pubs/gsw/html/gsw_contents.html}
#' for a list of the TEOS-10 functions and 
#' \url{http://teos-10.github.io/GSW-R/documentation.html} for a list
#' of the functions implemented in the present package.
#'
#' Each function is tested during the building of the package,
#' which means that results are guaranteed to match those of
#' the equivalent Matlab functions to at least 8 digits.
#'
#' A significant difference from the Matlab case is in the inspection
#' of the dimensions of arguments. The Matlab library has rules 
#' for expanding some arguments to match others. For example,
#' if Practical Salinity is a matrix and pressure is a single value,
#' then that single pressure is used throughout a calculation of 
#' Absolute Salinity. This convenience is only partly mimicked in the
#' present package.  Since the underlying C code works on vectors, 
#' the R functions in \code{gsw} start by transforming the arguments accordingly.
#' This involves using \code{\link{rep}} on each argument to get something
#' with length matching the first argument, and, after the computation
#' is complete, converting the return value into a matrix, if the first
#' argument was a matrix. There are some exceptions to this, however.
#' For example, \code{\link{gsw_SA_from_SP}} and similar functions
#' can handle the case in which the \code{SA} argument is a matrix and
#' \code{longitude} and \code{latitude} are vectors sized to match. 
#' This can be handy with gridded datasets. However, the careful
#' analyst will probably prefer to avoid this and other conveniences,
#' supplying properly-matched arguments from the outset.
#'
#' @docType package
#' @name gsw
NULL

#' Global SA lookup file
#'
#' @description
#' This dataset is not intended for users, but rather for internal use
#' within the \code{gsw} package. The dataset stores the 1.4M lookup
#' table defined in the 8.3M file \code{src/gsw_saar_data.c} in the C
#' library. (The .c file exceeds CRAN limitations on size.)
#'
#' @details
#' The data are designed to replace C elements defined as below
#' in \code{src/gsw_saar_data.c}:
#' \preformatted{
#'     static int	gsw_nx=91, gsw_ny=45, gsw_nz=45;
#'     static double	longs_ref[91];
#'     static double	lats_ref[45];
#'     static double	p_ref[45];
#'     static double	ndepth_ref[4095];
#'     static double	saar_ref[184275];
#'     static double	delta_sa_ref[184275];
#' }
#' 
#' R storage is in a list named \code{saar}, with elements named
#' as in the C code, i.e. \code{gsw_nx} etc.
#'
#' C storage for these variables is allocated as needed,
#' and the data are inserted, when \code{gsw} is launched.
#' Thus, the existing C library code "knows" about the data
#' as local storage, which keeps alterations to the C library to 
#' a minimum.
#'
#' The code used to create the RDA file (using the Fortran data
#' file, version 3.0.3) is given below.
#' \preformatted{
#'     gsw_nx <- 91
#'     gsw_ny <- 45
#'     gsw_nz <- 45
#'     f <- file("~/src/gsw_fortran_v3_03/gsw_data_v3_0.dat", "r")
#'     longs_ref <- scan(f, double(), n=gsw_nx)
#'     lats_ref <- scan(f, double(), n=gsw_ny)
#'     p_ref <- scan(f, double(), n=gsw_nz)
#'     ndepth_ref <- scan(f, double(), n=gsw_nx*gsw_ny)
#'     saar_ref <- scan(f, double(), n=gsw_nx*gsw_ny*gsw_nz)
#'     delta_sa_ref <- scan(f, double(), n=gsw_nx*gsw_ny*gsw_nz)
#'     saar <- list(gsw_nx=gsw_nx, gsw_ny=gsw_ny, gsw_nz=gsw_nz,
#'                  longs_ref=longs_ref, lats_ref=lats_ref, p_ref=p_ref, ndepth_ref=ndepth_ref,
#'                  saar_ref=saar_ref, delta_sa_ref=delta_sa_ref)
#'     save(saar, file="saar.rda")
#'     tools::resaveRdaFiles("saar.rda")
#'     close(f)
#'}
#'
#' @docType data
#' @name saar
NULL


## PART 2: utility functions

#' Reshape list elements to match the shape of the first element.
#'
#' This is mainly used within gsw, to ensure that arguments sent
#' to the C functions are of equal length.  This is a convenience, 
#' for processing data that often have this condition. For example, a
#' CTD profile is likely to have many values for SP, t, and p,
#' but just a single value for each of longitude and latitude.
#' It is important to call argfix() to handle such cases, because
#' otherwise the underlying C code will be looking past the end of
#' the vectors storing longitude and latitude, which can yield odd
#' results or even segmentation faults.
#'
#' @param list A list of elements, typically arguments that will be used in GSW functions.
#' @return A list with all elements of same shape (length or dimension).
argfix <- function(list)
{
    n <- length(list)
    if (n > 0) {
        length1 <- length(list[[1]])
        for (i in 2:n) {
            if (length(list[[i]]) != length1) {
                list[[i]] <- rep(list[[i]], length.out=length1)
            }
        }
        if (is.matrix(list[[1]])) {
            for (i in 2:n) {
                dim(list[[i]]) <- dim(list[[1]])
            }
        }
    }
    list
}



## PART 3: gsw (Gibbs SeaWater) functions, in alphabetical order (ignoring case)

#' Adiabatic lapse rate from Conservative Temperature
#'
#' Note that the unit is K/Pa, i.e. 1e-4 times K/dbar.
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return adiabatic lapse rate (note unconventional unit) [ K/Pa ]
#' @examples
#' gsw_adiabatic_lapse_rate_from_CT(34.7118, 28.7856, 10) # 2.40199646230069e-8
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_adiabatic_lapse_rate_from_CT.html}
gsw_adiabatic_lapse_rate_from_CT <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_adiabatic_lapse_rate_from_CT",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}
                                        
#' Thermal expansion coefficient with respect to Conservative Temperature. (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return thermal expansion coefficient with respect to Conservative Temperature [ 1/K ]
#' @examples
#' gsw_alpha(34.7118, 28.7856, 10) # 3.24480399390879e-3
#' @seealso The salinity analogue to this is \code{\link{gsw_beta}}; other related functions include \code{\link{gsw_beta_const_t_exact}}, \code{\link{gsw_alpha_wrt_t_exact}} and \code{\link{gsw_alpha_on_beta}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_alpha.html}
gsw_alpha <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_alpha",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Thermal expansion coefficient over haline contraction coefficient (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return ratio of thermal expansion coefficient to haline contraction coefficient [ (g/kg)/K ]
#' @examples
#' gsw_alpha_on_beta(34.7118, 28.8099, 10) # 0.452454540612631
#' @seealso This yields the ratio of the return values from \code{\link{gsw_alpha}} and \code{\link{gsw_beta}}, to within computational precision.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_alpha_on_beta.html}
gsw_alpha_on_beta <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_alpha_on_beta",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Thermal expansion coefficient with respect to in-situ temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90)  [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return thermal expansion coefficient with respect to in-situ temperature [ 1/K ]
#' @examples
#' gsw_alpha_wrt_t_exact(34.7118, 28.7856, 10) # 1e-3*0.325601747227247
#' @seealso \code{\link{gsw_alpha}}, \code{\link{gsw_beta}} and \code{\link{gsw_alpha_on_beta}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_alpha_wrt_t_exact.html}
gsw_alpha_wrt_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_alpha_wrt_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Saline contraction coefficient at constant Conservative Temperature. (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return saline contraction coefficient at constant Conservative Temperature [ kg/g ]
#' @examples
#' SA = c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT = c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p =  c(     10,      50,     125,     250,     600,    1000)
#' beta <- gsw_beta(SA,CT,p)
#' @seealso
#' The temperature analogue to this is \code{\link{gsw_alpha}}; other related functions
#' include \code{\link{gsw_alpha_wrt_t_exact}} and \code{\link{gsw_alpha_on_beta}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_beta.html}
gsw_beta <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_beta",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Saline contraction coefficient at constant in-situ temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return saline contraction coefficient at constant in-situ temperature [ kg/g ]
#' @examples
#' gsw_beta_const_t_exact(34.7118, 28.7856, 10) # 7.31120837010429e-4
#' @seealso
#' A related function is \code{\link{gsw_beta}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_beta_const_t_exact.html}
gsw_beta_const_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_beta_const_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Cabbeling coefficient (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return cabbeling coefficient with respect to Conservative Temperature [ 1/(K^2) ]
#' @examples
#' SA = c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT = c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p =  c(     10,      50,     125,     250,     600,    1000)
#' cabbeling <- gsw_cabbeling(SA,CT,p)
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_cabbeling.html}
gsw_cabbeling <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_cabbeling",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}
#' Electrical conductivity from Practical Salinity
#'
#' Note: the return value is not conductivity ratio, but rather
#' conductivity itself, in mS/cm.  To convert to conductivity
#' ratio, divide by gsw_C_from_SP(35, 15, 0).
#' 
#' @param SP Practical Salinity (PSS-78) [ unitless ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return electrical conductivity [ mS/cm ]
#' @examples 
#' gsw_C_from_SP(34.5487, 28.7856, 10) # 56.412599581571186
#' @seealso \code{\link{gsw_SP_from_C}} does the reverse.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_C_from_SP.html}
gsw_C_from_SP <- function(SP, t, p)
{
    l <- argfix(list(SP=SP, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_C_from_SP",
               SP=as.double(l$SP), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
    rval
}

#' Isobaric heat capacity
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return heat capacity [ J/(kg*K) ]
#' @examples 
#' gsw_cp_t_exact(34.7118, 28.7856, 10) # 4002.888003958537
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_cp_t_exact.html}
gsw_cp_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_cp_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}


#' Conservative temperature freezing point
#'
#' Note: as of 2014-12-23, this corresponds to the Matlab function
#' called \code{gsw_t_freezing_poly}. (The confusion arises from a
#' mismatch in release version between the Matlab and C libraries.)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param saturation_fraction saturation fraction of dissolved air in seawater
#' @return Conservative Temperature at freezing of seawater [ deg C ]. That is, the freezing temperature expressed in terms of Conservative Temperature (ITS-90). 
#' @examples 
#' gsw_CT_freezing(34.7118, 10) # -1.899657519404743
#' @seealso \code{\link{gsw_t_freezing}} is the analogue for in-situ temperature.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_CT_freezing.html}
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_CT_freezing_poly.html}
gsw_CT_freezing <- function(SA, p, saturation_fraction=1)
{
    l <- argfix(list(SA=SA, p=p, saturation_fraction=saturation_fraction))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_CT_freezing",
               SA=as.double(l$SA), p=as.double(l$p), saturation_fraction=as.double(l$saturation_fraction),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Conservative Temperature from potential temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param pt potential temperature (ITS-90) [ deg C ]
#' @return Conservative Temperature [ deg C ]
#' @examples 
#' gsw_CT_from_pt(34.7118, 28.7832) # 28.809923015982083
#' @seealso
#' \code{\link{gsw_CT_from_t}} calculates Conservative Temperature from in-situ temperature
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_CT_from_pt.html}
gsw_CT_from_pt <- function(SA, pt)
{
    l <- argfix(list(SA=SA, pt=pt))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_CT_from_pt",
               SA=as.double(l$SA), pt=as.double(l$pt),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Convert from temperature to conservative temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Conservative Temperature [ deg C ]
#' @examples 
#' gsw_CT_from_t(34.7118, 28.7856, 10) # 28.809919826700281
#' @seealso \code{\link{gsw_t_from_CT}} does the reverse
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_CT_from_t.html}
gsw_CT_from_t <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_CT_from_t",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Absolute Salinity Anomaly from Practical Salinity
#' 
#' @param SP Practical Salinity  (PSS-78) [ unitless ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
#' @param latitude latitude in decimal degrees north [ -90 ... +90 ]
#' @return deltaSA Absolute Salinity Anomaly  [ g/kg ]
#' @examples 
#' gsw_deltaSA_from_SP(34.7118, 10, 188, 4) # 0.000167203365230
#' @seealso
#' \code{\link{gsw_SA_from_SP}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_deltaSA_from_SP.html}
gsw_deltaSA_from_SP <- function(SP, p, longitude, latitude)
{
    if (missing(SP)) stop("must supply SP")
    if (missing(p)) stop("must supply p")
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    l <- argfix(list(SP=SP, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_deltaSA_from_SP",
               SP=as.double(l$SP), p=as.double(l$p),
               longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
    rval
}

#' Dynamic enthalpy of seawater (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return dynamic enthalpy [ J/kg ]
#' @examples 
#' gsw_dynamic_enthalpy(34.7118, 28.8099, 10) # 1e3*0.097864649180491
#' @seealso
#' \code{\link{gsw_enthalpy}} and \code{\link{gsw_enthalpy_t_exact}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_enthalpy.html}
gsw_dynamic_enthalpy <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_dynamic_enthalpy",
               SA=as.double(l$SA), t=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Specific enthalpy of seawater (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return specific enthalpy [ J/kg ]
#' @examples 
#' gsw_enthalpy(34.7118, 28.8099, 10) # 1.1510318130700132e5
#' @seealso
#' \code{\link{gsw_dynamic_enthalpy}} and \code{\link{gsw_enthalpy_t_exact}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_enthalpy.html}
gsw_enthalpy <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_enthalpy",
               SA=as.double(l$SA), t=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Specific enthalpy of seawater
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90)  [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return specific enthalpy [ J/kg ]
#' @examples 
#' gsw_enthalpy_t_exact(34.7118, 28.7856, 10) # 1.151032604783763e5
#' @seealso
#' \code{\link{gsw_enthalpy}} and \code{\link{gsw_dynamic_enthalpy}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_enthalpy_t_exact.html}
gsw_enthalpy_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_enthalpy_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Specific entropy as a function of in-situ temperature and pressure
#'
#' Calculates specific entropy given Absolute Salinity, in-situ
#' temperature and pressure.
#'
#' The related function gsw_entropy_from_CT() is not provided
#' in the C library, although it is available in the (later-
#' versioned) Matlab library.
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return specific entropy [ J/(kg*K) ]
#' @examples
#' gsw_entropy_from_t(34.7118, 28.7856, 10) # 400.3894252787245
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_entropy_from_t.html}
gsw_entropy_from_t <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_entropy_from_t",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Gravitational acceleration
#' 
#' @param latitude latitude in decimal degress north [ -90 ... +90 ]
#' @param p sea pressure [ dbar ]
#' @return gravitational acceleration [ m/s^2 ]
#' @examples
#' gsw_grav(c(-90, -60), 0) # 9.832186205884799, 9.819178859991149
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_grav.html}
gsw_grav <- function(latitude, p=0)
{
    l <- argfix(list(latitude=latitude, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_grav",
               latitude=as.double(l$latitude), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(latitude))
        dim(rval) <- dim(latitude)
    rval
}

#' Specific internal energy of seawater (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return specific internal energy [ J/kg ]
#' @examples
#' gsw_internal_energy(34.7118, 28.7856, 10) # 1.148091577452400e5
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_internal_energy.html}
gsw_internal_energy <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_internal_energy",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Ratio of vert. gradient of pot. density to vert grad of locally-referenced pot density
#'
#' Note that the C library had to be patched to get this working; a new
#' version of the library will address the bug directly.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @param p_ref reference pressure [ dbar ]
#' @return list containing IPV_vs_fNsquared_ratio [ unitless ] and mid-point pressure p_mid [ dbar ]
#' @examples 
#' SA <- c(34.7118, 34.8915)
#' CT <- c(28.8099, 28.4392)
#' p <-  c(     10,      50)
#' p_ref <- 0
#' r <- gsw_IPV_vs_fNsquared_ratio(SA, CT, p, p_ref)
#' r$IPV_vs_fNsquared_ratio # 0.999745283730840
#' r$p_mid                  # 30
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_IPV_vs_fNsquared_ratio.html}
gsw_IPV_vs_fNsquared_ratio <- function(SA, CT, p, p_ref=0)
{
    l <- argfix(list(SA=SA, CT=CT, p=p, p_ref=p_ref))
    n <- length(l[[1]])
    ## note: only use p_ref[1] since the C-library code says it must be constant
    r <- .C("wrap_gsw_IPV_vs_fNsquared_ratio",
            SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p), p_ref=as.double(l$p_ref[1]),
            n=as.integer(n),
            ratio=double(n-1), p_mid=double(n-1), NAOK=TRUE, package="gsw")
    if (is.matrix(SA))
        stop("gsw_IPV_vs_fNsquared_ratio() cannot handle matrix SA")
    list(IPV_vs_fNsquared_ratio=r$ratio, p_mid=r$p_mid)
}

#' Isentropic compressibility of seawater
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return isentropic compressibility [ 1/Pa ] (not 1/dbar)
#' @examples
#' gsw_kappa(34.7118, 28.7856, 10) # 4.11346577902628e-10
#' @seealso \code{\link{gsw_kappa_t_exact}} is an analogue in terms of in-situ temperature
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_kappa.html}
gsw_kappa <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_kappa",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Isentropic compressibility of seawater
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return isentropic compressibility [ 1/Pa ] (not 1/dbar)
#' @examples
#' gsw_kappa_t_exact(34.7118, 28.7856, 10) # 4.11245799180373e-10
#' @seealso \code{\link{gsw_kappa}} is an analogue in terms of Conservative Temperature
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_kappa_t_exact.html}
gsw_kappa_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_kappa_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Latent heat of evaporation
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return latent heat of evaporation [ J/kg ]
#' @examples
#' gsw_latentheat_evap_t(34.7118, 28.7856) # 2.429947107462561e6
#' @seealso \code{\link{gsw_latentheat_evap_t}} is an analouge in terms of in-situ temperature. For melting, see \code{\link{gsw_latentheat_melting}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_latentheat_evap_CT.html}
gsw_latentheat_evap_CT <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_latentheat_evap_CT",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Latent heat of evaporation
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @return latent heat of evaporation [ J/kg ]
#' @examples
#' gsw_latentheat_evap_t(34.7118, 28.7856) # 2.429882982734836e6
#' @seealso \code{\link{gsw_latentheat_evap_CT}} is an analogue in terms of Conservative Temperature. For melting, see \code{\link{gsw_latentheat_melting}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_latentheat_evap_t.html}
gsw_latentheat_evap_t <- function(SA, t)
{
    l <- argfix(list(SA=SA, t=t))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_latentheat_evap_t",
               SA=as.double(l$SA), t=as.double(l$t),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Latent heat of melting
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @return latent heat of freezing [ J/kg ]
#' @examples
#' gsw_latentheat_melting(34.7118, 10) # 3.299495187300804e5
#' @seealso \code{\link{gsw_latentheat_evap_CT}} and \code{\link{gsw_latentheat_evap_t}} are analogues for evaporation.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_latentheat_melting.html}
gsw_latentheat_melting <- function(SA, p)
{
    l <- argfix(list(SA=SA, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_latentheat_melting",
               SA=as.double(l$SA), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Calculate Brunt Vaisala Frequency squared
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return list containing N2 [ 1/s^ ] and mid-point pressure p_mid [ dbar ]
#' @examples 
#' SA <- c(34.7118, 34.8915)
#' CT <- c(28.8099, 28.4392)
#' p <- c(      10,      50)
#' latitude <- 4
#' gsw_Nsquared(SA, CT, p, latitude)$N2 # 6.0847042791371e-5
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html}
gsw_Nsquared <- function(SA, CT, p, latitude=0)
{
    l <- argfix(list(SA=SA, CT=CT, p=p, latitude=latitude))
    n <- length(l[[1]])
    r <- .C("wrap_gsw_Nsquared",
            SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p), latitude=as.double(l$latitude),
            n=as.integer(n), n2=double(n-1), p_mid=double(n-1), NAOK=TRUE, package="gsw")
    if (is.matrix(SA))
        stop("gsw_Nsquared() cannot handle matrix SA")
    list(N2=r$n2, p_mid=r$p_mid)
}

#' Pressure from z
#' 
#' @param z height, zero at surface (but note last 2 args) and positive upwards [ m ]
#' @param latitude latitude in decimal degrees north [ -90 ... +90 ]
#' @param geo_strf_dyn_height dynamic height anomaly [ m^2/s^2 ]
#' @param sea_surface_geopotential geopotential at zero sea pressure [ m^2/s^2 ]
#' @return sea pressure [ dbar ]
#' @examples
#' gsw_p_from_z(-10, 4) # 10.05572704136
#' @seealso
#' This is (almost) the reverse of \code{\link{gsw_z_from_p}}, apart from the last two arguments.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_p_from_z.html}
gsw_p_from_z <- function(z, latitude, geo_strf_dyn_height=0, sea_surface_geopotential=0)
{
    l <- argfix(list(z=z, latitude=latitude,
                     geo_strf_dyn_height=geo_strf_dyn_height,
                     sea_surface_geopotential=sea_surface_geopotential))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_p_from_z",
               z=as.double(l$z), latitude=as.double(l$latitude),
               geo_strf_dyn_height=as.double(l$geo_strf_dyn_height),
               sea_surface_geopotential=as.double(l$sea_surface_geopotential),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(z))
        dim(rval) <- dim(z)
    rval
}

#' Potential density
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @param p_ref reference pressure [ dbar ]
#' @return potential density [ kg/m^3 ]
#' @examples
#' gsw_pot_rho_t_exact(34.7118, 28.7856, 10, 0) # 1021.798145811089
#' @seealso
#' \code{\link{gsw_rho}} and \code{\link{gsw_rho_t_exact}} compute density; \code{\link{gsw_sigma0}} and related functions compute potential density at particular pressures.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_pot_rho_t_exact.html}
gsw_pot_rho_t_exact <- function(SA, t, p, p_ref)
{
    l <- argfix(list(SA=SA, t=t, p=p, p_ref=p_ref))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_pot_rho_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p), pref=as.double(l$p_ref),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential temperature referenced to the surface
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return potential temperature [ deg C ]
#' @examples
#' gsw_pt0_from_t(34.7118, 28.7856, 10) # 28.783196819670632
#' @seealso \code{\link{gsw_pt_from_CT}} and \code{\link{gsw_pt_from_t}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_pt0_from_t.html}
gsw_pt0_from_t <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_pt0_from_t",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential temperature from Conservative Temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential temperature [ deg C ]
#' @examples
#' gsw_pt_from_CT(34.7118, 28.8099) # 28.783177048624573 
#' @seealso \code{\link{gsw_pt0_from_t}} for the surface case and and \code{\link{gsw_pt_from_t}} for the general case.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_pt_from_CT.html}
gsw_pt_from_CT <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_pt_from_CT",
               SA=as.double(l$SA), t=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential temperature from in-situ temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90)  [ deg C ]
#' @param p sea pressure [ dbar ]
#' @param p_ref reference pressure [ dbar ]
#' @return potential temperature [ deg C ]
#' @examples
#' gsw_pt_from_t(34.7118, 28.7856, 10, 0) # 28.783196819670632
#' @seealso \code{\link{gsw_pt_from_CT}} is the analogue for Conservative Temperature
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_pt_from_t.html}
gsw_pt_from_t <- function(SA, t, p, p_ref=0)
{
    l <- argfix(list(SA=SA, t=t, p=p, p_ref=p_ref))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_pt_from_t",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p), p_ref=as.double(l$p_ref),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' In-situ density (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return in-situ density [ kg/m^3 ]
#' @seealso \code{\link{gsw_rho_t_exact}} is similar to this, but using in-situ temperature. SA and CT may be computed from UNESCO quantities using \code{\link{gsw_SA_from_SP}} and \code{\link{gsw_CT_from_t}}. For potential density anomalies, use \code{\link{gsw_sigma0}} and related functions.
#' @examples
#' gsw_rho(34.7118, 28.8099, 10) # 1021.8404465661
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_rho.html}
gsw_rho <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_rho",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' SA, CT and p partial derivatives of density (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return list containing drho_dSA [ kg^2/(g m^3) ], drho_dCT [ kg/(K m^3) ] and drho_dp [ kg/(Pa m^3) ]
#' @examples
#' gsw_rho_first_derivatives(34.7118, 28.8099, 10) #' # 0.73321 -0.33174 4.20305e-7
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_rho_first_derivatives.html}
gsw_rho_first_derivatives <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_rho_first_derivatives",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, drho_dSA=double(n), drho_dCT=double(n), drho_dp=double(n),
               NAOK=TRUE, package="gsw")
    if (is.matrix(SA))
        stop("gsw_rho_first_derivatives() cannot handle matrix SA")
    list(drho_dSA=rval$drho_dSA, drho_dCT=rval$drho_dCT, drho_dp=rval$drho_dp)
}

#' In-situ density
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return in-situ density [ kg/m^3 ]
#' @examples
#' gsw_rho_t_exact(34.7118, 28.7856, 10) # 1021.840173185531
#' @seealso \code{\link{gsw_rho}} is similar but uses SA and CT; SA may be computed from UNESCO quantities using \code{\link{gsw_SA_from_SP}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_rho_t_exact.html}
gsw_rho_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_rho_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Convert from density to absolute salinity
#'
#' @param rho seawater density [ kg/m^3 ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Absolute Salinity [ g/kg ]
#' @examples
#' gsw_SA_from_rho(1021.8482, 28.7856, 10) # 34.711382887931144
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_rho.html}
gsw_SA_from_rho <- function(rho, CT, p)
{
    l <- argfix(list(rho=rho, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SA_from_rho",
               SA=as.double(l$rho), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(rho))
        dim(rval) <- dim(rho)
    rval
}

#' Convert from practical salinity to absolute salinity
#'
#' Calculate Absolute Salinity from Practical Salinity, pressure,
#' longitude, and latitude.
#'
#' If SP is a matrix and if its dimensions correspond to the
#' lengths of longitude and latitude, then the latter are
#' converted to analogous matrices with \code{\link{expand.grid}}.
#' 
#' @param SP Practical Salinity (PSS-78) [ unitless ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 to 360 or -180 to 180]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return Absolute Salinity [ g/kg ]
#' @examples
#' gsw_SA_from_SP(34.5487, 10, 188, 4) # 34.711778344814114 
#' @seealso \code{\link{gsw_SP_from_SA}} does the reverse.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_SP.html}
gsw_SA_from_SP <- function(SP, p, longitude, latitude)
{
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    ## check for special case that SP is a matrix defined on lon and lat
    if (is.matrix(SP)) {
        dim <- dim(SP)
        if (length(longitude) == dim[1] && length(latitude) == dim[2]) {
            ll <- expand.grid(longitude=as.vector(longitude), latitude=as.vector(latitude))
            longitude <- ll$longitude
            latitude <- ll$latitude
        }
    }
    l <- argfix(list(SP=SP, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SA_from_SP",
               SP=as.double(l$SP), p=as.double(l$p),
               longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
    rval
}

#' Absolute Salinity from Preformed Salinity
#'
#' Calculate Absolute Salinity from Preformed Salinity, pressure,
#' longitude, and latitude.
#'
#' If Sstar is a matrix and if its dimensions correspond to the
#' lengths of longitude and latitude, then the latter are
#' converted to analogous matrices with \code{\link{expand.grid}}.
#' 
#' @param Sstar Preformed Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 to 360 or -180 to 180]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return Absolute Salinity [ g/kg ]
#' @examples
#' gsw_SA_from_Sstar(34.7115, 10, 188, 4) # 34.711724663585905
#' @seealso \code{\link{gsw_Sstar_from_SA}} does the reverse.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_Sstar.html}
gsw_SA_from_Sstar <- function(Sstar, p, longitude, latitude)
{
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    ## check for special case that Sstar is a matrix defined on lon and lat
    if (is.matrix(Sstar)) {
        dim <- dim(Sstar)
        if (length(longitude) == dim[1] && length(latitude) == dim[2]) {
            ll <- expand.grid(longitude=as.vector(longitude), latitude=as.vector(latitude))
            longitude <- ll$longitude
            latitude <- ll$latitude
        }
    }
    l <- argfix(list(Sstar=Sstar, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SA_from_Sstar",
               Sstar=as.double(l$Sstar), p=as.double(l$p),
               longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(Sstar))
        dim(rval) <- dim(Sstar)
    rval
}

#' Potential density anomaly referenced to 0 dbar
#'
#' This uses the 48-term density equation, and returns
#' potential density referenced to a pressure of 0 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly [ kg/m^3 ]
#' @examples
#' gsw_sigma0(34.7118, 28.8099) # 21.798411276610750
#' @seealso Use \code{\link{gsw_sigma1}} for 1000 dbar pressure, \code{\link{gsw_sigma2}} for 2000 dbar, \code{\link{gsw_sigma3}} for 3000 dbar, or \code{\link{gsw_sigma4}} for 4000 dbar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma0.html}
gsw_sigma0 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma0",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential density anomaly referenced to 1000 dbar
#'
#' This uses the 48-term density equation, and returns
#' potential density referenced to a pressure of 1000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly [ kg/m^3 ]
#' @examples
#' gsw_sigma1(34.7118, 28.8099) # 25.955891533636986
#' @seealso Use \code{\link{gsw_sigma0}} for 0 dbar pressure, \code{\link{gsw_sigma2}} for 2000 dbar, \code{\link{gsw_sigma3}} for 3000 dbar, or \code{\link{gsw_sigma4}} for 4000 dbar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma1.html}
gsw_sigma1 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma1",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential density anomaly referenced to 2000 dbar
#'
#' This uses the 48-term density equation, and returns
#' potential density referenced to a pressure of 2000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly [ kg/m^3 ]
#' @examples
#' gsw_sigma2(34.7118, 28.8099) # 30.022796416066058
#' @seealso Use \code{\link{gsw_sigma0}} for 0 dbar pressure, \code{\link{gsw_sigma1}} for 1000 dbar, \code{\link{gsw_sigma3}} for 3000 dbar, or \code{\link{gsw_sigma4}} for 4000 dbar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma2.html}
gsw_sigma2 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma2",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential density anomaly referenced to 3000 dbar
#'
#' This uses the 48-term density equation, and returns
#' potential density referenced to a pressure of 3000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly with reference pressure 3000 dbar [ kg/m^3 ]
#' @examples
#' gsw_sigma3(34.7118, 28.8099) # 34.002600253012133
#' @seealso Use \code{\link{gsw_sigma0}} for 0 dbar pressure, \code{\link{gsw_sigma1}} for 1000 dbar, \code{\link{gsw_sigma2}} for 2000 dbar, or \code{\link{gsw_sigma4}} for 4000 dbar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma3.html}
gsw_sigma3 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma3",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential density anomaly referenced to 4000 dbar
#'
#' This uses the 48-term density equation, and returns
#' potential density referenced to a pressure of 4000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly with reference pressure 4000 dbar [ kg/m^3 ]
#' @examples
#' gsw_sigma4(34.7118, 28.8099) # 37.898467323406976
#' @seealso Use \code{\link{gsw_sigma0}} for 0 dbar pressure, \code{\link{gsw_sigma1}} for 1000 dbar, \code{\link{gsw_sigma2}} for 2000 dbar, or \code{\link{gsw_sigma3}} for 3000 dbar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma4.html}
gsw_sigma4 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma4",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Sound speed with 48-term density
#'
#' This uses the 48-term density equation.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return sound speed [ m/s ]
#' @examples
#' gsw_sound_speed(34.7118, 28.7856, 10) # 1542.420534932182
#' @seealso \code{\link{gsw_sound_speed_t_exact}} for a precise formula using in-situ temperature
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sound_speed.html}
gsw_sound_speed<- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sound_speed",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Sound speed
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return sound speed [ m/s ]
#' @examples
#' gsw_sound_speed_t_exact(34.7118, 28.7856, 10) # 1542.420534932182
#' @seealso \code{\link{gsw_sound_speed}} for an approximate formula using CT
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sound_speed_t_exact.html}
gsw_sound_speed_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sound_speed_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Specific volume
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Specific volume (1/density)
#' @examples 
#' gsw_specvol(34.7118, 28.8099, 10) # 9.78626363206202e-4
#' @seealso With in-situ temperature, use \code{\link{gsw_specvol_t_exact}}; \code{\link{gsw_specvol_anom}} gives specific volume anomaly.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_specvol.html}
gsw_specvol  <- function(SA, CT, p)
{
    1 / gsw_rho(SA, CT, p)
}

#' Specific volume anomaly
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Specific volume anomaly [ m^3/kg ]
#' @examples 
#' gsw_specvol_anom(34.7118, 28.8099, 10) # 6.01005694856401e-6
#' @seealso Specific volume itself is given by \code{\link{gsw_specvol}} and \code{\link{gsw_specvol_t_exact}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_specvol_anom.html}
gsw_specvol_anom  <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_specvol_anom",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Specific volume
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90)  [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Specific volume [ m^3/kg ]
#' @examples 
#' gsw_specvol_t_exact(34.7118, 28.7856, 10) # 9.78626625025472e-4
#' @seealso With Conservative Temperature, use \code{\link{gsw_specvol}}; \code{\link{gsw_specvol_anom}} gives specific volume anomaly.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_specvol_t_exact.html}
gsw_specvol_t_exact  <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_specvol_t_exact",
               SA=as.double(l$SA), CT=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Convert from electrical conductivity to practical salinity
#' 
#' @param C conductivity [ mS/cm ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Practical Salinity (PSS-78) [ unitless ]
#' @examples 
#' gsw_SP_from_C(34.5487, 28.7856, 10) # 20.009869599086951
#' @seealso \code{\link{gsw_C_from_SP}} does the reverse.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_C.html}
gsw_SP_from_C <- function(C, t, p)
{
    l <- argfix(list(C=C, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SP_from_C",
               C=as.double(l$C), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(C))
        dim(rval) <- dim(C)
    rval
}

#' Convert from Absolute Salinity to Practical Salinity
#'
#' Calculate Practical Salinity from Absolute Salinity, pressure,
#' longitude, and latitude.
#'
#' If SP is a matrix and if its dimensions correspond to the
#' lengths of longitude and latitude, then the latter are
#' converted to analogous matrices with \code{\link{expand.grid}}.
#'
#' Note: unlike the corresponding Matlab function, this does not
#' return a flag indicating whether the location is in the ocean.
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 to 360 or -180 to 180]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return Practical Salinity (PSS-78) [ unitless ]
#' @examples 
#' gsw_SP_from_SA(34.7118, 10, 188, 4) # 34.548721553448317
#' @seealso \code{\link{gsw_SA_from_SP}} does the reverse, while \code{\link{gsw_SP_from_SK}}, \code{\link{gsw_SP_from_SR}} and \code{\link{gsw_SP_from_Sstar}} are similar to this.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_SA.html}
gsw_SP_from_SA <- function(SA, p, longitude, latitude)
{
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    ## check for special case that SP is a matrix defined on lon and lat
    if (is.matrix(SA)) {
        dim <- dim(SA)
        if (length(longitude) == dim[1] && length(latitude) == dim[2]) {
            ll <- expand.grid(longitude=as.vector(longitude), latitude=as.vector(latitude))
            longitude <- ll$longitude
            latitude <- ll$latitude
        }
    }
    l <- argfix(list(SA=SA, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SP_from_SA",
               SA=as.double(l$SA), p=as.double(l$p), longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, SP=double(n), NAOK=TRUE, package="gsw")$SP
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Calculate Practical Salinity from Knudsen Salinity
#'
#' @param SK Knudsen Salinity [ parts per thousand, ppt ]
#' @return Practical Salinity (PSS-78) [ unitless ]
#' @examples 
#' gsw_SP_from_SK(34.5487) # 34.548721553448317
#' @seealso \code{\link{gsw_SP_from_SA}}, \code{\link{gsw_SP_from_SR}} and \code{\link{gsw_SP_from_Sstar}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_SK.html}
gsw_SP_from_SK <- function(SK)
{
    if (missing(SK)) stop("must supply SK")
    n <- length(SK)
    rval <- .C("wrap_gsw_SP_from_SK",
               SA=as.double(SK), n=as.integer(n), SP=double(n), NAOK=TRUE, package="gsw")$SP
    if (is.matrix(SK))
        dim(rval) <- dim(SK)
    rval
}

#' Calculate Practical Salinity from Reference Salinity
#'
#' @param SR Reference Salinity [ g/kg ]
#' @return Practical Salinity (PSS-78) [ unitless ]
#' @examples 
#' gsw_SP_from_SR(34.5487) # 34.386552667080714
#' @seealso The reverse is \code{\link{gsw_SR_from_SP}}; also related are \code{\link{gsw_SP_from_SA}}, \code{\link{gsw_SP_from_SK}} and \code{\link{gsw_SP_from_Sstar}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_SR.html}
gsw_SP_from_SR <- function(SR)
{
    if (missing(SR)) stop("must supply SR")
    n <- length(SR)
    rval <- .C("wrap_gsw_SP_from_SR",
               SA=as.double(SR), n=as.integer(n), SP=double(n), NAOK=TRUE, package="gsw")$SP
    if (is.matrix(SR))
        dim(rval) <- dim(SR)
    rval
}

#' Practical Salinity from Preformed Salinity
#' 
#' @param Sstar Preformed Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
#' @param latitude latitude in decimal degrees north [ -90 ... +90 ]
#' @return Practical Salinity (PSS-78) [ unitless ]
#' @examples 
#' gsw_SP_from_Sstar(34.7115, 10, 188, 4) # 34.548646570969929
#' @seealso
#' \code{\link{gsw_Sstar_from_SP}} does the reverse; \code{\link{gsw_SA_from_Sstar}} is similar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_Sstar.html}
gsw_SP_from_Sstar <- function(Sstar, p, longitude, latitude)
{
    if (missing(Sstar)) stop("must supply Sstar")
    l <- argfix(list(Sstar=Sstar, p=p, longitude=longitude, latitude=latitude))
    if (is.null(l$p)) stop("must supply p")
    if (is.null(l$longitude)) stop("must supply longitude")
    if (is.null(l$latitude)) stop("must supply latitude")
    n <- length(Sstar)
    rval <- .C("wrap_gsw_SP_from_Sstar",
               Sstar=as.double(l$Sstar), p=as.double(l$p), longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(Sstar))
        dim(rval) <- dim(Sstar)
    rval
}

#' Calculate Reference Salinity from Practical Salinity
#'
#' @param SP Practical Salinity (PSS-78) [ unitless ]
#' @return Reference Salinity [ g/kg ]
#' @examples 
#' gsw_SR_from_SP(34.5487) # 34.711611927085727
#' @seealso The reverse is \code{\link{gsw_SP_from_SR}}; also related are \code{\link{gsw_SP_from_SA}}, \code{\link{gsw_SP_from_SK}} and \code{\link{gsw_SP_from_Sstar}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SR_from_SP.html}
gsw_SR_from_SP <- function(SP)
{
    if (missing(SP)) stop("must supply SP")
    n <- length(SP)
    rval <- .C("wrap_gsw_SR_from_SP",
               SP=as.double(SP), n=as.integer(n), SR=double(n), NAOK=TRUE, package="gsw")$SR
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
    rval
}

#' Convert from Absolute Salinity to Preformed Salinity
#'
#' Calculate Preformed Salinity from Absolute Salinity, pressure,
#' longitude, and latitude.
#'
#' If SA is a matrix and if its dimensions correspond to the
#' lengths of longitude and latitude, then the latter are
#' converted to analogous matrices with \code{\link{expand.grid}}.
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 to 360 or -180 to 180]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return Preformed Salinity [ g/kg ]
#' @examples
#' gsw_Sstar_from_SA(34.7118, 10, 188, 4) # 34.711575335926490
#' @seealso \code{\link{gsw_SA_from_Sstar}} does the reverse; \code{\link{gsw_Sstar_from_SP}} is similar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_Sstar_from_SA.html}
gsw_Sstar_from_SA <- function(SA, p, longitude, latitude)
{
    if (missing(SA)) stop("must supply SA")
    if (missing(p)) stop("must supply p")
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    ## check for special case that SP is a matrix defined on lon and lat
    if (is.matrix(SA)) {
        dim <- dim(SA)
        if (length(longitude) == dim[1] && length(latitude) == dim[2]) {
            ll <- expand.grid(longitude=as.vector(longitude), latitude=as.vector(latitude))
            longitude <- ll$longitude
            latitude <- ll$latitude
        }
    }
    l <- argfix(list(SA=SA, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_Sstar_from_SA",
               SA=as.double(l$SA), p=as.double(l$p), longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Convert from Practical Salinity to Preformed Salinity
#'
#' Calculate Preformed Salinity from Practical Salinity, pressure,
#' longitude, and latitude.
#'
#' If SP is a matrix and if its dimensions correspond to the
#' lengths of longitude and latitude, then the latter are
#' converted to analogous matrices with \code{\link{expand.grid}}.
#' 
#' @param SP Practical Salinity (PSS-78) [ unitless ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 to 360 or -180 to 180]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return Preformed Salinity [ g/kg ]
#' @examples
#' gsw_Sstar_from_SP(34.5487, 10, 188, 4) # 34.711553680880769
#' @seealso \code{\link{gsw_SP_from_Sstar}} does the reverse; \code{\link{gsw_Sstar_from_SA}} is similar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_Sstar_from_SP.html}
gsw_Sstar_from_SP <- function(SP, p, longitude, latitude)
{
    if (missing(SP)) stop("must supply SP")
    if (missing(p)) stop("must supply p")
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    ## check for special case that SP is a matrix defined on lon and lat
    if (is.matrix(SP)) {
        dim <- dim(SP)
        if (length(longitude) == dim[1] && length(latitude) == dim[2]) {
            ll <- expand.grid(longitude=as.vector(longitude), latitude=as.vector(latitude))
            longitude <- ll$longitude
            latitude <- ll$latitude
        }
    }
    l <- argfix(list(SP=SP, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_Sstar_from_SP",
               SP=as.double(l$SP), p=as.double(l$p), longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
    rval
}

#' Freezing temperature
#'
#' Note: as of 2014-12-23, this corresponds to the Matlab function
#' called \code{gsw_t_freezing_poly}. (The confusion arises from a
#' mismatch in release version between the Matlab and C libraries.)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param saturation_fraction saturation fraction of dissolved air in seawater
#' @return in-situ freezing temperature (ITS-90) [ deg C ]
#' @examples 
#' gsw_t_freezing(34.7118, 10) # -1.902704434299200
#' @seealso \code{\link{gsw_CT_freezing}} is the analogue for Conservative Temperature.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_t_freezing.html}
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_t_freezing_poly.html}
gsw_t_freezing <- function(SA, p, saturation_fraction=1)
{
    l <- argfix(list(SA=SA, p=p, saturation_fraction=saturation_fraction))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_t_freezing",
               SA=as.double(l$SA), p=as.double(l$p), saturation_fraction=as.double(l$saturation_fraction),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' In situ temperature from Conservative Temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return in-situ temperature (ITS-90) [ deg C ]
#' @examples 
#' gsw_t_from_CT(34.7118, 28.8099, 10) # 28.785580227725703
#' @seealso \code{\link{gsw_CT_from_t}} does the reverse.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_t_from_CT.html}
gsw_t_from_CT <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_t_from_CT",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Thermobaric coefficient (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return thermobaric coefficient wrt Conservative Temperature [ 1/(K Pa) ]
#' @examples 
#' gsw_thermobaric(34.7118, 28.8099, 10) # 1.40572143831373e-12
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_thermobaric.html}
gsw_thermobaric <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_thermobaric",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Turner angle and density ratio
#'
#' This uses the 48-term density equation. The values of Turner Angle
#' Tu and density ratio Rrho are calculated at mid-point pressures, p_mid.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return list containing Tu [ degrees ], Rsubrho [ unitless ], and p_mid [ dbar ]
#' @examples
#' SA = c(34.7118, 34.8915)
#' CT = c(28.8099, 28.4392)
#' p =  c(     10,      50)
#' r <- gsw_Turner_Rsubrho(SA, CT, p) # -2.064830032393999, -0.9304018848608, 30
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_Turner_Rsubrho.html}
gsw_Turner_Rsubrho <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    r <- .C("wrap_gsw_Turner_Rsubrho",
            SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
            n=n, Tu=double(n-1), Rsubrho=double(n-1), p_mid=double(n-1))
    Tu <- r$Tu
    Rsubrho <- r$Rsubrho
    p_mid <- r$p_mid
    if (is.matrix(SA)) {
        stop("gsw_Turner_Rsubrho() cannot handle matrix SA")
        ## dim(Tu) <- dim(SA)
        ## dim(Rsubrho) <- dim(SA)
        ## dim(p_mid) <- dim(SA)
    }
    list(Tu=Tu, Rsubrho=Rsubrho, p_mid=p_mid)
}

#' Height from pressure (48-term equation)
#' 
#' @param p sea pressure [ dbar ]
#' @param latitude latitude in decimal degrees north [ -90 ... +90 ]
#' @return height [ m ]
#' @examples
#' gsw_z_from_p(10, 4) # -9.9445831334188
#' @seealso
#' This is (almost) the reverse of \code{\link{gsw_p_from_z}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_z_from_p.html}
gsw_z_from_p <- function(p, latitude)
{
    l <- argfix(list(p=p, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_z_from_p",
               p=as.double(l$p), lat=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(p))
        dim(rval) <- dim(p)
    rval
}

