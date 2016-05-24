.onAttach <- function(libname, pkgname){
    options(locatorBell = FALSE)# no effect for R versions prior to 1.8.0
    ##FIXME(?): should we warn when changing a global option
}

## running in R or S+ ? -- the following even works in S(not plus):
if(!exists("version") || is.null(vl <- version$language) || vl != "R")
{ ## we are not in R :
    is.R <- function () FALSE
    if(!exists("version") || version$major < 6)
        "%in%" <- function(x, table) !is.na(match(x, table, nomatch = NA))

    dev.interactive <- function () {
	interactive() && exists(".Device") &&
	.Device %in%
        c("motif", "graphsheet", "java.graph", "X11", "openlook")
    }

    which.min <- function(x) sort.list(x)[1]
}
if(exists("version")) rm(vl)

##----------------------------------------------------------------------------
if(!exists("cov2cor", mode = "function"))## for S+ ; is in R 1.8.0 and later
cov2cor <- function(V)
{
    ## Purpose: Covariance matrix |--> Correlation matrix -- efficiently
    ## ----------------------------------------------------------------------
    ## Arguments: V: a covariance matrix (i.e. symmetric and positive definite)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 12 Jun 2003, 11:50
    p <- (d <- dim(V))[1]
    if(!is.numeric(V) || length(d) != 2 || p != d[2])
	stop("`V' is not a square numeric matrix")
    Is <- sqrt(1/diag(V)) # diag( 1/sigma_i )
    if(any(!is.finite(Is)))
	warning("diagonal has non-finite entries")
    r <- V # keep dimnames
    r[] <- Is * V * rep(Is, each = p)
    ##	== D %*% V %*% D  where D = diag(Is)
    r[cbind(1:p,1:p)] <- 1 # exact in diagonal
    r
}
