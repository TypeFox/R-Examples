#' ecd: A package for the elliptic distribution.
#'
#' The ecd package provides the core class and functions to calculate 
#' the elliptic distribution. They are either generic or use \code{ecd.} namespace.
#' The lambda distribution is using \code{ecld.} namespace. SGED is considered part of ecld.
#' The option pricing API is using \code{ecop.} namespace.
#'
#' @author Stephen H-T. Lihn
#'
#' @docType package
#' @name ecd-package
#' @import xts methods polynom graphics moments parallel yaml RSQLite
#' 
NULL

# Some areas of this package require multi-core capability
cores <- switch( Sys.info()[['sysname']],
    Windows = 1,
    Linux   = parallel::detectCores(),
    Darwin  = parallel::detectCores(),
    parallel::detectCores()
    )

if (is.null(getOption("mc.cores"))) {
    options("mc.cores"=cores)
}

# MPFR default settings

if (is.null(getOption("ecd.precBits"))) {
    options("ecd.precBits"=120L)
}

# MPFR default Inf conversion, number of sigma as replacement for +/- Inf
# for integrateR and imgf

.ecd.mpfr.N.sigma <- 300

# end

