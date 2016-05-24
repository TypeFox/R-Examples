#' Enrich a basic \code{ecdattr} object
#' 
#' It takes a basic \code{ecdattr} object, enrich it with ecd attributes.
#' This function is computationally heavy. So the objects are often wrapped in a list 
#' and computed via \code{parallel::mclapply}.
#'
#' @param p a basic \code{ecdattr} object
#' 
#' @return an enriched \code{ecdattr} object
#'
#' @keywords ecdattr
#'
#' @export
#'
#' @importFrom utils str
#'
### <======================================================================>
ecdattr.enrich <- function(p) {
    if (class(p) != "ecdattr") {
        stop(paste("pair is not a ecdattr:", str(p)))
    }
    
    d <- NULL
    sigma <- 1
    if(p@use.mpfr) sigma <- ecd.mpfr(1)
    
    if (p@cusp > 0) {
        d <- ecd.cusp(alpha=p@alpha, sigma=sigma)
    } else {
        d <- ecd(alpha=p@alpha, gamma=p@gamma, sigma=sigma)
    }
    
    if (length(d@stats)==0) {
        stop("stats is not computed in ecd object (d)")
    }
    st <- d@stats
    
    j <- jinv(d)
    j <- ifelse(is.nan(j) | is.infinite(j), "NULL", ecd.mp2f(j))
    
    el <- ellipticity(d)
    xe <- "NULL"
    if (is.list(el) & "avg" %in% names(el)) {
        xe <- ellipticity(d)$avg
        xe <- ifelse(is.nan(xe) | is.infinite(xe), "NULL", ecd.mp2f(xe))
    }

    a <- list()
    a$stdev    <- ecd.mp2f(st$stdev)
    a$kurtosis <- ecd.mp2f(st$kurtosis)
    a$discr    <- ecd.mp2f(discr(d))
    a$jinv     <- j
    a$ellipticity <- xe
    a$const    <- ecd.mp2f(d@const)
    a$time_stamp  <- as.integer(Sys.time())
    
    p@ecd <- d
    p@attr <- a
    p@enriched <- TRUE
    
    return(p)
}


