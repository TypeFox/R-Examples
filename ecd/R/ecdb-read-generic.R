#' Read API for the ecdb
#' 
#' Read ecdb into data.frame. This can be accomplished by either specifying
#' the range of \code{alpha,gamma} or the cartesian product of \code{alpha,gamma}
#' point by point, or both. If both are specified, it follows a similar logic as plot 
#' how \code{x,y} is scoped by \code{xlim,ylim}. 
#'
#' @method read ecdb
#' 
#' @param object an object of ecdb class
#' @param alpha,gamma numeric vectors of points for cartesian product
#' @param alim,glim length-two numeric vectors of min and max range 
#' @param cusp numeric. Type of cusp. Only 0 and 1 are allowed.
#'                If cusp=1, read cusp data on the critical line.
#'                Reading cusp data must be done from the \code{alpha} side.
#'                Default: 0.
#' @param polar_ext logical, for polar coordinate extension: 
#'                  \code{R, theta, angle}.
#'                  Default: \code{FALSE}.
#' 
#' @return The data.frame from ECDATTR table.
#'
#' @keywords ecdb 
#'
#' @export read
#' @export read.ecdb
#' 
#' @import RSQLite
#'
### <======================================================================>
"read.ecdb" <- function(object, alpha=NULL, gamma=NULL, 
                        alim=NULL, glim=NULL, 
                        cusp=0, polar_ext=FALSE)
{
    if (!is.null(alpha) & is.null(alim)) {
        alim <- c(min(alpha), max(alpha))
    }
    if (!is.null(gamma) & is.null(glim)) {
        glim <- c(min(gamma), max(gamma))
    }
    if (cusp != 0 & cusp != 1) stop("Cusp must be 0 or 1")
    # --
	if (length(alim)==1) {
	    alim <- c(alim, alim) 
	}
	if (length(glim)==1) {
	    glim <- c(glim, glim) 
	}
	if (length(alim) != 2) {
		stop("Length of alim must be 1 or 2")
	}
	if (length(glim) != 2) {
	    if(cusp==0) stop("Length of glim must be 1 or 2")
	}
    # --
    if (! is.numeric(alim)) {
        stop("alim must be numeric")
    }
    if (! is.numeric(glim)) {
        if(cusp==0) stop("glim must be numeric")
    }
    
    C <- 1000000
	alim_m <- round(alim*C)
	glim_m <- round(glim*C)
	alim_str <- paste(alim_m, collapse=" AND ")
	glim_str <- paste(glim_m, collapse=" AND ")

    # use alim and glim to read data from ecdb
    select <- "SELECT 
			a.alpha_m/1000000.0 as alpha,
			a.gamma_m/1000000.0 as gamma,
            a.cusp,
			a.stdev,
			a.kurtosis,
			a.discr,
			a.jinv,
            a.ellipticity,
            a.const,
			a.time_stamp
		FROM ECDATTR a"

    sql <- paste(select,
           "WHERE a.alpha_m BETWEEN ", alim_str,
           "  AND a.gamma_m BETWEEN ", glim_str
            )
    if (cusp==1) {
        sql <- paste(select,
               "WHERE a.alpha_m BETWEEN ", alim_str, " AND a.cusp=1"
        )
        
    }

    df <- RSQLite::dbGetQuery(object@conn, sql)
    
    # subsetting
    if (!is.null(alpha)) {
        alpha_m <- round(alpha*C)
        df <- subset(df, round(alpha*C) %in% alpha_m)
    }
    if (!is.null(gamma)) {
        gamma_m <- round(gamma*C)
        df <- subset(df, round(gamma*C) %in% gamma_m)  
    }
    
    if (polar_ext) {
        a <- df[["alpha"]]
        r2 <- ecd.adj_gamma(df[["gamma"]])
        R <- sqrt(a^2+ r2^2)        
        t2 <- acos(a/R)
        df[["theta"]] <- ifelse(r2>=0, t2, 2*pi-t2)
        df[["R"]] <- R
        df[["angle"]] <- round(df[["theta"]]/pi*180)
    }
    
    return(df)
}
### <---------------------------------------------------------------------->
#' @rdname read.ecdb
setGeneric("read", function(object, alpha=NULL, gamma=NULL,
                            alim=NULL, glim=NULL, 
                            cusp=0, polar_ext=FALSE) standardGeneric("read"))
#' @rdname read.ecdb
setMethod("read", signature("ecdb"), read.ecdb)
### <---------------------------------------------------------------------->
