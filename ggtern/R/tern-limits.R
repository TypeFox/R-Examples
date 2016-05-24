#' Restrict Ternary Limits
#' 
#' \code{tern_limits} (or its aliasses) appends new \code{T}, \code{L} and \code{R} ternary continuous scales, 
#' where the maximum scale value is specified, and, where the minimums for each are solved.
#' 
#' The contra value (ie minimum value) for the \code{T}, \code{L} and \code{R} species is solved using
#' linear equations, therefore, if the solution is degenerate, or, the solution results in a zero range in either
#' of the proposed scales, then a warning message will be reported and an empty list returned. Note that 
#' \code{limits_tern(\dots), limit_tern(\dots)} and \code{tern_limit(\dots)} are all aliasses for 
#' the main function, \code{tern_limits(\dots)} and can be used interchangeably.
#' 
#' @return Either an empty list (when no solution can be found), or a list containing one of each 
#' of \code{scale_X_continuous} (\code{X = T, L, R})
#' 
#' @param T,L,R numeric value (scalar) of the maximum \code{T,L,R} species limit for each scale respectively
#' @param ... other arguments to pass to ALL of \code{scale_X_continuous} (\code{X = T, L, R})
#' 
#' @seealso \code{\link{scale_T_continuous}}, \code{\link{scale_L_continuous}} and \code{\link{scale_R_continuous}}
#' @author Nicholas Hamilton
#' @examples 
#' df = data.frame(x=runif(10),y=runif(10),z=runif(10))
#' ggtern(df,aes(x,y,z)) + geom_point() + tern_limits(0.7,0.3,0.4)
#' @name tern_limits
#' @rdname tern_limits
NULL

#' @rdname tern_limits
#' @export
tern_limits <- function(T=1,L=1,R=1,...){
  ret <- list()
  if(!all(sapply(list(T,L,R),function(x){length(x) == 1 && is.numeric(x)})))
    stop("Arguments T, L and R must be numeric and scalar",call.=FALSE)
  tryCatch({
    s    <- round(solve(diag(-1,3,3) + 1, c(1-T,1-L,1-R)),3)
    o    <- function(x){x[order(x)]}
    T    <- o(c(s[1],T)); L = o(c(s[2],L)); R = o(c(s[3],R))
    lims <- list(T,L,R)
    if(any(sapply(lims,function(x){diff(x) == 0})))
      stop("Invalid limits, solution produces zero ranges on some scales",call.=FALSE)
    if(any(sapply(lims,function(x){min(x) < 0 | max(x) > 1})))
      warning("Solution to limits produces range outside of [0,1] for some scales",call.=FALSE)
    ret <- list( scale_T_continuous(limits=T,...),
                 scale_L_continuous(limits=L,...),
                 scale_R_continuous(limits=R,...))
  },error=function(e){ warning(e)  })
  invisible(ret)
}

#'@rdname tern_limits
#'@export 
limits_tern <- function(...){tern_limits(...)}

#'@rdname tern_limits
#'@export 
limit_tern <- function(...){tern_limits(...)}

#'@rdname tern_limits
#'@export 
tern_limit <- function(...){tern_limits(...)}


