"jointCOP" <-
function(t, cop=NULL, para=NULL, type=c("and", "or"), ...) {

  type <- match.arg(type)

  if(length(t) > 1) {
     warning("only the first value in 't' will be used")
     t <- t[1]
  }

  "orfunc" <- function(v, LHS=NULL, cop=NULL, para=NULL, ...) {
      u <- v # duplication is made so that the heredity of the next line is seen
      RHS <- duCOP(u, v, cop = cop, para = para, ...)
      #message("u=",u,",  v=",v, ",  LHS=",LHS, ",  RHS=",RHS, ", DELTA=", LHS-RHS)
      return(LHS - RHS)
  }
  my.rt <- NULL
  if(type == "and") {
     u <- v <- diagCOPatf(t, cop=cop, para=para)
     zz <- c(u, v, t)
     names(zz) <- c("U", "V", "jointANDprob")
     return(zz)
  } else if(type == "or") {
     lo <- duCOP(0,0, cop=cop, para=para, ...)
     if(is.na(lo) | ! is.finite(lo)) lo <- .Machine$double.eps
     try(my.rt <- uniroot(orfunc, interval=c(lo,t), LHS=t,
              cop=cop, para=para, tol=.Machine$double.eps/10, ...), silent=FALSE)
     if(is.null(my.rt)) {
        u <- v <- NA
     } else if(length(my.rt$root) != 0) {
        u <- v <- my.rt$root
     } else {
        u <- v <- NA
     }
     zz <- c(u, v, t)
     names(zz) <- c("U", "V", "jointORprob")
     return(zz)
  } else {
     stop("should not be here in logic")
  }
}
