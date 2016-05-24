"joint.curvesCOP" <-
function(cop=NULL, para=NULL, type=c("and", "or"),
         probs=c(0.5, 0.8, 0.90, 0.96, 0.98, 0.99, 0.995, 0.998),
         zero2small=TRUE, small=1E-6, divisor=100, delu=0.001, ...) {

  type <- match.arg(type)

  "orfunc" <- function(v, u=NULL, LHS=NULL, cop=NULL, para=NULL, ...) {
      RHS <- duCOP(u, v, cop = cop, para = para, ...)
      #message("u=",u,",  v=",v, ",  LHS=",LHS, ",  RHS=",RHS, ", DELTA=", LHS-RHS)
      return(LHS - RHS)
  }

  zz <- new.env()
  for(t in probs) {
    if(type == "and") {
       # Concerning the tmp handling, there seems to be cases in which empirical copulas
       # are being used with small sample sizes that can trigger a 'wrong sign in 'by' argument
       # error on the following sequence, if so, let us try reversing the sequence and if that
       # fails then bail out enterly.
       delta <- diff(range(c(t+delu, 1-delu)))/divisor
       tmp <- NULL
       try(tmp <- seq(t+delu, 1-delu, by=delta), silent=FALSE)
       if(is.null(tmp)) {
          warning("trying to compensate for 'by' error by reversing")
          try(tmp <- seq(t+delu, 1-delu, by=-delta), silent=FALSE)
       }
       u <- c(t, t+delu/100, t+delu/10, t+delu/5, t+delu/2,
               tmp,
              1-delu/2, 1-delu/5, 1-delu/10, 1-delu/100, 1)
       u  <- sort(unique(u))
       v  <- sapply(1:length(u), function(i) { COPinv(cop=cop, u[i], t, para=para) })
       if(zero2small) v[v == 0] <-     small
       if(zero2small) v[v == 1] <- 1 - small
       uv <- data.frame(U=u, V=v)
       assign(as.character(t), uv, envir=zz)
    } else if(type == "or") {
       delta <- diff(range(c(0, t)))/divisor
       u <- c(0, delu/100, delu/10, delu/5, delu/2,
             seq(delu,t, by=delta),
             t-delu/2, t-delu/5, t-delu/10, t-delu/100, t)
       u <- sort(unique(u))
       v <- vector(mode="numeric", length(u))
       for(i in 1:length(u)) {
          lo <- duCOP(0,0, cop=cop, para=para, ...)
          if(is.na(lo) | ! is.finite(lo)) lo <- .Machine$double.eps
          my.rt <- NULL
          try(my.rt <- uniroot(orfunc, interval=c(lo,t), u=u[i], LHS=t,
                               cop=cop, para=para, ...), silent=TRUE)
          if(is.null(my.rt)) {
             v[i] <- NA
          } else if(length(my.rt$root) != 0) {
             v[i] <- my.rt$root
          } else {
             v[i] <- NA
          }
       }
       if(zero2small) v[v == 0] <-     small
       if(zero2small) v[v == 1] <- 1 - small
       uv <- data.frame(U=u, V=v)
       assign(as.character(t), uv, envir=zz)
    } else {
       stop("should not be here in logic")
    }
  }
  zzz <- as.list(zz); rm(zz)
  return(zzz)
}
