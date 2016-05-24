"joint.curvesCOP2" <-
function(cop=NULL, para=NULL, type=c("and", "or"),
         probs=c(0.5, 0.8, 0.90, 0.96, 0.98, 0.99, 0.995, 0.998),
         zero2small=TRUE, small=1E-6, divisor=100, delv=0.001, ...) {

  type <- match.arg(type)

  "orfunc" <- function(u, v=NULL, LHS=NULL, cop=NULL, para=NULL, ...) {
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
       delta <- diff(range(c(t+delv, 1-delv)))/divisor
       tmp <- NULL
       try(tmp <- seq(t+delv, 1-delv, by=delta), silent=FALSE)
       if(is.null(tmp)) {
          warning("trying to compensate for 'by' error by reversing")
          try(tmp <- seq(t+delv, 1-delv, by=-delta), silent=FALSE)
       }
       v <- c(t, t+delv/100, t+delv/10, t+delv/5, t+delv/2,
           tmp,
           1-delv/2, 1-delv/5, 1-delv/10, 1-delv/100, 1)
       v  <- sort(unique(v))
       u  <- sapply(1:length(v), function(i) { COPinv(cop=cop, v[i], t, para=para) })
       if(zero2small) u[u == 0] <-     small
       if(zero2small) u[u == 1] <- 1 - small
       uv <- data.frame(U=u, V=v)
       assign(as.character(t), uv, envir=zz)
    } else if(type == "or") {
       delta <- diff(range(c(0, t)))/divisor
       v <- c(0, delv/100, delv/10, delv/5, delv/2,
             seq(delv,t, by=delta),
             t-delv/2, t-delv/5, t-delv/10, t-delv/100, t)
       v <- sort(unique(v))
       u <- vector(mode="numeric", length(v))
       for(i in 1:length(v)) {
          lo <- duCOP(0,0, cop=cop, para=para, ...)
          if(is.na(lo) | ! is.finite(lo)) lo <- .Machine$double.eps
          my.rt <- NULL
          try(my.rt <- uniroot(orfunc, interval=c(lo,t), v=v[i], LHS=t,
                               cop=cop, para=para, ...), silent=TRUE)
          if(is.null(my.rt)) {
             u[i] <- NA
          } else if(length(my.rt$root) != 0) {
             u[i] <- my.rt$root
          } else {
             u[i] <- NA
          }
       }
       if(zero2small) u[u == 0] <-     small
       if(zero2small) u[u == 1] <- 1 - small
       uv <- data.frame(U=u, V=v)
       assign(as.character(t), uv, envir=zz)
    } else {
       stop("should not be here in logic")
    }
  }
  zzz <- as.list(zz); rm(zz)
  return(zzz)
}
