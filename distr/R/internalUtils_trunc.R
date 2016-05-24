.trunc.up <- function(object, upper){
   ep <- .Machine$double.eps^2
   plN <- p(object)(upper, lower.tail = TRUE, log.p=TRUE)
   rnew <- function(n){
           q(object)(plN-rexp(n), lower.tail = TRUE, log.p=TRUE)
   }
   pnew <- function(q, lower.tail = TRUE, log.p = FALSE){
           indNA <- is.na(q)
           q[indNA] <- mean(q[!indNA])
           ind <- (q > upper-ep)
           p0 <- q*0
           p0[ind] <- if(lower.tail) 1 else 0
           if(log.p) p0[ind] <- log(p0[ind])
           q1 <- q[!ind]
           p1 <- p(object)(q1, lower.tail = TRUE,
                           log.p = TRUE) - plN
           p0[!ind] <- if(!log.p || !lower.tail) exp(p1) else p1
           if(!lower.tail) p0[!ind] <- 1-p0[!ind]
           if(log.p && !lower.tail) p0[!ind] <- log(p0[!ind])
           p0[indNA] <- NA
           return(p0)
   }
   dnew <- function(x, log = FALSE){
           indNA <- is.na(x)
           x[indNA] <- mean(x[!indNA])
           ind <- (x > upper-ep)
           d0 <- x*0
           d0[ind] <- 0
           if(log) d0[ind] <- log(d0[ind])
           x1 <- x[!ind]
           d1 <- d(object)(x1,  log = TRUE)-plN
           d0[!ind] <- if(log) d1 else exp(d1)
           d0[indNA] <- NA
           return(d0)
   }
   qnew <- function(p, lower.tail = TRUE, log.p = FALSE){
           indNA <- is.na(p)
           p[indNA] <- 0.5
           p0 <- if(log.p) exp(p) else p
           ind1 <- (p0 < -ep) | p0>1+ep
           indis0 <- .isEqual(p0,0,ep)
           indis1 <- .isEqual(p0,1,ep)
           in01 <- !(ind1|indis0|indis1)
           q0 <- 0*p
           q0[ind1] <- NA
           q0[indis1] <- if(lower.tail)
                                  upper else q(object)(0)
           q0[indis0] <- if(lower.tail)
                                  q(object)(0) else upper
           p1 <- p[in01]
           if(log.p && lower.tail) p1l <- plN + p1
           else{ if(log.p) p1 <- exp(p1)
                 p1l <- plN + if(lower.tail) log(p1) else log(1-p1)
                }
           q0[in01] <- q(object)(p1l, log.p = TRUE)
           q0[indNA] <- NA
           return(q0)

   }
   return(list(r=rnew,p=pnew,d=dnew,q=qnew))
}

.trunc.low <- function(object, lower){
   ep <- .Machine$double.eps^2
   if(is(object,"DiscreteDistribution")){
      Pl <- p.l(object)
      Qr <- q.r(object)
   }else{
      Pl <- p(object)
      Qr <- q(object)
   }
   plN <- Pl(lower,  lower.tail = FALSE, log.p = TRUE)
   rnew <- function(n){
           Qr(plN-rexp(n), lower.tail = FALSE, log.p = TRUE)
   }
   pnew <- function(q, lower.tail = TRUE, log.p = FALSE){
           indNA <- is.na(q)
           q[indNA] <- mean(q[!indNA])
           ind <- (q < lower+ep)
           p0 <- q*0
           p0[ind] <- if(lower.tail) 0 else 1
           if(log.p) p0[ind] <- log(p0[ind])
           q1 <- q[!ind]
           p1 <- p(object)(q1, lower.tail=FALSE,
                           log.p = TRUE)-plN
           p0[!ind] <- if(!log.p || lower.tail) exp(p1) else p1
           if(lower.tail) p0[!ind] <- 1-p0[!ind]
           if(log.p && lower.tail) p0[!ind] <- log(p0[!ind])
           p0[indNA] <- NA
           return(p0)
   }
   dnew <- function(x, log = FALSE){
           indNA <- is.na(x)
           x[indNA] <- mean(x[!indNA])
           ind <- (x < lower+ep)
           d0 <- x*0
           d0[ind] <- 0
           if(log) d0[ind] <- log(d0[ind])
           x1 <- x[!ind]
           d1 <- d(object)(x1,  log = TRUE)-plN
           d0[!ind] <- if(log) d1 else exp(d1)
           d0[indNA] <- NA
           return(d0)
   }
   qnew <- function(p, lower.tail = TRUE, log.p = FALSE){
           indNA <- is.na(p)
           p[indNA] <- 0.5
           p0 <- if(log.p) exp(p) else p
           ind1 <- (p0 < -ep) | p0>1+ep
           indis0 <- .isEqual(p0,0,ep)
           indis1 <- .isEqual(p0,1,ep)
           in01 <- !(ind1|indis0|indis1)
           q0 <- 0*p
           q0[ind1] <- NA
           q0[indis1] <- if(lower.tail)
                                  q(object)(1) else lower
           q0[indis0] <- if(lower.tail)
                                  lower else q(object)(1)
           p1 <- p[in01]
           if(log.p && !lower.tail) p1l <- plN + p1
           else{ if(log.p) p1 <- exp(p1)
                 p1l <- plN + if(lower.tail) log(1-p1) else log(p1) }
           q0[in01] <- q(object)(p1l, lower.tail = FALSE, log.p = TRUE)
           q0[indNA] <- NA
           return(q0)
   }
   return(list(r=rnew,p=pnew,d=dnew,q=qnew))
}
