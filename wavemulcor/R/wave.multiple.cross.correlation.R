wave.multiple.cross.correlation <-
function(xx, lag.max=NULL, ymaxr=NULL) { 
  sum.of.squares <- function(x) { sum(x^2, na.rm=TRUE) / sum(!is.na(x)) } 
  sum.of.not.squares <- function(x) { sum(x, na.rm=TRUE) / sum(!is.na(x)) } 
  d <- length(xx) 
  dd <- d*(d-1)/2 
  l <- length(xx[[1]])
  if(is.null(lag.max)) {lag.max=trunc(sqrt(length(xx[[1]][[l]]))/2)}
  lm <- min(length(xx[[1]][[l]])-1, lag.max, na.rm=TRUE) 
  x.var <- vector("list", d) 
  for(j in 1:d) { 
  x.var[[j]] <- unlist(lapply(xx[[j]], sum.of.squares)) } 
  xy.cor <- vector("list", dd) 
  xy <- vector("list", l) 
  jk <- 0 
  for(k in 1:(d-1)) { 
     for(j in (k+1):d) { 
     jk <- jk+1 
     for(i in 1:l) { 
        xy[[i]] <- as.vector(xx[[j]][[i]] * xx[[k]][[i]]) 
     } 
     xy.cov <- unlist(lapply(xy, sum.of.not.squares)) 
     xy.cor[[jk]] <- xy.cov / sqrt(x.var[[j]] * x.var[[k]]) 
  }} 
  xy.cor.vec <- matrix(unlist(xy.cor),l,dd) 
  xy.mulcor <- matrix(0, l, 2*lm+1)
  YmaxR <- vector("numeric",l) 
  for(i in 1:l) {
     r <- xy.cor.vec[i,] 
     P <- diag(d)/2 
     P[lower.tri(P)] <- r 
     P <- P+t(P) 
     Pidiag <- diag(solve(P))
     if(is.null(ymaxr)) { 
       YmaxR[i] <- Pimax <- which.max(Pidiag) ## detect i | x[i] on rest x gives max R2
     } else {YmaxR[i] <- Pimax <- ymaxr}
     xy.mulcor[i,lm+1] <- sqrt(1-1/Pidiag[Pimax]) ## lag=0: this must be same as in wave.multiple.correlation 
     if(lm>0) {
       y <- xx[[Pimax]][[i]]
       z <- xx[[Pimax]][[i]]
       vlength <- length(y)
       x <- xx[-Pimax][[1]][[i]]
       if(d>2) { for(k in 2:(d-1)) { x <- cbind(x,xx[-Pimax][[k]][[i]]) } }
       for(j in 1:lm) {       ## now we obtain R2 of var[Pimax] with lagged values 
         y <- c(y[2:vlength], NA) 
         z <- c(NA, z[1:(vlength-1)]) 
         xy.mulcor[i,lm+1+j] <- sqrt( summary(lm(formula = y ~ x))$r.squared )
         xy.mulcor[i,lm+1-j] <- sqrt( summary(lm(formula = z ~ x))$r.squared )
     }}
  }
#browser()
Lst <- list(xy.mulcor=as.matrix(xy.mulcor),YmaxR=YmaxR)
return(Lst)
}

