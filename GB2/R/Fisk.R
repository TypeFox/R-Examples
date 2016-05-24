
fisk <- function(z, w=rep(1, length(z))){

wmv <- function(z, w){
# weighted mean and variance of log(z)
  sw <- sum(w)
  logz <- log(z)
  mlz <- sum(w*logz)/sw
  vlz <- sum(w*(logz-mlz)^2)/sw
  return(list(mlz,vlz))
}

ab0 <- wmv(z, w)
# Initial values under Fisk
x0 <- c(pi/sqrt(3*ab0[[2]]),exp(ab0[[1]]),1,1)
return(x0)
}

fiskh <- function(z, w=rep(1, length(z)), hs=rep(1, length(z))){

wmv <- function(z, w, hs){
# weighted mean and variance of log(z)
  sw <- sum(w*hs)
  logz <- log(z)
  mlz <- sum(w*hs*logz)/sw
  vlz <- sum(w*hs*(logz-mlz)^2)/sw
  return(list(mlz,vlz))
}

ab0 <- wmv(z, w, hs)
# Initial values under Fisk
x0 <- c(pi/sqrt(3*ab0[[2]]),exp(ab0[[1]]),1,1)
return(x0)
}