epane <- function(X){
XX <- X^2
ep <- (3/4)*(1-XX)*(XX<=1)
return(ep)
}
gaussien <- function(X){
XX <- X^2
ga <- 1/sqrt(2*pi)*exp(-0.5*XX)
return(ga)
}
quartic <- function(X){
XX <- X^2
qu <- (15/16)*((1-XX)^2)*(XX<=1)
return(qu)
}
uniform <- function(X){
XX<-X^2
un<-1*(XX<=1)
return(un)
}
