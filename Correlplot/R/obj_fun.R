obj_fun <- function(x,aux) {
   FUN <- match.fun(aux$FUN)
   theta <- c(0,x)
   CC <- matrix(rep(theta,aux$n),nrow=aux$n)
   CC <- FUN(CC-t(CC))
   y <- sum((aux$R-CC)^2)
   return(y)
}
