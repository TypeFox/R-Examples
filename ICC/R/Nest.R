Nest<-function(est.type = c("hypothetical", "pilot"), w, ICC = NULL, k = NULL, x = NULL, y = NULL, data = NULL, alpha = 0.05){
  icall <- list(y = substitute(y), x = substitute(x))
  z <- qnorm(1-alpha/2, mean = 0, sd = 1)
  type <- match.arg(est.type)
  if(type == "hypothetical") {
    n.est <- matrix(nrow = length(ICC), ncol = length(k))   
    for (i in 1:length(ICC)){
      for(j in 1:length(k)){
        n.est1 <- 8*(z^2)*(((1-ICC[i])^2)*((1+(k[j]-1)*ICC[i])^2))/(k[j]*(k[j]-1)*(w^2))+1 
        n.est2 <- ceiling(n.est1)
        n.est[i,j] <- n.est2
      }
    }
    n.est.table <- data.frame(n.est, row.names = ICC)
    names(n.est.table) <- k
    return(n.est.table)
    } else {
         square<-function(z){z^2}
         if(is.character(icall$y)){
           warning("passing a character string to 'y' is deprecated since ICC vesion 2.3.0 and will not be supported in future versions. The argument to 'y' should either be an unquoted column name of 'data' or an object")
           if(missing(data)) stop("Supply either the unquoted name of the object containing 'y' or supply both 'data' and then 'y' as an unquoted column name to 'data'")
           icall$y <- eval(as.name(y), data, parent.frame())
         } 
         if(is.name(icall$y)) icall$y <- eval(icall$y, data, parent.frame())
         if(is.call(icall$y)) icall$y <- eval(icall$y, data, parent.frame())
         if(is.character(icall$y)) icall$y <- eval(as.name(icall$y), data, parent.frame())


         if(is.character(icall$x)){
           warning("passing a character string to 'x' is deprecated since ICC vesion 2.3.0 and will not be supported in future versions. The argument to 'x' should either be an unquoted column name of 'data' or an object")
           if(missing(data)) stop("Supply either the unquoted name of the object containing 'x' or supply both 'data' and then 'x' as an unquoted column name to 'data'")
           icall$x <- eval(as.name(x), data, parent.frame())
         } 
         if(is.name(icall$x)) icall$x <- eval(icall$x, data, parent.frame())
         if(is.call(icall$x)) icall$x <- eval(icall$x, data, parent.frame())
         if(is.character(icall$x) && length(icall$x) == 1) icall$x <- eval(as.name(icall$x), data, parent.frame())


         tdata <- data.frame(icall)
         ICC.results <- ICCest(x = x, y = y, data = tdata)
         n.est1b <- 8*(z^2)*(((1-ICC.results$ICC)^2)*((1+(ICC.results$k-1)*ICC.results$ICC)^2))/(ICC.results$k*(ICC.results$k-1)*(w^2))+1
         n.est2b <- ceiling(n.est1b)
         return(n.est2b)
      }
}
