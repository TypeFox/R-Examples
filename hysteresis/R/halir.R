halir <- function(g,sigma.sq=0,co=1.345,e=0.001,iter=1) {
  absres <- abs(g$residuals)
  robust_sigma <- 1.4826*median(absres)
 outlier <- (absres > co*robust_sigma)
  w <- co*robust_sigma/absres*outlier+!outlier
  iteration <- direct_renorm(g$x,g$y,w,sigma.sq)
  sigma.sq <- c(sigma.sq+abs(iteration$fit[[4]]/iteration$fit[[3]]))
  if (sum(((iteration$fit[[1]]-g$fit[[1]])^2)) < e & iter > 2)
  result <- list(iteration,"weights"=w,"iter"=iter,"sigma.sq"=sigma.sq)
    else result <- halir(list("residuals"=iteration$residuals,"fit"=iteration$fit,
                                  "x"=g$x,"y"=g$y),sigma.sq,co,e,iter+1)
  

result
}
