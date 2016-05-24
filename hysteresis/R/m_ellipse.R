m_ellipse <- function(g,co=1.345,e=0.001,iter=1) {
  absres <- abs(g$residuals)
  robust_sigma <- 1.4826*median(absres)
 outlier <- (absres > co*robust_sigma)
  if (sum(outlier) > 0){
  w <- co*robust_sigma/absres*outlier+!outlier
  iteration <- direct(g$x,g$y,w)
  if (sum(((iteration$fit[[1]]-g$fit[[1]])^2)) < e)
  result <- list(iteration,"weights"=w,"iter"=iter)
    else result <- m_ellipse(list("residuals"=iteration$residuals,"fit"=iteration$fit,
                                  "x"=g$x,"y"=g$y),co,e,iter+1)
  }
  else result <- m_ellipse(list("residuals"=g$residuals,"fit"=g$residuals,
                                "x"=g$x,"y"=g$y),co,e,iter)
result
}
