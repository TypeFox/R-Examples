"lik.ci" <-
function(psi, logL, conf=c(0.975,0.025))
{
  fit <- smooth.spline(psi, logL)
  p <- predict(fit, psi, deriv=1)
  psihat <- predict(smooth.spline(p$y, p$x), 0)$y
  logL.max <- predict(fit,psihat)$y
  logL.2 <- -predict(fit,psihat,deriv=2)$y
 logL <- sign(psihat-fit$x)*sqrt( 2*(logL.max-fit$y) ) 
  fit <- smooth.spline(logL,fit$x)
  p <- predict(fit, qnorm(conf))
  print(c("MLE and SE:      ",psihat,1/sqrt(logL.2)))
  print(" ")
  print(c("Confidence level:",conf))
  print(c("LR limits:       ",p$y))
  print(c("Normal limits:   ",psihat-qnorm(conf)/sqrt(logL.2)))
}

