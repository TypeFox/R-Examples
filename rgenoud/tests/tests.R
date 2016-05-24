suppressMessages(library(rgenoud))

#maximize the sin function
sin1 <- genoud(sin, nvars=1, max=TRUE, unif.seed=912821, int.seed=93058, print.level=0)
sin1$value <- signif(sin1$value,6)
sin1$par <- signif(sin1$par,6)
sin1$gradients <- signif(sin1$gradients,6)
print(sin1)

#minimize the sin function
sin2 <- genoud(sin, nvars=1, max=FALSE, unif.seed=912821, int.seed=93058, print.level=0)
sin2$value <- signif(sin2$value,6)
sin2$par <- signif(sin2$par,6)
sin2$gradients <- signif(sin2$gradients,6)
print(sin2)

#maximize a univariate normal mixture which looks like a claw
claw <- function(xx) {
  x <- xx[1]
  y <- (0.46*(dnorm(x,-1.0,2.0/3.0) + dnorm(x,1.0,2.0/3.0)) +
        (1.0/300.0)*(dnorm(x,-0.5,.01) + dnorm(x,-1.0,.01) + dnorm(x,-1.5,.01)) +
        (7.0/300.0)*(dnorm(x,0.5,.07) + dnorm(x,1.0,.07) + dnorm(x,1.5,.07))) 
  return(y)
}
claw1   <- genoud(claw, nvars=1,pop.size=3000,max=TRUE, unif.seed=912821, int.seed=93058, print.level=0)
claw1$value <- signif(claw1$value,6)
claw1$par <- signif(claw1$par,6)
claw1$gradients <- signif(claw1$gradients,6)
print(claw1)

# Maximize a bivariate normal mixture which looks like a claw.
biclaw <- function(xx) {
  mNd2 <- function(x1, x2, mu1, mu2, sigma1, sigma2, rho)
    {
      z1 <- (x1-mu1)/sigma1
      z2 <- (x2-mu2)/sigma2
      w <- (1.0/(2.0*pi*sigma1*sigma2*sqrt(1-rho*rho))) 
      w <- w*exp(-0.5*(z1*z1 - 2*rho*z1*z2 + z2*z2)/(1-rho*rho)) 
      return(w)
    }
  x1 <- xx[1]+1
  x2 <- xx[2]+1
  
  y <- (0.5*mNd2(x1,x2,0.0,0.0,1.0,1.0,0.0) +
        0.1*(mNd2(x1,x2,-1.0,-1.0,0.1,0.1,0.0) +
             mNd2(x1,x2,-0.5,-0.5,0.1,0.1,0.0) +
             mNd2(x1,x2,0.0,0.0,0.1,0.1,0.0) +
             mNd2(x1,x2,0.5,0.5,0.1,0.1,0.0) +
             mNd2(x1,x2,1.0,1.0,0.1,0.1,0.0)))
  
  return(y)
}
biclaw1 <- genoud(biclaw, default.domains=20, nvars=2,pop.size=5000,max=TRUE, unif.seed=912821, int.seed=93058, print.level=0)
biclaw1$value <- signif(biclaw1$value,6)
biclaw1$par <- signif(biclaw1$par,6)
biclaw1$gradients <- signif(biclaw1$gradients,6)
print(biclaw1)
