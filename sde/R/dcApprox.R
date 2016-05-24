dcEuler <- function(x, t, x0, t0, theta, d, s, log=FALSE){
 Dt <- t-t0
 lik <- dnorm(x, mean = x0 + d(t0, x0, theta)*Dt, sd= sqrt(Dt)*s(t0,x0,theta),log=log)
 lik[is.infinite(lik)] <- NA
 lik
}

dcElerian <- function(x, t, x0, t0, theta, d, s, sx, log=FALSE){
 Dt <- t-t0
 A <- s(t0, x0, theta)*sx(t0, x0, theta)*Dt/2
 B <- -s(t0, x0,theta)/(2*sx(t0, x0,theta)) + x0 + d(t0, x0, theta)*Dt - A
 z <- (x-B)/A
 z[z<0] <- NA
 C <- 1/((s(t0, x0,theta)^2)*Dt)
 tmp <- sqrt(C*z)
 tmp2 <- numeric(length(tmp))
 idx <- which(abs(tmp)>500)
 tmp2[idx] <- tmp[idx] - log(2)
 tmp2[-idx] <- log(cosh(tmp[-idx]))
 lik <- -(C+z)/2 + tmp2 -0.5*log(z) - log(abs(A)) -log(2*pi)
 if(!log)
  lik <- exp(lik)
 lik[is.infinite(lik)] <- NA
 idx <- which(A==0)
 lik[idx] <- dcEuler(x, t, x0, t0, theta, d, s, log)[idx]
 lik
}

dcKessler <- function(x, t, x0, t0, theta, d, dx, dxx, s, sx, sxx, log=FALSE){
 Dt <- t-t0 
 mu <- d(t0,x0,theta)
 mu1 <- dx(t0,x0,theta)
 mu2 <- dxx(t0,x0,theta)
 sg <- s(t0,x0,theta)
 sg1 <- sx(t0,x0,theta)
 sg2 <- sxx(t0,x0,theta)
 Ex <- (x0 + mu*Dt + (mu*mu1+0.5*(sg^2 * mu2))*(Dt^2)/2)
 Vx <- ( x0^2 +(2*mu*x0 + sg^2)*Dt + (2*mu*(mu1*x0+mu+sg*sg1)+
	   sg^2*(mu2*x0 + 2*mu1 + sg1^2 + sg*sg2))*(Dt^2)/2 - Ex^2) 
 Vx[Vx<0] <- NA
 dnorm(x, mean = Ex, sd=sqrt(Vx),log=log)
}

dcOzaki <- function(x, t, x0, t0, theta, d, dx, s, log=FALSE){
  Dt <- t-t0
  Lx <- dx(t0,x0,theta)
  Kx <- log(1+d(t0,x0,theta)*(exp(Lx*Dt)-1)/(x0*Lx))/Dt
  Ex <- x0 + d(t0,x0,theta)/Lx*(exp(Lx*Dt)-1)
  Vx <- s(t0,x0,theta)^2 * (exp(2*Kx*Dt) -1)/(2*Kx)
  dnorm(x, mean=Ex, sd=sqrt(Vx),log=log) 
}

dcShoji <- function(x, t, x0, t0, theta, d, dx, dxx, dt, s, log=FALSE){
    Dt <- t-t0
    Lx <- dx(t0,x0,theta)
    Mx <- s(t0,x0,theta)^2 * dxx(t0,x0,theta)/2 + dt(t0,x0,theta)
    Ex <- (x0 + d(t0,x0,theta)*(exp(Lx*Dt)-1)/Lx + 
	       Mx*(exp(Lx*Dt) -1 -Lx*Dt)/Lx^2) 
    Vx <- s(t0,x0,theta)^2*(exp(2*Lx*Dt)-1)/(2*Lx)
    dnorm(x, mean=Ex, sd=sqrt(Vx),log=log) 
}
