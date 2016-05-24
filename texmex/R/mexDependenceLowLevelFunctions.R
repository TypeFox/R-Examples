# check constraints on parameters under constrained Laplace estimation
ConstraintsAreSatisfied <- function(a,b,z,zpos,zneg,v){
   C1e <- a <= min(1, 1 - b*min(z)*v^(b-1), 1 - v^(b-1)*min(z) + min(zpos)/v) &
  	      a <= min(1, 1 - b*max(z)*v^(b-1), 1 - v^(b-1)*max(z) + max(zpos)/v)
				   
   C1o <- a <= 1 & 
          a > 1 - b * min(z) * v^(b-1) & 
          a > 1 - b * max(z) * v^(b-1) &
   		   (1 - 1/b)*(b*min(z))^(1/(1-b)) * (1-a)^(-b/(1 - b)) + min(zpos) > 0 &
			   (1 - 1/b)*(b*max(z))^(1/(1-b)) * (1-a)^(-b/(1 - b)) + max(zpos) > 0

	 C2e <- -a <= min(1, 1 + b*v^(b-1)*min(z), 1 + v^(b-1)*min(z) - min(zneg)/v) &
		      -a <= min(1, 1 + b*v^(b-1)*max(z), 1 + v^(b-1)*max(z) - max(zneg)/v)

   C2o <- -a <= 1 & 
          -a > 1 + b*v^(b-1)*min(z) & 
          -a > 1 + b*v^(b-1)*max(z) &
		     (1-1/b)*(-b*min(z))^(1/(1-b))*(1+a)^(-b/(1-b)) - min(zneg) > 0 &
		     (1-1/b)*(-b*max(z))^(1/(1-b))*(1+a)^(-b/(1-b)) - max(zneg) > 0

   if (any(is.na(c(C1e, C1o, C2e, C2o)))) {
       warning("Strayed into impossible area of parameter space")
       C1e <- C1o <- C2e <- C2o <- FALSE
   }

   (C1e | C1o) && (C2e | C2o)
}

# positive dependence Gumbel and pos or neg dependence Laplace neg likelihood function
PosGumb.Laplace.negloglik <- function(yex, ydep, a, b, m, s, constrain, v, aLow) {
  BigNumber <- 10^40
  WeeNumber <- 10^(-10)

  if(a < aLow[1] | s < WeeNumber | a > 1-WeeNumber  | b > 1-WeeNumber) {
    res <- BigNumber
  } else {
    mu <- a * yex + m * yex^b
    sig <- s * yex^b

    res <- sum(0.5 * log(2*pi) + log(sig) + 0.5 * ((ydep - mu)/sig)^2)

    if (is.infinite(res)){
        if (res < 0){ 
          res <- -BigNumber 
        } else {
          res <- BigNumber
        }
        warning("Infinite value of Q in mexDependence")
    } else if (constrain){
		   #v <- v * max(yex)
			 zpos <- range(ydep - yex) # q0 & q1
			 z <- range((ydep - yex * a) / (yex^b)) # q0 & q1 
			 zneg <- range(ydep + yex) # q0 & q1
				   
			 if (!ConstraintsAreSatisfied(a,b,z,zpos,zneg,v)){
          res <- BigNumber
       }
		} 
  }    
  res
}

PosGumb.Laplace.negProfileLogLik <- function(yex, ydep, a, b, constrain, v, aLow) { 
  Z <- (ydep - yex * a) / (yex^b)

  m <- mean(Z)
  s <- sd(Z)

  res <- PosGumb.Laplace.negloglik(yex,ydep,a,b,m=m,s=s,constrain,v,aLow=aLow)
  res <- list(profLik=res,m=m, s=s)
  res
}
