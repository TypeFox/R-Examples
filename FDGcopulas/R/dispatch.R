## Generators ##

g.sinus <- function(t,a,b){
  if(round(a,5)==round(b,5)){ 
    ( 4*sin(t*a)^2+t*a*(2*(1-t)*a+sin(2*a)-sin(2*t*a)) ) / ( 4*sin(a)^2 )
  }else{    
    term1 <- sin(a*t)*sin(b*t)/(sin(a)*sin(b))
    term2 <- t*a*b/(2*(b^2-a^2)*sin(a)*sin(b))
    term3 <- (a+b)*(sin((a-b)*t)+sin(b-a))+(b-a)*(sin(a+b)-sin((a+b)*t))
    term1 + term2*term3
  }
}
g.exponential <- function(t,a,b){
  if(b<=0){
    if(a<=0){
      1
    }else{
      exp((t^a-1)/a)
    }
  }else{
    if(a<=0){
      exp((t^b-1)/b)
    }else{
      h <- function(t,c){ (t^c-1)/c }
      integrand <- function(x){ exp(h(x,a)+h(x,b))*x^(a+b-2) }
      exp(h(t,a)+h(t,b))+t*integrate(integrand,t,1)$value
    }
  }
}
  


## Pairwise dependence coefficients
# Spearman's rho
rho.frechet <- function(a,b) a*b
rho.cuadrasauge <- function(a,b){ p <- a*b; s <- a+b; 3*p/(5-s) }
rho.sinus<-function(a,b){
    integrand<-function(x,a,b){
        x^2*sin(a*x)*sin(b*x)+.25*a*b*x^4*cos(a*x)*cos(b*x)
    }
    12*integrate(integrand,0,1,a=a,b=b)$value/(sin(a)*sin(b))-3
}
rho.exponential<-function(a,b){
  if(b<=0){
    if(a<=0){
      1
    }else{
      12*integrate(function(x)x^2*exp((x^a-1)/a),0,1)$value-3
    }
  }else{
    if(a<=0){
      12*integrate(function(x)x^2*exp((x^b-1)/b),0,1)$value-3
    }else{
      integrand<-function(x,a,b){
        (x^2+.25*x^(2+a+b))*exp((x^a-1)/a+(x^b-1)/b)
      }
      12*integrate(integrand,0,1,a=a,b=b)$value-3
    }
  }
}
rho.EV.caf <- function(a,b){
  c <- a*b
  3*c/(4-c)
}
rho.EV.sinus <- function(a,b){
  da <- a/tan(a)
  db <- b/tan(b)
  c <- (1-da)*(1-db)
  3*c/(4-c)
}
rho.EV.exponential <- function(a,b){ 0 }

# Kendall's tau
tau.frechet <- function(a,b){ p <- a*b; p*(p+2)/3 }
tau.cuadrasauge <- function(a,b){ p <- a*b; s <- a+b; p*(p+6-2*s)/(s^2-8*s+15) }
tau.sinus <- function(a,b){
  4*integrate(function(x)x*g.sinus(x,a,b)^2, 0, 1)$value-1
}
tau.exponential <- function(a,b){
  4*integrate(function(x)x*g.exponential(x,a,b)^2, 0, 1)$value-1
}

tau.EV.caf <- function(a,b){ p <- a*b; p/(2-p) }
tau.EV.sinus <- function(a,b){
  da <- a/tan(a)
  db <- b/tan(b)
  c <- (1-da)*(1-db)
  c/(2-c)
}
tau.EV.exponential <- function(a,b){ 0 }

# Upper tail dependence coefficient
utdc.sinus <- function(a,b){ (1-a/tan(a))*(1-b/tan(b)) }
utdc.frechet <- function(a,b) a*b
utdc.cuadrasauge <- function(a,b) a*b
utdc.exponential <- function(a,b) 0
utdc.EV.caf <- function(a,b) a*b
utdc.EV.sinus <- function(a,b) (1-a/tan(a))*(1-b/tan(b))
utdc.EV.exponential <- function(a,b) 0
  
# Lower tail dependence coefficient
ltdc.sinus <- function(a,b) 0
ltdc.frechet <- function(a,b) a*b
ltdc.cuadrasauge <- function(a,b) 0
ltdc.exponential <- function(a,b) exp(-1/a-1/b)
ltdc.EV.sinus <- function(a,b) 0
ltdc.EV.caf <- function(a,b) 0
ltdc.EV.exponential <- function(a,b) 0



## Partial derivative of r(a,b) where
# r(a,b) is the dependence coefficient viewed as a function of the
# parameters a and b

# Spearman's rho 
tfun.frechet.rho <- function(a,b) b
tfun.cuadrasauge.rho <- function(a,b) (5-b)*3*b/(5-a-b)^2
tfun.sinus.rho <- function(a,b){    
    integrand1<-function(x,a,b){
        x^2*sin(a*x)*sin(b*x)+.25*a*b*x^4*cos(a*x)*cos(b*x)
    }
    integrand2<-function(x,a,b){
        x^3*cos(a*x)*sin(b*x)+.25*x^4*b*cos(b*x)*(cos(a*x)-a*x*sin(a*x))
    }
    term1 <- integrate(integrand1,0,1,a=a,b=b)$value*cos(a)*sin(b)
    term2 <- integrate(integrand2,0,1,a=a,b=b)$value*sin(b)*sin(a)
    12*(term2-term1)/(sin(a)*sin(b))^2
}
tfun.exponential.rho <- function(a,b){
  ## integrand<-function(x,a,b){
  ##     ((a*x^a*log(x)-x^a+1)*(x^2+.25*x^(2+a+b))/a^2+.25*x^(2+a+b)*log(x))*exp((x^a-1)/a+(x^b-1)/b)
  ## }
  ## 12*integrate(integrand,0,1,a=a,b=b)$value
  rho.exponential.wrapper <- function(x) rho.exponential(x[1],x[2])
  numDeriv::grad(rho.exponential.wrapper,c(a,b))[1]
}
tfun.EV.caf.rho <- function(a,b) 12*b/(4-a*b)^2
tfun.EV.sinus.rho <- function(a,b){
  x <- (1-a/tan(a))*(1-b/tan(b))
  y <- (1-b/tan(b))*(a*(1+tan(a)^2)-tan(a))/tan(a)^2
  12*y/(4-x)^2
  ## num <- (1-b/tan(b))*12*(tan(a)-a*(1+tan(a)^2))
  ## denom <- (tan(a)*(4-(1-a/tan(a))*(1-b/tan(b))))^2
  ## num/denom
}
tfun.EV.exponential.rho <- function(a,b) 0



# Kendall's tau
tfun.frechet.tau <- function(a,b) 2*(a*b+1)*b/3
tfun.cuadrasauge.tau <- function(a,b){
  s <- a+b
  p <- a*b
  t1 <- (6+2*p-4*a-2*b)*b/(s^2-8*s+15)
  t2 <- (p*(6+p-2*s)*(2*s-8))/(s^2+15-8*s)^2
  t1+t2
}
tfun.sinus.tau <- function(a,b){
  tau.sinus.wrapper <- function(x) tau.sinus(x[1],x[2])
  numDeriv::grad(tau.sinus.wrapper,c(a,b))[1]
}
tfun.exponential.tau <- function(a,b){
  tau.exponential.wrapper <- function(x) tau.exponential(x[1],x[2])
  numDeriv::grad(tau.exponential.wrapper,c(a,b))[1]
}
tfun.EV.caf.tau <- function(a,b) 2/(2-a*b)^2
tfun.EV.sinus.tau <- function(a,b){
  c <- (1-a/tan(a))*(1-b/tan(b))
  (1-b/tan(b))*(2/(2-c)^2)
}
tfun.EV.exponential.tau <- function(a,b) 0
  
  
# Upper tail dependence coefficient
tfun.frechet.utdc <- function(a,b) b
tfun.cuadrasauge.utdc <- function(a,b) b
tfun.sinus.utdc <- function(a,b){
  (1-b/tan(b))*(a+a/tan(a)^2-1/tan(a))
}
tfun.exponential.utdc <- function(a,b) 0 
tfun.EV.caf.utdc <- function(a,b)  b 
tfun.EV.sinus.utdc <- function(a,b){
  (1-b/tan(b))*(a+a/tan(a)^2-1/tan(a))
}
tfun.EV.exponential.utdc <- function(a,b) 0 

# Lower tail dependence coefficient
tfun.frechet.ltdc <- function(a,b) b
tfun.cuadrasauge.ltdc <- function(a,b) b
tfun.sinus.ltdc <- function(a,b){
  (1-b/tan(b))*(a+a/tan(a)^2-1/tan(a))
}
tfun.exponential.ltdc <- function(a,b) 0 
tfun.EV.caf.ltdc <- function(a,b) 0
tfun.EV.sinus.ltdc <- function(a,b) 0
tfun.EV.exponential.ltdc <- function(a,b) 0



## Function that picks the right internal funcions according to the
# type of dependence coefficient and family
dispatch <- function(copula, depcoefType){
  if(copula@extremevalue){  
    if(depcoefType=="spearman"){
      switch(copula@family,
             "frechet" = list(tfun=tfun.EV.caf.rho,rho=rho.EV.caf),
             "cuadrasauge" = list(tfun=tfun.EV.caf.rho,rho=rho.EV.caf),
             "sinus" = list(tfun=tfun.EV.sinus.rho,rho=rho.EV.sinus),
             "exponential" = list(tfun=tfun.EV.exponential.rho,
                 rho=rho.EV.exponential))
    }else if(depcoefType=="kendall"){
      switch(copula@family,
             "frechet" = list(tfun=tfun.EV.caf.tau,tau=tau.EV.caf) ,
             "cuadrasauge" = list(tfun=tfun.EV.caf.tau,tau=tau.EV.caf),
             "sinus" = list(tfun=tfun.EV.sinus.tau,tau=tau.EV.sinus),
             "exponential" = list(tfun=tfun.EV.exponential.tau,
                 tau=tau.EV.exponential))
    }else if(depcoefType=="utdc"){
      switch(copula@family,
             "frechet" = list(tfun=tfun.EV.caf.utdc,utdc=utdc.EV.caf) ,
             "cuadrasauge" = list(tfun=tfun.EV.caf.utdc,utdc=utdc.EV.caf),
             "sinus" = list(tfun=tfun.EV.sinus.utdc,utdc=utdc.EV.sinus),
             "exponential" = list(tfun=tfun.EV.exponential.utdc,
                 utdc=utdc.EV.exponential))
    }else if(depcoefType=="ltdc"){
      switch(copula@family,
             "frechet" = list(tfun=tfun.EV.caf.ltdc,ltdc=ltdc.EV.caf) ,
             "cuadrasauge" = list(tfun=tfun.EV.caf.ltdc,ltdc=ltdc.EV.caf),
             "sinus" = list(tfun=tfun.EV.sinus.ltdc,ltdc=ltdc.EV.sinus),
             "exponential" = list(tfun=tfun.EV.exponential.ltdc,
                 ltdc=ltdc.EV.exponential))
    }else{
      stop("Wrong type of dependence coefficient")
    }
  }else{
    if(depcoefType=="spearman"){
      switch(copula@family,
             "frechet" = list(tfun=tfun.frechet.rho,rho=rho.frechet),
             "cuadrasauge" = list(tfun=tfun.cuadrasauge.rho,
                 rho=rho.cuadrasauge),
             "sinus" = list(tfun=tfun.sinus.rho,rho=rho.sinus),
             "exponential" = list(tfun=tfun.exponential.rho,
                 rho=rho.exponential))
    }else if(depcoefType=="kendall"){
      switch(copula@family,
             "frechet" = list(tfun=tfun.frechet.tau,tau=tau.frechet) ,
             "cuadrasauge" = list(tfun=tfun.cuadrasauge.tau,
                 tau=tau.cuadrasauge),
             "sinus" = list(tfun=tfun.sinus.tau,tau=tau.sinus),
             "exponential" = list(tfun=tfun.exponential.tau,
                 tau=tau.exponential))
    }else if(depcoefType=="utdc"){
      switch(copula@family,
             "frechet" = list(tfun=tfun.frechet.utdc,utdc=utdc.frechet),
             "cuadrasauge" = list(tfun=tfun.cuadrasauge.utdc,
               utdc=utdc.cuadrasauge),
             "sinus" = list(tfun=tfun.sinus.utdc,utdc=utdc.sinus),
             "exponential" = list(tfun=tfun.exponential.utdc,
                 utdc=utdc.exponential))
    }else if(depcoefType=="ltdc"){
      switch(copula@family,
             "frechet" = list(tfun=tfun.frechet.ltdc,ltdc=ltdc.frechet),
             "cuadrasauge" = list(tfun=tfun.cuadrasauge.ltdc,
               ltdc=ltdc.cuadrasauge),
             "sinus" = list(tfun=tfun.sinus.ltdc,ltdc=ltdc.sinus),
             "exponential" = list(tfun=tfun.exponential.ltdc,
                 ltdc=ltdc.exponential))
    }else{
      stop("Wrong type of dependence coefficient")   
    }
  }
}
