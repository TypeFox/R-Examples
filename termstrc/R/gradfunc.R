### Gradient of Svensson grid loss function for bonds
grad_sv_bonds_grid <- function(beta, tau, m, cf, w, p) {
  .Call("grad_sv_bonds_gridCpp", beta, tau, m, cf, w, p)
    }

### Gradient of Svensson loss function for bonds
grad_sv_bonds <- function(beta,m,cf,w,p){

  emt1 <- exp(-m/beta[4])
  emt2 <- exp(-m/beta[6])
  oemt1 <- (1-emt1)
  t1emt1 <- beta[4]*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt2tm <- (-emt2 + beta[6]*(1 - emt2)/m)
  emt1t <- emt1/beta[4]
  emt2t <- emt2/beta[6]

  a <- exp((-beta[1] - beta[3]*emt1tm - beta[5]*emt2tm - (beta[2]*t1emt1)/m)*m/100)

 
  acf <- a*cf
  b <- -2*w*(p-cSums(acf,na.rm=TRUE))
  d <- acf/100
  dm <- d*m
  
  gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
  gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
  gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))
  gbeta5 <- sum(b*(-cSums(dm*emt2tm, na.rm=TRUE)))     
  gbeta4 <- sum(b*(-cSums( dm*(-beta[2]*emt1t + beta[2]*oemt1/m + beta[3]*(-emt1t+ oemt1/m - emt1t*m/beta[4])),na.rm=TRUE)))
  gbeta6 <- sum(b*(-cSums(dm*beta[5]*( -emt2t + (1-emt2)/m - emt2t*m/beta[6]),na.rm=TRUE)))

  c(gbeta1,gbeta2,gbeta3,gbeta4,gbeta5,gbeta6)
  
}

### Gradient of Nelson/Siegel grid loss function for bonds
grad_ns_bonds_grid <- function(beta, tau, m, cf, w, p){
   emt1 <- exp(-m/tau)
   oemt1 <- (1-emt1)
   t1emt1 <- tau*oemt1
   emt1tm <- (-emt1 + t1emt1/m)
   emt1t <- emt1/tau
 
   a <- exp(-((beta[1] + beta[3]*(-emt1 +t1emt1/m) + 
         (beta[2]*t1emt1)/m)*m)/100)

   acf <- a*cf
   b <- -2*w*(p-cSums(acf,na.rm=TRUE))
   d <- acf/100
   dm <- d*m

   gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
   gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
   gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))    

   c(gbeta1,gbeta2,gbeta3)
}

### Gradient of Nelson/Siegel grid loss function for bonds
grad_ns_bonds <- function(beta, m, cf, w, p){
   emt1 <- exp(-m/beta[4])
   oemt1 <- (1-emt1)
   t1emt1 <- beta[4]*oemt1
   emt1tm <- (-emt1 + t1emt1/m)
   emt1t <- emt1/beta[4]
 
   a <- exp(-((beta[1] + beta[3]*(-emt1 +t1emt1/m) + 
         (beta[2]*t1emt1)/m)*m)/100)

   acf <- a*cf
   b <- -2*w*(p-cSums(acf,na.rm=TRUE))
   d <- acf/100
   dm <- d*m

   gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
   gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
   gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))    
   gbeta4 <- sum(b*(-cSums( dm*(-beta[2]*emt1t + beta[2]*oemt1/m + beta[3]*(-emt1t+ oemt1/m - emt1t*m/beta[4])),na.rm=TRUE)))
   
   c(gbeta1,gbeta2,gbeta3,gbeta4)
}

### Gradient of Adjusted Svensson grid loss function for bonds
grad_asv_bonds_grid <- function(beta, tau, m, cf, w, p){
  emt1 <- exp(-m/tau[1])
  emt2 <- exp(-m/tau[2])
  oemt1 <- (1-emt1)
  t1emt1 <- tau[1]*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt2tm <- (-exp(-2*m/tau[2]) + (tau[2]*(1 -emt2))/m)
  emt1t <- emt1/tau[1]
 
  a <- exp(-(beta[1] + beta[3]*emt1tm + 
         beta[4]*emt2tm + 
         (beta[2]*t1emt1)/m)*m/100)

 
  acf <- a*cf
  b <- -2*w*(p-cSums(acf,na.rm=TRUE))
  d <- acf/100
  dm <- d*m

  gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
  gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
  gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))
  gbeta5 <- sum(b*(-cSums(dm*emt2tm, na.rm=TRUE)))     
  
  c(gbeta1, gbeta2, gbeta3, gbeta5)
  
}

### Gradient of Adjusted Svensson loss function for bonds
grad_asv_bonds <-  function(beta, m, cf, w, p){
  emt1 <- exp(-m/beta[4])
  emt2 <- exp(-m/beta[6])
  oemt1 <- (1-emt1)
  t1emt1 <- beta[4]*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt2tm <- (-exp(-2*m/beta[6]) + (beta[6]*(1 -emt2))/m)
  emt1t <- emt1/beta[4]
 
  a <- exp(-(beta[1] + beta[3]*emt1tm + 
         beta[5]*emt2tm + 
         (beta[2]*t1emt1)/m)*m/100)

 
  acf <- a*cf
  b <- -2*w*(p-cSums(acf,na.rm=TRUE))
  d <- acf/100
  dm <- d*m

  gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
  gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
  gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))
  gbeta4 <- sum(b*(-cSums( dm*(-beta[2]*emt1t + beta[2]*oemt1/m + beta[3]*(-emt1t+ oemt1/m - emt1t*m/beta[4])),na.rm=TRUE)))
  gbeta5 <- sum(b*(-cSums(dm*emt2tm, na.rm=TRUE)))     
  gbeta6 <- sum(b*(-cSums(dm*beta[5]*(-emt2/beta[6] + (1-emt2)/m - 2*exp(-2*m/beta[6])*m/beta[6]^2 ),na.rm=TRUE)))

  c(gbeta1, gbeta2, gbeta3, gbeta4, gbeta5, gbeta6)

}

### Gradient of Diebold/Li loss function for bonds
grad_dl_bonds <- function(beta, lambda, m, cf, w, p){
   tau <- 1/lambda 
   emt1 <- exp(-m/tau)
   oemt1 <- (1-emt1)
   t1emt1 <- tau*oemt1
   emt1tm <- (-emt1 + t1emt1/m)
   emt1t <- emt1/tau
 
   a <- exp(-((beta[1] + beta[3]*(-emt1 +t1emt1/m) + 
         (beta[2]*t1emt1)/m)*m)/100)

   acf <- a*cf
   b <- -2*w*(p-cSums(acf,na.rm=TRUE))
   d <- acf/100
   dm <- d*m

   gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
   gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
   gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))    

   c(gbeta1,gbeta2,gbeta3)
}

### was a Faster version of column sums
cSums <- function (x, na.rm = FALSE, dims = 1L) {

  colSums(x, na.rm = na.rm, dims = dims )

##    dn <- dim(x)
##    n <- prod(dn[1L:dims])
##    dn <- dn[-(1L:dims)]
##    z <-  .Internal(colSums(x, n, prod(dn), na.rm))
##    z
}




### Gradient of Nelson/Siegel loss function for yields
grad_ns <- function(beta, m, y)
      {
        c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
          (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum((-2*beta[4]*(1 - exp(-m/beta[4]))*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
           (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))/m),

          sum(-2*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m)*
           (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
           (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum(2*(beta[2]/(beta[4]*exp(m/beta[4])) - (beta[2]*(1 - exp(-m/beta[4])))/m - 
                 beta[3]*(-(1/(beta[4]*exp(m/beta[4]))) + (1 - exp(-m/beta[4]))/m - m/(beta[4]^2*exp(m/beta[4]))))
              *(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
                (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))
          )
      }


grad_ns_grid <- function(beta, tau, m, y)
      {
        c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
          (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y)),

          sum((-2*tau[1]*(1 - exp(-m/tau[1]))*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
           (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))/m),

          sum(-2*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m)*
           (-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
           (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))
          )
      }



### Gradient of Svensson loss function for yields
grad_sv <- function(beta, m, y)
      {
        c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
      beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
      (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum((-2*beta[4]*(1 - exp(-m/beta[4]))*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
        beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
        (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))/m),

          sum(-2*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m)*
    (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
      beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
      (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum(2*(beta[2]/(beta[4]*exp(m/beta[4])) - (beta[2]*(1 - exp(-m/beta[4])))/m - 
      beta[3]*(-(1/(beta[4]*exp(m/beta[4]))) + (1 - exp(-m/beta[4]))/m - m/(beta[4]^2*exp(m/beta[4]))))*
    (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
      beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
      (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum(-2*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m)*
    (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
      beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
      (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum(-2*beta[5]*(-(1/(beta[6]*exp(m/beta[6]))) + (1 - exp(-m/beta[6]))/m - m/(beta[6]^2*exp(m/beta[6])))*
    (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
      beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
      (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))
          )
      }


grad_sv_grid <- function(beta, tau, m, y)
      {
        .Call("grad_sv_gridCpp", beta, tau, m, y)
        
##         c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
##       beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
##       (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y)),

##           sum((-2*tau[1]*(1 - exp(-m/tau[1]))*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
##         beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
##         (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))/m),

##           sum(-2*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m)*
##     (-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
##       beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
##       (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y)),


##           sum(-2*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m)*
##     (-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
##       beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
##       (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))
##           )
      }


### Gradient of adjusted Svensson loss function for yields
grad_asv <- function(beta, m, y)
      {
        c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
          beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
          (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum((-2*beta[4]*(1 - exp(-m/beta[4]))*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
          beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
          (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))/m),
          
          sum(-2*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m)*
              (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
               beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
               (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
          
          sum(2*(beta[2]/(beta[4]*exp(m/beta[4])) - (beta[2]*(1 - exp(-m/beta[4])))/m - 
                 beta[3]*(-(1/(beta[4]*exp(m/beta[4]))) + (1 - exp(-m/beta[4]))/m - m/(beta[4]^2*exp(m/beta[4]))))*
              (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
               beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
               (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum(-2*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m)*
              (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
               beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
               (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
          
          sum(-2*beta[5]*(-(1/(beta[6]*exp(m/beta[6]))) + (1 - exp(-m/beta[6]))/m - 
                          (2*m)/(beta[6]^2*exp((2*m)/beta[6])))*
              (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
               beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
               (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))
          )
      }



grad_asv_grid <- function(beta, tau, m, y)
      {
        c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
          beta[4]*(-exp((-2*m)/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
          (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y)),

          sum((-2*tau[1]*(1 - exp(-m/tau[1]))*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
          beta[4]*(-exp((-2*m)/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
          (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))/m),
          
          sum(-2*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m)*
              (-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
               beta[4]*(-exp((-2*m)/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
               (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y)),
    

          sum(-2*(-exp((-2*m)/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m)*
              (-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
               beta[4]*(-exp((-2*m)/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
               (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))
          )
      }



grad_dl <- function(beta,lambda, m, y)
      {
        tau <- 1/lambda
        c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/tau) + (tau*(1 - exp(-m/tau)))/m) - 
          (beta[2]*tau*(1 - exp(-m/tau)))/m + y)),

          sum((-2*tau*(1 - exp(-m/tau))*(-beta[1] - beta[3]*(-exp(-m/tau) + (tau*(1 - exp(-m/tau)))/m) - 
           (beta[2]*tau*(1 - exp(-m/tau)))/m + y))/m),

          sum(-2*(-exp(-m/tau) + (tau*(1 - exp(-m/tau)))/m)*
           (-beta[1] - beta[3]*(-exp(-m/tau) + (tau*(1 - exp(-m/tau)))/m) - 
           (beta[2]*tau*(1 - exp(-m/tau)))/m + y))
          )
      }

