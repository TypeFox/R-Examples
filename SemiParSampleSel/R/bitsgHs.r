bitsgHs <- function(cop,margin,i1,F1,F2,f2,eta1,ph,teta,e2=0,sqv=0,ver=0){


# tofyDEL2

if(margin=='N'){


   if(cop=="C0"){

    lnF1 <- log(F1)
    lnF2 <- log(F2)
    teta1 <- teta+1	
    F1mt <- F1^(-teta) 			
    F2mt <- F2^(-teta) 			
    u <- F1mt+F2mt-1   			
    lnu <- log(u)	
    z <- F2^(-teta1)*u^(-teta1/teta)
    z <- ifelse(z<0.99999999,z,0.99999999)
    zt <- z/(1-z)
    h <- teta1*zt*f2*F2^(-teta1)*(u^(-1)-F2^teta)
    E <- f2*F2^(-teta1)*(F2^teta-teta/u) - h/z + e2/sqv 
    if(ver!=2){
       p <- teta1*zt*F1^{-teta1}*u^{-1}
       A <- -i1*p*(h/z+teta*f2*F2^(-teta1)*u^(-1))*ph
       P <- p*( F1^(-1)*ph*(teta1-F1mt*u^(-1)*(teta+teta1*(1-z)^(-1))) - eta1 )
    }
    if(ver!=1){
      C <- F1mt*lnF1+F2mt*lnF2
      Ct <- F1mt*lnF1^2+F2mt*lnF2^2
      b <- zt*(teta*lnF2-lnu/teta-C*teta1/u)
      B <- i1*h*( teta/teta1 -b/z + teta*(F2^teta*u-1)^(-1)*(lnF2-C/u) )
      if(ver!=2) h14 <- i1*teta*p*( 1/teta1 + C/u - b/(teta*z) - lnF1 )*ph
      h44 <- i1*(zt*( teta*lnF2 + lnu/teta + (1-teta)*C/u + teta*teta1*(Ct/u-C^2/u^2) ) - b^2/z)
    }
  }	
  
   if(cop=="C270"){

    F2b <- 1-F2
    lnF1 <- log(F1)
    lnF2b <- log(F2b)
    teta1 <- teta+1	
    F1mt <- F1^(-teta) 			
    F2bmt <- F2b^(-teta) 			
    u <- F1mt+F2bmt-1   			
    lnu <- log(u)	
    z <- F2b^(-teta1)*u^(-teta1/teta)
    z <- ifelse(z<0.99999999,z,0.99999999)
    zt <- z/(1-z)
    h <- - teta1*zt*f2*F2b^(-teta1)*(u^(-1)-F2b^teta)
    E <- -f2*F2b^(-teta1)*(F2b^teta-teta/u) - h/z + e2/sqv 
    if(ver!=2){
       p <- teta1*zt*F1^{-teta1}*u^{-1}
       A <- i1*p*( -h/z + teta*f2*F2b^(-teta1)*u^(-1) )*ph
       P <- p*( F1^(-1)*ph*(teta1-F1mt*u^(-1)*(teta+teta1*(1-z)^(-1))) - eta1 )
    }
    if(ver!=1){
      C <- F1mt*lnF1+F2bmt*lnF2b
      Ct <- F1mt*lnF1^2+F2bmt*lnF2b^2
      b <- zt*(teta*lnF2b-lnu/teta-C*teta1/u)
      B <- i1*h*( teta/teta1 -b/z + teta*(F2b^teta*u-1)^(-1)*(lnF2b-C/u) )
      if(ver!=2) h14 <- i1*teta*p*( 1/teta1 + C/u - b/(teta*z) - lnF1 )*ph
      h44 <- i1*(zt*( teta*lnF2b + lnu/teta + (1-teta)*C/u + teta*teta1*(Ct/u-C^2/u^2) ) - b^2/z)
    }
  }
	
  if(cop=="C90"){
  
    F1b <- 1-F1
    lnF1b <- log(F1b)
    lnF2 <- log(F2)
    teta1 <- teta+1
    u <- F1b^(-teta)+F2^(-teta)-1
    lnu <- log(u)
    z <- F2^(-teta1)*u^(-teta1/teta)
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    h <- teta1*f2*F2^(-teta1)*(F2^teta-u^(-1))
    E <- e2/sqv+f2*F2^(-teta1)*(F2^teta-teta/u)
    if(ver!=2){
       p <- teta1*F1b^(-teta1)*u^(-1)
       A <- -i1*teta*p*u^(-1)*f2*F2^(-teta1)*ph
       P <- -p*( ph*(teta1-F1b^(-teta)*teta/u)/F1b + eta1 )
    }
    if(ver!=1){    
       C <- F1b^(-teta)*lnF1b+F2^(-teta)*lnF2
       Ct <- F1b^(-teta)*lnF1b^2+F2^(-teta)*lnF2^2
       b <- lnu/teta+teta1*C*u^(-1)-teta*lnF2
       B <- i1*teta*f2*F2^(-teta1)*(F2^teta-teta1*u^(-1)*(1/teta1+C/u-lnF2))
       if(ver!=2) h14 <- i1*teta*p*(1/teta1+C/u-lnF1b)*ph
       h44 <- i1*((teta-1)*C/u-lnu/teta-teta*lnF2-teta*teta1*(Ct/u-C^2/u^2))
     }
     z <- 1-z

  }	
  
  if(cop=="C180"){
  
    F1b <- 1-F1
    F2b <- 1-F2
    lnF1b <- log(F1b)
    lnF2b <- log(F2b)
    teta1 <- teta+1
    u <- F1b^(-teta)+F2b^(-teta)-1
    lnu <- log(u)
    z <- F2b^(-teta1)*u^(-teta1/teta)
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    h <- - teta1*f2*F2b^(-teta1)*(F2b^teta-u^(-1))
    E <- e2/sqv - f2*F2b^(-teta1)*(F2b^teta-teta/u)
    if(ver!=2){
       p <- teta1*F1b^(-teta1)*u^(-1)
       A <- i1*teta*p*u^(-1)*f2*F2b^(-teta1)*ph
       P <- -p*( ph*(teta1-F1b^(-teta)*teta/u)/F1b + eta1 )
    }
    if(ver!=1){    
       C <- F1b^(-teta)*lnF1b+F2b^(-teta)*lnF2b
       Ct <- F1b^(-teta)*lnF1b^2+F2b^(-teta)*lnF2b^2
       b <- lnu/teta+teta1*C*u^(-1)-teta*lnF2b
       B <- - i1*teta*f2*F2b^(-teta1)*(F2b^teta-teta1*u^(-1)*(1/teta1+C/u-lnF2b))
       if(ver!=2) h14 <- i1*teta*p*(1/teta1+C/u-lnF1b)*ph
       h44 <- i1*((teta-1)*C/u-lnu/teta-teta*lnF2b-teta*teta1*(Ct/u-C^2/u^2))
     }
     z <- 1-z

  }	

  if(cop=="J0"){
	
    F1b <- 1-F1
    F2b <- 1-F2
    lnF1b <- log(F1b)
    lnF2b <- log(F2b)
    F1bt <- F1b^teta
    F2bt <- F2b^teta
    teta1 <- teta-1
    u <- F1bt+F2bt-F1bt*F2bt
    ut <- u^(1/teta-1)
    z <- (1-F1bt)*F2b^teta1*ut
    z <- ifelse(z<0.99999999,z,0.99999999)
    zt <- (z-1)^(-1)
    h <- teta1*z*zt*f2*F2b^(-1)*F1bt*u^(-1)
    E <- f2*F2b^(-1)*( teta*u^(-1)*F1bt-teta-1 ) - h/z + e2/sqv
    if(ver!=2){
       p <- -zt*(F1b*F2b)^teta1*ut*u^(-1)*(u+teta1)
       P <- ph*F1b^(-1)*( p*(teta1-p*F1b) - teta1*u^(-1)*(u-F2bt)*(2*p-z*zt*F1b^teta1*u^(-1)*(1-F2bt)) ) - p*eta1 
       A <- i1*(teta1*p*f2*F2b^(-1)*(1-z*zt*F1bt*u^(-1)) + F2bt*(2*teta+u-1)*u^(-1)*h*F1b^(-1))*ph
    }
    if(ver!=1){      
       C <- F1bt*lnF1b+F2bt*lnF2b-F1bt*F2bt*(lnF1b+lnF2b)
       Ct <- F1bt*lnF1b^2+F2bt*lnF2b^2-F1bt*F2bt*(lnF1b+lnF2b)^2
       b <- -teta1*z*zt*( log(u)/teta^2-lnF2b+teta1*teta^(-1)*C/u+F1bt*(1-F1bt)^(-1)*lnF1b )
       B <- i1*h*( 1-b/z+teta1*(lnF1b-C/u) )
       if(ver!=2)  h14 <- i1*p*ph*( teta1*( (1-teta1*C*u^(-1))*(u+teta1)^(-1)  + lnF1b*(1-F1bt)^(-1) ) - b/z )
       h44 <- i1*(  b*(1-b) + teta1*( (lnF2b-teta1*teta^(-1)*C*u^(-1)-log(u)/teta^2)*(b-teta1*z*zt*lnF1b*(1-F1bt)^(-1))
		+ b*lnF1b + teta1*z*zt*( -teta1*teta^(-1)*(Ct/u-C^2/u^2) - 2*teta^(-2)*(C/u-log(u)/teta) )
	      )  )
    }

  }

  if(cop=="J270"){
	
    F1b <- 1-F1
    lnF1b <- log(F1b)
    lnF2 <- log(F2)
    F1bt <- F1b^teta
    F2t <- F2^teta
    teta1 <- teta-1
    u <- F1bt+F2t-F1bt*F2t
    ut <- u^(1/teta-1)
    z <- (1-F1bt)*F2^teta1*ut
    z <- ifelse(z<0.99999999,z,0.99999999)
    zt <- (z-1)^(-1)
    h <- -teta1*z*zt*f2*F2^(-1)*F1bt*u^(-1)
    E <- -f2*F2^(-1)*( teta*u^(-1)*F1bt-teta-1 ) - h/z + e2/sqv
    if(ver!=2){
       p <- -zt*(F1b*F2)^teta1*ut*u^(-1)*(u+teta1)
       P <- ph*F1b^(-1)*( p*(teta1-p*F1b) - teta1*u^(-1)*(u-F2t)*(2*p-z*zt*F1b^teta1*u^(-1)*(1-F2t)) ) - p*eta1 
       A <- -i1*(teta1*p*f2*F2^(-1)*(1-z*zt*F1bt*u^(-1)) - F2t*(2*teta+u-1)*u^(-1)*h*F1b^(-1))*ph
    }
    if(ver!=1){      
       C <- F1bt*lnF1b+F2t*lnF2-F1bt*F2t*(lnF1b+lnF2)
       Ct <- F1bt*lnF1b^2+F2t*lnF2^2-F1bt*F2t*(lnF1b+lnF2)^2
       b <- -teta1*z*zt*( log(u)/teta^2-lnF2+teta1*teta^(-1)*C/u+F1bt*(1-F1bt)^(-1)*lnF1b )
       B <- i1*h*( 1-b/z+teta1*(lnF1b-C/u) )
       if(ver!=2)  h14 <- i1*p*ph*( teta1*( (1-teta1*C*u^(-1))*(u+teta1)^(-1)  + lnF1b*(1-F1bt)^(-1) ) - b/z )
       h44 <- i1*(  b*(1-b) + teta1*( (lnF2-teta1*teta^(-1)*C*u^(-1)-log(u)/teta^2)*(b-teta1*z*zt*lnF1b*(1-F1bt)^(-1))
		+ b*lnF1b + teta1*z*zt*( -teta1*teta^(-1)*(Ct/u-C^2/u^2) - 2*teta^(-2)*(C/u-log(u)/teta) )
	      )  )
    }

  }	
  
  if(cop=="J90"){

    F2b <- 1-F2
    lnF1 <- log(F1)
    lnF2b <- log(F2b)
    teta1 <- teta-1
    u <- F1^teta+F2b^teta-(F1*F2b)^teta
    ut <- u^(1/teta-1)
    z <- (1-F1^teta)*(F2b^teta1)*ut
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    h <- teta1*f2*F2b^(-1)*F1^teta*u^(-1)
    E <- f2*F2b^(-1)*( teta*u^(-1)*F1^teta-teta-1 ) + e2/sqv
    if(ver!=2){
       p <- F1^teta1*(1-F1^teta)^(-1)*(u+teta1)*u^(-1)
       A <- -i1*teta*teta1*u^(-2)*F1^teta1*F2b^teta1*f2*ph
       P <- p*( ph*F1^(-1)*( teta1*teta*(u-F2b^teta)*u^(-1)*(teta1+u)^(-1) - (teta1+F1^teta)*(1-F1^teta)^(-1) ) - eta1 )
    }
    if(ver!=1){ 
       C <- (F1^teta)*lnF1+(F2b^teta)*lnF2b-((F1*F2b)^teta)*(lnF1+lnF2b)
       b <- teta1*(lnF2b+(1/teta-1)*(C/u)-log(u)/(teta^2)-F1^teta*(1-F1^teta)^(-1)*lnF1)
       Ct <- (F1^teta)*(lnF1^2)+(F2b^teta)*(lnF2b^2)-((F1*F2b)^teta)*((lnF1+lnF2b)^2)
       B <- i1*h*(1+teta1*(lnF1-C/u))
       if(ver!=2)  h14 <- i1*teta1*(p*(lnF1*(1-F1^teta)^(-1)-C/u) + (1+C)*u^(-1)*F1^teta1*(1-F1^teta)^(-1))*ph 
       h44 <- i1*( b + teta1^2*(-teta1*teta^(-1)*(Ct/u-C^2/u^2) - (2/teta^2)*(C/u) + 
                              2*log(u)*teta^(-3) - lnF1^2*F1^teta*(1-F1^teta)^(-2) ) )
    }
    z <- 1-z

  }
  
  if(cop=="J180"){

    lnF1 <- log(F1)
    lnF2 <- log(F2)
    teta1 <- teta-1
    u <- F1^teta+F2^teta-(F1*F2)^teta
    ut <- u^(1/teta-1)
    z <- (1-F1^teta)*(F2^teta1)*ut
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    h <- - teta1*f2*F2^(-1)*F1^teta*u^(-1)
    E <- - f2*F2^(-1)*( teta*u^(-1)*F1^teta-teta-1 ) + e2/sqv
    if(ver!=2){
       p <- F1^teta1*(1-F1^teta)^(-1)*(u+teta1)*u^(-1)
       A <- i1*teta*teta1*u^(-2)*F1^teta1*F2^teta1*f2*ph
       P <- p*( ph*F1^(-1)*( teta1*teta*(u-F2^teta)*u^(-1)*(teta1+u)^(-1) - (teta1+F1^teta)*(1-F1^teta)^(-1) ) - eta1 )
        }
    if(ver!=1){ 
       C <- (F1^teta)*lnF1+(F2^teta)*lnF2-((F1*F2)^teta)*(lnF1+lnF2)
       b <- teta1*(lnF2+(1/teta-1)*(C/u)-log(u)/(teta^2)-F1^teta*(1-F1^teta)^(-1)*lnF1)
       Ct <- (F1^teta)*(lnF1^2)+(F2^teta)*(lnF2^2)-((F1*F2)^teta)*((lnF1+lnF2)^2)
       B <- i1*h*(1+teta1*(lnF1-C/u))
       if(ver!=2)  h14 <- i1*teta1*(p*(lnF1*(1-F1^teta)^(-1)-C/u) + (1+C)*u^(-1)*F1^teta1*(1-F1^teta)^(-1))*ph 
       h44 <- i1*( b + teta1^2*(-teta1*teta^(-1)*(Ct/u-C^2/u^2) - (2/teta^2)*(C/u) + 
                              2*log(u)*teta^(-3) - lnF1^2*F1^teta*(1-F1^teta)^(-2) ) )
    }
    z <- 1-z

  }

  if(cop=="FGM"){

    u <- 1-teta*F1*(1-2*F2)
    z <- 1-(1-F1)*(1-teta*F1*(1-2*F2))
    h <- -2*teta*F1*f2*u^(-1) 
    E <- u^(-1)*(e2*u*sqv^(-1) - h*u)
    if(ver!=2){
       p <- (1-F1)^(-1) + teta*(1-2*F2)/u
       P <- -( ph*( (1-F1)^(-2) + (teta*(1-2*F2)/u)^2 ) + eta1*((1-F1)^(-1) + teta*(1-2*F2)/u ) )
       A <- i1*2*teta*u^(-2)*f2*ph
    }
    if(ver!=1){ 
       b <- (teta^2-1)*F1*(1-2*F2)*u^(-1)
       B <- -i1*2*F1*f2*u^(-2)*(1-teta^2)
       if(ver!=2)  h14 <- i1*(1-2*F2)*u^(-2)*ph*(1-teta^2)
       h44 <- -i1*b*(2*teta+b)
    }

  }
  
  if(cop=="F"){
  
    u <- exp(teta*(F1+F2))-exp(teta*(1+F2))
    z <- 1-u/(u-exp(teta*(1+F1))+exp(teta))
    z <- ifelse(z<0.99999999,z,0.99999999)
    h <- (1-z)*u^(-1)*f2*teta*exp(teta)*(exp(teta*F1)-1)
    E <- teta*(1-z)*f2+e2/sqv
    if(ver!=2){
       p <- (1-exp(teta))*teta*(z-1)*exp(teta*(1+F1+F2))*u^(-2)
       A <- i1*(1-exp(teta))^(-1)*f2*p^2*u*(exp(-teta*F1)-exp(-teta))*ph
       P <- p*( ph*(exp(teta)-1)^(-1)*p*( exp(teta*F1)*(exp(teta*(F2-1))-1) - exp(teta*(1-F1))*(exp(teta*F2)-1) ) - eta1 )
    }
    if(ver!=1){ 
       b <- u^(-1)*( z*( u*(F1+F2-1) + (F1-1)*exp(teta*(1+F2)) ) + (1-z)*F1*exp(teta*(1+F1)) )
       B <- i1*f2*exp(teta)*(1-z)*u^(-1)*( (1+teta+teta*F1)*exp(teta*F1)-teta-1-
		teta*(exp(teta*F1)-1)*(1-z)*u^(-1)*(u*(F1+F2)+(F1-1)*exp(teta*(1+F2))-(1+F1)*exp(teta*(1+F1))+exp(teta) ) )
       if(ver!=2)  h14 <- i1*ph*teta^(-1)*(exp(teta)-1)^(-1)*p*( 
		     (1+teta)*exp(teta) - 1 - p*( (F1+F2-1)*exp(teta*(F1+F2-1)) - 2*F2*exp(teta*F2) -
		     F1*exp(teta*F1) + (1-F1+F2)*exp(teta*(1-F1+F2)) - (1-F1)*exp(teta*(1-F1)) + exp(teta) )
	            )
       h44 <- i1*(  b*( b - F1 - F2 - u^(-1)*exp(teta*(1+F2))*(F1-1) ) + 
	  (1-z)*u^(-1)*F1*(1+F1)*exp(teta*(1+F1)) +
	  (F1+F2-1)*(z*(F1+F2)-b)  +  u^(-1)*exp(teta*(1+F2))*(F1-1)*(z*(F1+2*F2)-b)
	)
    }
  }
  
  if(cop=="AMH"){

    u <- 1-teta*(1-F1)*(1-F2)
    z <- F1*(1-teta+teta*F1)*u^(-2);
    z <- ifelse(z<0.9999999,z,0.9999999)
    zt <- (z-1)^(-1)
    h <- 2*teta*z*zt*(1-F1)*f2*u^(-1)
    E <- teta*u^(-1)*(1-F1)*f2*(z-3)*zt+e2/sqv
    if(ver!=2){
       p <- zt*(2*teta*(z*u*(1-F2)-F1)+teta-1)*u^(-2)
       P <- ph*(2*teta*u^(-2)*((1-F2)*(2*p*u-z*zt*teta*(1-F2))+zt)-p^2)-p*eta1
       A <- i1*2*teta*zt*u^(-1)*f2*(z*u^(-1)-p*(1-F1))*ph
    }
    if(ver!=1){ 
       b <- zt*(1-teta^2)*(1-F1)*(2*z*(1-F2)*u^(-1)-F1*u^(-2))
       B <- i1*2*(1-F1)*f2*u^(-1)*zt*((1-teta^2)*z*u^(-1)-teta*b)
       if(ver!=2)  h14 <- i1*( (1-teta^2)*u^(-2)*( 2*teta^(-1)*p*u*(1-u) + z*zt*2*(2*u-1)*(1-F2) - (2*F1-1)*zt )  
                  + 2*teta*(1-F2)*b*u^(-1) - b*p
   	      )*ph
       h44 <- i1*( -b*(b+2*teta) + 2*(1-teta^2)*teta^(-1)*(1-u)*u^(-1)*(2*b-(1-teta^2)*teta^(-1)*(1-u)*u^(-1)*z*zt) )
    }
  }
  
  if(cop=="G0"){
  
    lnF1 <- -log(F1)
    lnF2 <- -log(F2)
    u <- lnF1^teta + lnF2^teta
    u <- ifelse(u>0.0000001,u,0.000001)
    ut <- u^(1/teta)
    z <- exp(-ut)*(F2^(-1))*ut*(u^(-1))*lnF2^(teta-1);	
    z <- ifelse(z<0.999999,z,0.999999)
    zt <- z/(z-1)
    p <- zt*lnF1^(teta-1)*F1^(-1)*u^(-1)*(1-teta-ut)
    h <- zt*F2^(-1)*f2*( lnF2^(-1)*(-F1*lnF1*p*zt^(-1)-ut) + 1 )
    E <- e2/sqv - h/z + f2*F2^(-1)*(
	   1 + (teta-1)*lnF2^(-2)*zt*f2*F2^(-1)*h^(-1)*(
	    lnF2^teta*( (1-teta-ut)*u^(-1) + (teta+ut)*u^(-2)*lnF2^teta ) - 1  )   )
    if(ver!=2){
       P <- p*( F1^(-1)*ph*( 1+lnF1^(-1)*(teta-1-lnF1^teta*u^(-1)*(teta+ut*(1-teta-ut)^(-1))) - p*F1*z^(-1) ) - eta1 )
       A <- i1*ph*( (teta-1)*f2*F2^(-1)*p*(teta+ut)*u^(-1)*lnF2^(teta-1)*(1-teta-ut)^(-1) - h*p*z^(-1) )
    }
    if(ver!=1){
       C <- lnF1^teta*log(lnF1)+lnF2^teta*log(lnF2)
       Ct <- lnF1^teta*(log(lnF1))^2+lnF2^teta*(log(lnF2))^2
       b <- zt*(teta-1)*( teta^(-1)*(C/u-log(u)/teta)*(1-ut) + log(lnF2) - C/u )
       B <- i1*(  (teta-1)*zt*F2^(-1)*f2*( lnF2^(teta-1)*u^(-1)*
   	   (log(lnF2)*(1-teta-ut)+(teta-1)*teta^(-1)*C*u^(-1)*(teta+ut)-1+ut*log(u)*teta^(-2)) - 
   	   b*h*(teta-1)^(-1)*(z-1)*z^(-2)*f2^(-1)*F2 + lnF2^(-1) )  )
       if(ver!=2)  h14 <- i1*p*(
	   (teta-1)*(
	    log(lnF1)-(1-teta-ut)^(-1)*((1-teta)*teta^(-1)*C*u^(-1)*(teta+ut)+1-ut*log(u)*teta^(-2))
	   ) - b*z^(-1) )*ph
       h44 <- i1*( (teta-1)^2*zt*( (1-ut)*teta^(-1)*(Ct/u-C^2/u^2-2*teta^(-1)*(C/u-log(u)/teta)) -
		  ut*teta^(-2)*(C/u-log(u)/teta)^2 - Ct/u + C^2/u^2 ) +  b - b^2/z )
    }
  }
  

  if(cop=="G270"){
  
    lnF1 <- -log(F1)
    F2b <- 1-F2
    lnF2b <- -log(F2b)
    u <- lnF1^teta + lnF2b^teta
    u <- ifelse(u>0.0000001,u,0.000001)
    ut <- u^(1/teta)
    z <- exp(-ut)*(F2b^(-1))*ut*(u^(-1))*lnF2b^(teta-1);	
    z <- ifelse(z<0.999999,z,0.999999)
    zt <- z/(z-1)
    p <- zt*lnF1^(teta-1)*F1^(-1)*u^(-1)*(1-teta-ut)
    h <- - zt*F2b^(-1)*f2*( lnF2b^(-1)*(-F1*lnF1*p*zt^(-1)-ut) + 1 )
    E <-  e2/sqv - h/z - f2*F2b^(-1)*(
	   1 - (teta-1)*lnF2b^(-2)*zt*f2*F2b^(-1)*h^(-1)*(
	    lnF2b^teta*( (1-teta-ut)*u^(-1) + (teta+ut)*u^(-2)*lnF2b^teta ) - 1  ) )
    if(ver!=2){
       P <- p*( F1^(-1)*ph*( 1+lnF1^(-1)*(teta-1-lnF1^teta*u^(-1)*(teta+ut*(1-teta-ut)^(-1))) - p*F1*z^(-1) ) - eta1 )
       A <- - i1*ph*( (teta-1)*f2*F2b^(-1)*p*(teta+ut)*u^(-1)*lnF2b^(teta-1)*(1-teta-ut)^(-1) + h*p*z^(-1) )
    }
    if(ver!=1){
       C <- lnF1^teta*log(lnF1)+lnF2b^teta*log(lnF2b)
       Ct <- lnF1^teta*(log(lnF1))^2+lnF2b^teta*(log(lnF2b))^2
       b <- zt*(teta-1)*( teta^(-1)*(C/u-log(u)/teta)*(1-ut) + log(lnF2b) - C/u )
       B <- - i1*(  (teta-1)*zt*F2b^(-1)*f2*( lnF2b^(teta-1)*u^(-1)*
   	   (log(lnF2b)*(1-teta-ut)+(teta-1)*teta^(-1)*C*u^(-1)*(teta+ut)-1+ut*log(u)*teta^(-2)) + 
   	   b*h*(teta-1)^(-1)*(z-1)*z^(-2)*f2^(-1)*F2b + lnF2b^(-1) )  )
       if(ver!=2)    h14 <- i1*p*(
	   (teta-1)*(
	    log(lnF1)-(1-teta-ut)^(-1)*((1-teta)*teta^(-1)*C*u^(-1)*(teta+ut)+1-ut*log(u)*teta^(-2))
	   ) - b*z^(-1) )*ph
       h44 <- i1*( (teta-1)^2*zt*( (1-ut)*teta^(-1)*(Ct/u-C^2/u^2-2*teta^(-1)*(C/u-log(u)/teta)) -
		  ut*teta^(-2)*(C/u-log(u)/teta)^2 - Ct/u + C^2/u^2 ) +  b - b^2/z )
    }
  }

  if(cop=="G90"){
  
    F1b <- 1-F1
    lnF1b <- -log(F1b)
    lnF2 <- -log(F2)
    u <- lnF1b^teta + lnF2^teta
    u <- ifelse(u>0.0000001,u,0.0000001)
    ut <- u^(1/teta)
    z <- exp(-ut)*(F2^(-1))*ut*(u^(-1))*lnF2^(teta-1)	
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    p <- -lnF1b^(teta-1)*F1b^(-1)*u^(-1)*(1-teta-ut)
    h <- F2^(-1)*f2*(lnF2^(-1)*(F1b*lnF1b*p-ut)+1)
    E <- e2/sqv + F2^(-1)*f2*(
         1 + (teta-1)*F2^(-1)*f2*(lnF2)^(-2)*h^(-1)*(lnF2^teta*((1-teta-ut)/u+lnF2^teta*(teta+ut)*u^(-2))-1)
        )
    if(ver!=2){
       A <- -i1*( (teta-1)*(teta+ut)*u^(-2)*(lnF2*lnF1b)^(teta-1)*F1b^(-1)*F2^(-1)*f2 )*ph
       P <- -p*(  F1b^(-1)*ph*( (teta-1)*lnF1b^(-1) + 1 + (teta-1)*lnF1b^(teta-1)*u^(-1)*(ut+teta)*(1-teta-ut)^(-1))
		+eta1  )
    }
    if(ver!=1){
       C <- (lnF1b^teta)*log(lnF1b)+(lnF2^teta)*log(lnF2)
       b <- (teta-1)*(teta^(-1)*(C/u-log(u)/teta)*(1-ut)+log(lnF2)-C/u)
       Ct <- (lnF1b^teta)*(log(lnF1b))^2+(lnF2^teta)*(log(lnF2))^2
       B <- i1*(teta-1)*f2*F2^(-1)*lnF2^(-1)*(
         1+lnF2^teta*((1-teta-ut)*u^(-1)*(log(lnF2)-C/u) - 1/u -teta^(-1)*ut*u^(-1)*(C/u-log(u)/teta))  )
       if(ver!=2)  h14 <- i1*(teta-1)*p*( log(lnF1b)-C/u-(1-teta-ut)^(-1)*(1+ut*teta^(-1)*(C/u-log(u)*teta^(-1))) )*ph
       h44 <- i1*( b +(teta-1)^2*(
                 (1-teta-ut)*teta^(-1)*(Ct/u-C^2/u^2) - teta^(-2)*(C/u-log(u)/teta)*(2*(1-ut)+ut*(C/u-log(u)/teta)) 
		) )
    }
    z <- 1-z

  }

  if(cop=="G180"){
  
    F1b <- 1-F1
    F2b <- 1-F2
    lnF1b <- -log(F1b)
    lnF2b <- -log(F2b)
    u <- lnF1b^teta + lnF2b^teta
    u <- ifelse(u>0.0000001,u,0.0000001)
    ut <- u^(1/teta)
    z <- exp(-ut)*(F2b^(-1))*ut*(u^(-1))*lnF2b^(teta-1)	
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    p <- -lnF1b^(teta-1)*F1b^(-1)*u^(-1)*(1-teta-ut)
    h <- - F2b^(-1)*f2*(lnF2b^(-1)*(F1b*lnF1b*p-ut)+1)
    E <- e2/sqv - F2b^(-1)*f2*(
         1 - (teta-1)*F2b^(-1)*f2*(lnF2b)^(-2)*h^(-1)*(lnF2b^teta*((1-teta-ut)/u+lnF2b^teta*(teta+ut)*u^(-2))-1)
        )
    if(ver!=2){
       A <-  i1*( (teta-1)*(teta+ut)*u^(-2)*(lnF2b*lnF1b)^(teta-1)*F1b^(-1)*F2b^(-1)*f2 )*ph
       P <- -p*(  F1b^(-1)*ph*( (teta-1)*lnF1b^(-1) + 1 + (teta-1)*lnF1b^(teta-1)*u^(-1)*(ut+teta)*(1-teta-ut)^(-1))
		+eta1  )
    }
    if(ver!=1){
       C <- (lnF1b^teta)*log(lnF1b)+(lnF2b^teta)*log(lnF2b)
       b <- (teta-1)*(teta^(-1)*(C/u-log(u)/teta)*(1-ut)+log(lnF2b)-C/u)
       Ct <- (lnF1b^teta)*(log(lnF1b))^2+(lnF2b^teta)*(log(lnF2b))^2
       B <- - i1*(teta-1)*f2*F2b^(-1)*lnF2b^(-1)*(
         1+lnF2b^teta*((1-teta-ut)*u^(-1)*(log(lnF2b)-C/u) - 1/u -teta^(-1)*ut*u^(-1)*(C/u-log(u)/teta))  )
       if(ver!=2)  h14 <- i1*(teta-1)*p*( log(lnF1b)-C/u-(1-teta-ut)^(-1)*(1+ut*teta^(-1)*(C/u-log(u)*teta^(-1))) )*ph
       h44 <- i1*( b +(teta-1)^2*(
                 (1-teta-ut)*teta^(-1)*(Ct/u-C^2/u^2) - teta^(-2)*(C/u-log(u)/teta)*(2*(1-ut)+ut*(C/u-log(u)/teta)) 
		) )
    }
    z <- 1-z

  }
}    

#################################################################################################################

if(margin=='G'){


  if(cop=="N"){
  
    q <- qnorm(F2)
    tetat <- 1/sqrt(1-teta^2)
    d  <- (eta1+teta*q)*tetat
    l1 <- ifelse(eta1<37.5,dnorm(-eta1)/pnorm(-eta1),eta1)
    l2 <- ifelse(d>-37.5,dnorm(d)/pnorm(d),-d)
    h <- -teta*tetat*l2*(dnorm(q))^(-1)  
    E <- (teta*tetat*d-q)*(dnorm(q))^(-1)-h
    if(ver!=2){
       A <- -i1*tetat*h*(d+l2)
       p <- l2*tetat*ph^(-1)
       P <- -ph^(-1)*tetat^2*l2*(d+l2)
    }
    if(ver!=1){
       b  <- l2*(teta*eta1+q)*tetat
       B <- i1*h*( -b*(l2^(-1)*d+1) + teta^(-1) )
       h14 <- i1*tetat*( -b*(d+l2) + l2*teta )
       h44 <- i1*( b*(-b*(l2^(-1)*d+1)+teta) + l2*eta1/tetat )
    }
    z <- 1-pnorm(d)

  }	
  
  if(cop=="C0"){
  
    lnF1 <- log(F1)
    lnF2 <- log(F2)
    teta1 <- teta+1
    u <- F1^(-teta)+F2^(-teta)-1
    z <- F2^(-teta1)*u^(-teta1/teta);  
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    zt <- z/(1-z)
    h <- -zt*teta1*F2^(-teta1)*(F2^teta-u^(-1))
    E <- F2^(-1)-F2^(-teta1)*teta*u^(-1)-h/z
    if(ver!=2){
       p <- teta1*F1^(-teta1)*u^(-1)*zt
       P <- i1*p*( ph*(teta1-F1^(-teta)*(teta+teta1/(1-z))/u)/F1 - eta1 )
       A <- -i1*p*( F2^(-1) - E )*ph 
    }
    if(ver!=1){
       C <- F1^(-teta)*lnF1+F2^(-teta)*lnF2
       b <- teta*zt*( lnF2 - log(u)/teta^2 - C*teta1*teta^(-1)*u^(-1) )
       B <- -teta*zt*i1*F2^(-teta1)*( (1-teta1*b*teta^(-1)*z^(-1))*(F2^teta-u^(-1)) + teta1*(lnF2/u-C*u^(-2)) )
       if(ver!=2)  h14 <- i1*ph*p*( teta1^(-1) + C/u - lnF1 - b*teta^(-1)*z^(-1) )*teta
       h44 <- i1*zt*( z^(-2)*b^2*(z-1) + teta*lnF2 + log(u)/teta + (1-teta)*C*u^(-1) +
	   teta*teta1*((F1^(-teta)*lnF1^2+F2^(-teta)*lnF2^2)*u^(-1)-C^2*u^(-2)) )
    }

  }

  if(cop=="C270"){
  
    F2b <- 1-F2
    lnF1 <- log(F1)
    lnF2b <- log(F2b)
    teta1 <- teta+1
    u <- F1^(-teta)+F2b^(-teta)-1
    z <- F2b^(-teta1)*u^(-teta1/teta);  
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    zt <- z/(1-z)
    h <- zt*teta1*F2b^(-teta1)*(F2b^teta-u^(-1))
    E <- - F2b^(-1) + F2b^(-teta1)*teta*u^(-1) - h/z
    if(ver!=2){
       p <- teta1*F1^(-teta1)*u^(-1)*zt
       P <- i1*p*( ph*(teta1-F1^(-teta)*(teta+teta1/(1-z))/u)/F1 - eta1 )
       A <- i1*p*( F2b^(-1) - F2b^(-1) + F2b^(-teta1)*teta*u^(-1) - h/z )*ph
    }
    if(ver!=1){
       C <- F1^(-teta)*lnF1+F2b^(-teta)*lnF2b
       b <- teta*zt*( lnF2b - log(u)/teta^2 - C*teta1*teta^(-1)*u^(-1) )
       B <- teta*zt*i1*F2b^(-teta1)*( (1-teta1*b*teta^(-1)*z^(-1))*(F2b^teta-u^(-1)) + teta1*(lnF2b/u-C*u^(-2)) )
       if(ver!=2)  h14 <- i1*ph*p*( teta1^(-1) + C/u - lnF1 - b*teta^(-1)*z^(-1) )*teta
       h44 <- i1*zt*( z^(-2)*b^2*(z-1) + teta*lnF2b + log(u)/teta + (1-teta)*C*u^(-1) +
	   teta*teta1*((F1^(-teta)*lnF1^2+F2b^(-teta)*lnF2b^2)*u^(-1)-C^2*u^(-2)) )
    }

  }	

  if(cop=="C90"){
  
    F1b <- 1-F1
    lnF1b <- log(F1b)
    lnF2 <- log(F2)
    teta1 <- teta+1
    u <- F1b^(-teta)+F2^(-teta)-1
    z <- F2^(-teta1)*u^(-teta1/teta)  
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    h <- teta1*F2^(-teta1)*(F2^teta-u^(-1))
    E <- F2^(-teta1)*(F2^(teta)-teta/u)
    if(ver!=2){
       p <- F1b^(-teta1)*u^(-1)*teta1
       P <- -p*(ph*(teta1-F1b^(-teta)*teta/u)/F1b+eta1)
       A <- -i1*teta*p*u^(-1)*F2^(-teta1)*ph  
    }
    if(ver!=1){
       C <- F1b^(-teta)*lnF1b+F2^(-teta)*lnF2
       Ct <- F1b^(-teta)*lnF1b^2+F2^(-teta)*lnF2^2
       b <- log(u)/teta+teta1*C*u^(-1)-teta*lnF2
       B <- i1*teta*( h/teta1+teta1*u^(-1)*F2^(-teta1)*(lnF2-C/u) )
       if(ver!=2)  h14 <- i1*teta*p*(1/teta1+C/u-lnF1b)*ph
       h44 <- i1*((teta-1)*C/u-log(u)/teta-teta*lnF2-teta*teta1*(Ct/u-C^2/u^2))
    }
    z <- 1-z

  }	

  if(cop=="C180"){
  
    F1b <- 1-F1
    F2b <- 1-F2
    lnF1b <- log(F1b)
    lnF2b <- log(F2b)
    teta1 <- teta+1
    u <- F1b^(-teta)+F2b^(-teta)-1
    z <- F2b^(-teta1)*u^(-teta1/teta)  
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    h <- - teta1*F2b^(-teta1)*(F2b^teta-u^(-1))
    E <- - F2b^(-teta1)*(F2b^(teta)-teta/u)
    if(ver!=2){
       p <- F1b^(-teta1)*u^(-1)*teta1
       P <- -p*(ph*(teta1-F1b^(-teta)*teta/u)/F1b+eta1)
       A <- i1*teta*p*u^(-1)*F2b^(-teta1)*ph  
    }
    if(ver!=1){
       C <- F1b^(-teta)*lnF1b+F2b^(-teta)*lnF2b
       Ct <- F1b^(-teta)*lnF1b^2+F2b^(-teta)*lnF2b^2
       b <- log(u)/teta+teta1*C*u^(-1)-teta*lnF2b
       B <- - i1*teta*( - h/teta1 + teta1*u^(-1)*F2b^(-teta1)*(lnF2b-C/u) )
       if(ver!=2)  h14 <- i1*teta*p*(1/teta1+C/u-lnF1b)*ph
       h44 <- i1*((teta-1)*C/u-log(u)/teta-teta*lnF2b-teta*teta1*(Ct/u-C^2/u^2))
    }
    z <- 1-z

  }	
  
  if(cop=="J0"){
  
    F1b <- 1-F1
    F2b <- 1-F2
    lnF1b <- log(F1b)
    lnF2b <- log(F2b)
    teta1 <- teta-1
    u <- F1b^teta+F2b^teta-(F1b*F2b)^teta
    ut <- u^(1/teta-1)
    z <- (1-F1b^teta)*(F2b^teta1)*ut;  		
    z <- ifelse(z<0.99999999, z,0.99999999)
    z <- ifelse(z>0.000000001,z,0.000000001)
    zt <- (z-1)^(-1)
    h <- z*zt*teta1*F2b^(-1)*F1b^teta*u^(-1)
    E <- -F2b^(-1)*( 1 + teta*(1-u^(-1)*F1b^teta) + h*F2b*z^(-1) ) 

    if(ver!=2){
       p <- F2b^teta1*ut*(u+teta1)*teta^(-1)*u^(-1)
       P <- -teta*zt*F1b^(teta-2)*(
   	   ph*( p*(teta1+teta*p*F1b^teta*zt) - teta1*u^(-1)*(u-F2b^teta)*(2*p+z*(teta*u)^(-1)*(1-F2b^teta)) )
	   -F1b*p*eta1
	   )
       A <- teta*i1*h*F1b^teta1*( p*zt*z^(-1) + F2b^teta*F1b^(-teta)*u^(-1) )*ph
    }
    if(ver!=1){
       C <- (F1b^teta)*lnF1b+(F2b^teta)*lnF2b-((F1b*F2b)^teta)*(lnF1b+lnF2b)
       b <- lnF2b-(teta1/teta)*(C/u)-log(u)/(teta^2)
       B <- z*b-(F2b^teta1)*ut*(F1b^teta)*lnF1b
       Ct <- (F1b^teta)*(lnF1b^2)+(F2b^teta)*(lnF2b^2)-((F1b*F2b)^teta)*((lnF1b+lnF2b)^2)
       Bt <- i1*h*(1-teta1*(C/u-lnF1b+B*z^(-1)*zt))
       if(ver!=2)  h14 <- -i1*teta1*zt*F1b^teta1*p*(
		        1 + teta*(lnF1b-B*zt+b+((1-u)/teta-teta1*C*u^(-1))/(u+teta1))
	              )*ph
       h44 <- i1*zt*teta1*( B - teta1*(
	       zt*B^2-B*b+z*(teta1*teta^(-1)*(Ct/u-C^2/u^2) + 2*teta^(-2)*(C/u-log(u)/teta)) + (z*b-B)*(lnF1b+b)
	      ) )

       b <- zt*B*teta1
       B <- Bt
    }

    if(ver!=2) p <- -teta*zt*F1b^teta1*p


  }	
  
  if(cop=="J270"){
  
    F1b <- 1-F1
    lnF1b <- log(F1b)
    lnF2 <- log(F2)
    teta1 <- teta-1
    u <- F1b^teta+F2^teta-(F1b*F2)^teta
    ut <- u^(1/teta-1)
    z <- (1-F1b^teta)*(F2^teta1)*ut;  		
    z <- ifelse(z<0.99999999, z,0.99999999)
    z <- ifelse(z>0.000000001,z,0.000000001)
    zt <- (z-1)^(-1)
    h <- - z*zt*teta1*F2^(-1)*F1b^teta*u^(-1)
    E <- F2^(-1)*( 1 + teta*(1-u^(-1)*F1b^teta) - h*F2*z^(-1) ) 

    if(ver!=2){
       p <- F2^teta1*ut*(u+teta1)*teta^(-1)*u^(-1)
       P <- -teta*zt*F1b^(teta-2)*(
   	   ph*( p*(teta1+teta*p*F1b^teta*zt) - teta1*u^(-1)*(u-F2^teta)*(2*p+z*(teta*u)^(-1)*(1-F2^teta)) )
	   -F1b*p*eta1
	   )
       A <- teta*i1*h*F1b^teta1*( p*zt*z^(-1) + F2^teta*F1b^(-teta)*u^(-1) )*ph
    }
    if(ver!=1){
       C <- (F1b^teta)*lnF1b+(F2^teta)*lnF2-((F1b*F2)^teta)*(lnF1b+lnF2)
       b <- lnF2-(teta1/teta)*(C/u)-log(u)/(teta^2)
       B <- z*b - (F2^teta1)*ut*(F1b^teta)*lnF1b
       Ct <- (F1b^teta)*(lnF1b^2)+(F2^teta)*(lnF2^2)-((F1b*F2)^teta)*((lnF1b+lnF2)^2)
       Bt <- i1*h*(1-teta1*(C/u-lnF1b+B*z^(-1)*zt))
       if(ver!=2)      h14 <- -i1*teta1*zt*F1b^teta1*p*(
		        1 + teta*(lnF1b-B*zt+b+((1-u)/teta-teta1*C*u^(-1))/(u+teta1))
	              )*ph
       h44 <- i1*zt*teta1*( B - teta1*(
	       zt*B^2-B*b+z*(teta1*teta^(-1)*(Ct/u-C^2/u^2) + 2*teta^(-2)*(C/u-log(u)/teta)) + (z*b-B)*(lnF1b+b)
	      ) )
       b <- zt*B*teta1
       B <- Bt
    }

    if(ver!=2) p <- -teta*zt*F1b^teta1*p


  }

  if(cop=="J90"){
  
    F2b <- 1-F2
    lnF1 <- log(F1)
    lnF2b <- log(F2b)
    teta1 <- teta-1
    u <- F1^teta+F2b^teta-(F1*F2b)^teta
    ut <- u^(1/teta-1)
    z <- (1-F1^teta)*(F2b^teta1)*ut		
    z <- ifelse(z<0.99999999, z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.000000001)
    h <- teta1*F2b^(-1)*F1^teta*u^(-1)
    E <- F2b^(-1)*( - 1 - teta + teta*u^(-1)*F1^teta )
    if(ver!=2){
       p <- F1^teta1*(1-F1^teta)^(-1)*(u+teta1)*u^(-1)
       P <- p*( ph*F1^(-1)*( teta1*teta*(u-F2b^teta)*u^(-1)*(teta1+u)^(-1) - 
             (teta1+F1^teta)*(1-F1^teta)^(-1) ) - eta1 )
       A <- -i1*teta*teta1*u^(-2)*F1^teta1*F2b^teta1*ph
    }
    if(ver!=1){
       C <- F1^teta*lnF1 + F2b^teta*lnF2b - (F1*F2b)^teta*(lnF1+lnF2b)
       b <- teta1*( lnF2b + (1/teta-1)*(C/u) - log(u)/(teta^2) - F1^teta*(1-F1^teta)^(-1)*lnF1 )
       Ct <- F1^teta*lnF1^2 + F2b^teta*lnF2b^2 - (F1*F2b)^teta*(lnF1+lnF2b)^2
       B <- h*( 1 + teta1*(lnF1-C/u) )
       if(ver!=2)  h14 <- i1*teta1*(p*(lnF1*(1-F1^teta)^(-1)-C/u) + (1+C)*u^(-1)*F1^teta1*(1-F1^teta)^(-1))*ph  
       h44 <- i1*( b + teta1^2*(-teta1*teta^(-1)*(Ct/u-C^2/u^2) - (2/teta^2)*(C/u) + 
                              2*log(u)*teta^(-3) - lnF1^2*F1^teta*(1-F1^teta)^(-2) ) )
    }

    z <- 1-z

  }
  
  if(cop=="J180"){
  
    lnF1 <- log(F1)
    lnF2 <- log(F2)
    teta1 <- teta-1
    u <- F1^teta+F2^teta-(F1*F2)^teta
    ut <- u^(1/teta-1)
    z <- (1-F1^teta)*(F2^teta1)*ut		
    z <- ifelse(z<0.99999999, z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.000000001)
    h <- - teta1*F2^(-1)*F1^teta*u^(-1)
    E <- - F2^(-1)*( - 1 - teta + teta*u^(-1)*F1^teta )
    if(ver!=2){
       p <- F1^teta1*(1-F1^teta)^(-1)*(u+teta1)*u^(-1)
       P <- p*( ph*F1^(-1)*( teta1*teta*(u-F2^teta)*u^(-1)*(teta1+u)^(-1) - 
             (teta1+F1^teta)*(1-F1^teta)^(-1) ) - eta1 )
       A <- i1*teta*teta1*u^(-2)*F1^teta1*F2^teta1*ph
    }
    if(ver!=1){
       C <- F1^teta*lnF1 + F2^teta*lnF2 - (F1*F2)^teta*(lnF1+lnF2)
       b <- teta1*( lnF2 + (1/teta-1)*(C/u) - log(u)/(teta^2) - F1^teta*(1-F1^teta)^(-1)*lnF1 )
       Ct <- F1^teta*lnF1^2 + F2^teta*lnF2^2 - (F1*F2)^teta*(lnF1+lnF2)^2
       B <- h*( 1 + teta1*(lnF1-C/u) )
       if(ver!=2)  h14 <- i1*teta1*(p*(lnF1*(1-F1^teta)^(-1)-C/u) + (1+C)*u^(-1)*F1^teta1*(1-F1^teta)^(-1))*ph  
       h44 <- i1*( b + teta1^2*(-teta1*teta^(-1)*(Ct/u-C^2/u^2) - (2/teta^2)*(C/u) + 
                              2*log(u)*teta^(-3) - lnF1^2*F1^teta*(1-F1^teta)^(-2) ) )
    }

    z <- 1-z

  }

  if(cop=="FGM"){
  
    u <- 1-teta*F1*(1-2*F2)
    z <- 1-(1-F1)*u
    h <- -2*teta*F1*u^(-1)
    E <- -h
    if(ver!=2){
       p <- (1-F1)^(-1) + teta*(1-2*F2)*u^(-1) 
       P <- -( ph*( (1-F1)^(-2) + (teta*(1-2*F2)/u)^2 ) + eta1*( (1-F1)^(-1) + teta*(1-2*F2)/u ) )
       A <- i1*2*teta*u^(-2)*ph
    }
    if(ver!=1){
       b <- (1-u^(-1))*(1-teta^2)*teta^(-1)
       B <- i1*2*F1*u^(-2)*(teta^2-1)
       if(ver!=2)  h14 <- i1*(1-2*F2)*u^(-2)*ph*(1-teta^2)
       h44 <- i1*(u-1)*u^(-1)*( (1-teta^2)*u^(-1) - 1 - teta^2 )*(1-teta^2)*teta^(-2)
    }

 }
  
  if(cop=="F"){
  
    u <- exp(teta*(F1+F2)) - exp(teta*(1+F2))
    z <- u/(u-exp(teta*(1+F1))+exp(teta))
    z <- ifelse(z>0.000001,z,0.000001)
    h <- z*u^(-1)*teta*exp(teta)*(exp(teta*F1)-1)
    E <- z*teta
    if(ver!=2){
       p <- -(1-exp(teta))*teta*z*exp(teta*(1+F1+F2))*u^(-2)
       P <- p*( -ph*p*(1-exp(teta))^(-1)*(exp(teta*F1)*(exp(teta*(F2-1))-1)-exp(teta*(1-F1))*(exp(teta*F2)-1)) - eta1 )
       A <- i1*(1-exp(teta))^(-1)*p^2*u*(exp(-teta*F1)-exp(-teta))*ph
    }
    if(ver!=1){
       b <- z*u^(-1)*( (1-z)*z^(-1)*(u*(F1+F2-1)+(F1-1)*exp(teta*(1+F2))) + F1*exp(teta*(1+F1)) )
       B <- i1*exp(teta)*z*u^(-1)*( (1+teta+teta*F1)*exp(teta*F1) - teta - 1 -
		teta*(exp(teta*F1)-1)*z*u^(-1)*(u*(F1+F2)+(F1-1)*exp(teta*(1+F2))-(1+F1)*exp(teta*(1+F1))+exp(teta) ) )
       if(ver!=2)  h14 <- -i1*teta^(-1)*(1-exp(teta))^(-1)*p*( (1+teta)*exp(teta) - 1 - p*(
		(F1+F2-1)*exp(teta*(F1+F2-1))-2*F2*exp(teta*F2)-F1*exp(teta*F1)+
		(1-F1+F2)*exp(teta*(1-F1+F2))-(1-F1)*exp(teta*(1-F1))+exp(teta)
   	      ) )*ph
       h44 <- i1*( b*( b - u^(-1)*(u*(F1+F2)+exp(teta*(1+F2))*(F1-1)) )
		+ z*u^(-1)*F1*(1+F1)*exp(teta*(1+F1)) + (F1+F2-1)*((1-z)*(F1+F2)-b)
		+ u^(-1)*exp(teta*(1+F2))*(F1-1)*( (1-z)*(F1+2*F2)-b )
		)
    }
    z <- 1-z

  }
  
  if(cop=="AMH"){
  
    u <- 1-teta*(1-F1)*(1-F2)
    u <- ifelse(u>0.0000001,u,0.0000001)
    z <- F1*(1-teta+teta*F1)*u^(-2);
    z <- ifelse(z<0.9999999,z,0.9999999)
    zt <- (z-1)^(-1)
    h <- 2*teta*z*zt*(1-F1)*u^(-1)
    if(ver!=2){
       p <- zt*( 2*teta*(z*u*(1-F2)-F1) + teta - 1 )*u^(-2)
       P <- ph*(2*teta*u^(-2)*((1-F2)*(2*p*u-z*zt*teta*(1-F2))+zt)-p^2) - p*eta1
       A <- i1*2*teta*zt*u^(-1)*(z*u^(-1)-p*(1-F1))*ph
    }
    if(ver!=1){
       b <- (1-teta^2)*zt*(1-F1)*( 2*z*(1-F2)*u^(-1) - F1*u^(-2) )
       E <- teta*u^(-1)*(1-F1)*(z-3)*zt
       B <- i1*2*(1-F1)*u^(-1)*zt*( (1-teta^2)*z*u^(-1) - teta*b )
       if(ver!=2)  h14 <- i1*( 
	      (1-teta^2)*zt*u^(-2)*(2*(1-F2)*(p*zt^(-1)*u*(1-F1)+z*(2*u-1)+teta*b*zt^(-1)*(1-teta^2)^(-1)*u)-2*F1+1) - b*p
	      )*ph
       h44 <- i1*( -b^2 +
	   (1-teta^2)*2*teta^(-1)*(u^(-1)-1)*(b+(1-teta^2)*zt*(1-F1)*(z*u^(-1)*(1-F2)-F1*u^(-2)))
	   - 2*teta*b 
	   )
    }
  }
  
  if(cop=="G0"){

    lnF1 <- -log(F1)
    lnF2 <- -log(F2)
    u <- lnF1^teta + lnF2^teta
    ut <- u^(1/teta)
    z <- exp(-ut)*F2^(-1)*ut*u^(-1)*lnF2^(teta-1);
    z <- ifelse(z<0.9999999,z,0.9999999)
    zt <- z/(z-1)
    h <- zt*F2^(-1)
    p <- lnF1^(teta-1)*F1^(-1)*u^(-1)*(1-teta-ut)
    d <- lnF2^(-1)*(-F1*lnF1*p-ut)+1
    E <- (teta-1)*F2^(-1)*lnF2^(-2)*( lnF2^teta*u^(-1)*(1-teta-ut+u^(-1)*lnF2^teta*(teta+ut)) - 1 )
    if(ver!=2){ 
       P <- zt*p*(
           F1^(-1)*ph*( 1+lnF1^(-1)*(teta-1-lnF1^teta*u^(-1)*(teta+ut*(1-teta-ut)^(-1))) - p*F1*(z-1)^(-1) )
	   - eta1 )
       A <- i1*h*p*( (teta-1)*(teta+ut)*u^(-1)*lnF2^(teta-1)*(1-teta-ut)^(-1) - d/(z-1) )*ph
    }
    if(ver!=1){ 
       C <- lnF1^teta*log(lnF1)+lnF2^teta*log(lnF2)
       b <- teta^(-1)*(C/u-log(u)/teta)*(1-ut)+log(lnF2)-C/u
       Ct <- (lnF1^teta)*(log(lnF1))^2+(lnF2^teta)*(log(lnF2))^2
       B <- i1*h*(teta-1)*( lnF2^(teta-1)*u^(-1)*
	( log(lnF2)*(1-teta-ut) + (teta-1)*teta^(-1)*C*u^(-1)*(teta+ut) - 1 + ut*log(u)*teta^(-2) )
	 - b*d*(z-1)^(-1) + lnF2^(-1) )
       if(ver==0)  h14 <- i1*(teta-1)*zt*p*(
	        log(lnF1)-(1-teta-ut)^(-1)*((1-teta)*teta^(-1)*C*u^(-1)*(teta+ut)+1-ut*log(u)*teta^(-2)) - b/(z-1)
	       )*ph
       h44 <- i1*zt*(teta-1)*( b + (teta-1)*(
		teta^(-1)*(1-ut)*(Ct/u-C^2/u^2-2*teta^(-1)*(C/u-log(u)/teta)) - 
		ut*teta^(-2)*(C/u-log(u)/teta)^2 - b^2/(z-1) - Ct/u + C^2/u^2   
            ) )
       b <- zt*b*(teta-1)
    }

    E <- E/d+F2^(-1)*(1-d/(z-1))
    p <- zt*p    
    h <- d*h

  }

  if(cop=="G270"){

    F2b <- 1-F2
    lnF1 <- -log(F1)
    lnF2b <- -log(F2b)
    u <- lnF1^teta + lnF2b^teta
    ut <- u^(1/teta)
    z <- exp(-ut)*F2b^(-1)*ut*u^(-1)*lnF2b^(teta-1);
    z <- ifelse(z<0.9999999,z,0.9999999)
    zt <- z/(z-1)
    h <- - zt*F2b^(-1)
    p <- lnF1^(teta-1)*F1^(-1)*u^(-1)*(1-teta-ut)
    d <- lnF2b^(-1)*(-F1*lnF1*p-ut)+1
    E <- (teta-1)*F2b^(-1)*lnF2b^(-2)*( lnF2b^teta*u^(-1)*(1-teta-ut+u^(-1)*lnF2b^teta*(teta+ut)) - 1 )
    if(ver!=2){ 
       P <- zt*p*(
           F1^(-1)*ph*( 1+lnF1^(-1)*(teta-1-lnF1^teta*u^(-1)*(teta+ut*(1-teta-ut)^(-1))) - p*F1*(z-1)^(-1) )
	   - eta1 )
       A <- i1*h*p*( (teta-1)*(teta+ut)*u^(-1)*lnF2b^(teta-1)*(1-teta-ut)^(-1) - d/(z-1) )*ph
    }
    if(ver!=1){ 
       C <- lnF1^teta*log(lnF1)+lnF2b^teta*log(lnF2b)
       b <- teta^(-1)*(C/u-log(u)/teta)*(1-ut)+log(lnF2b)-C/u
       Ct <- (lnF1^teta)*(log(lnF1))^2+(lnF2b^teta)*(log(lnF2b))^2
       B <- i1*h*(teta-1)*( lnF2b^(teta-1)*u^(-1)*
	( log(lnF2b)*(1-teta-ut) + (teta-1)*teta^(-1)*C*u^(-1)*(teta+ut) - 1 + ut*log(u)*teta^(-2) )
	 - b*d*(z-1)^(-1) + lnF2b^(-1) )
   h14 <- i1*(teta-1)*zt*p*(
	        log(lnF1)-(1-teta-ut)^(-1)*((1-teta)*teta^(-1)*C*u^(-1)*(teta+ut)+1-ut*log(u)*teta^(-2)) - b/(z-1)
	       )*ph
       h44 <- i1*zt*(teta-1)*( b + (teta-1)*(
		teta^(-1)*(1-ut)*(Ct/u-C^2/u^2-2*teta^(-1)*(C/u-log(u)/teta)) - 
		ut*teta^(-2)*(C/u-log(u)/teta)^2 - b^2/(z-1) - Ct/u + C^2/u^2   
            ) )
       b <- zt*b*(teta-1)
    }

    E <- - E/d - F2b^(-1)*(1-d/(z-1))
    p <- zt*p    
    h <- d*h

  }

  if(cop=="G90"){
  
    F1b <- 1-F1
    lnF1b <- -log(F1b)
    lnF2 <- -log(F2)
    u <- lnF1b^teta + lnF2^teta
    ut <- u^(1/teta)
    z <- exp(-ut)*F2^(-1)*ut*u^(-1)*lnF2^(teta-1)
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    p <- -lnF1b^(teta-1)*F1b^(-1)*u^(-1)*(1-teta-ut)
    h <- F2^(-1)*( lnF2^(-1)*(F1b*lnF1b*p-ut) + 1 )
    E <- F2^(-1)*( 1+h^(-1)*F2^(-1)*(teta-1)*lnF2^(-2)*( -1+lnF2^teta*( (1-teta-ut)*u^(-1)+(teta+ut)*u^(-2)*lnF2^teta ) ) )
    if(ver!=2){
       P <- -p*( F1b^(-1)*ph*( (teta-1)*lnF1b^(-1) + 1 + (teta-1)*lnF1b^(teta-1)*u^(-1)*(ut+teta)*(1-teta-ut)^(-1)) + eta1 )
       A <- i1*(teta-1)*p*lnF2^(teta-1)*F2^(-1)*u^(-1)*(teta+ut)*(1-teta-ut)^(-1)*ph
    }
    if(ver!=1){
       C <- lnF1b^teta*log(lnF1b) + lnF2^teta*log(lnF2)
       b <- (teta-1)*( teta^(-1)*(C/u-log(u)/teta)*(1-ut) + log(lnF2) - C/u )
       Ct <- lnF1b^teta*(log(lnF1b))^2 + lnF2^teta*(log(lnF2))^2
       B <- i1*(teta-1)*F2^(-1)*lnF2^(-1)*( 1+lnF2^teta*( (1-teta-ut)*u^(-1)*(log(lnF2)-C/u)-u^(-1)-teta^(-1)*u^(1/teta-1)*(C/u-log(u)*teta^(-1)) ) )
       if(ver!=2)  h14 <- i1*(teta-1)*p*(
	      log(lnF1b)-C/u-(1-teta-ut)^(-1)*(1+ut*teta^(-1)*(C/u-log(u)*teta^(-1))) 
	      )*ph
       h44 <- i1*( b +(teta-1)^2*(
                          (1-teta-ut)*teta^(-1)*(Ct/u-C^2/u^2) - teta^(-2)*(C/u-log(u)/teta)*(2*(1-ut)+ut*(C/u-log(u)/teta)) 
		) )
    }

    z <- 1-z

  }

  if(cop=="G180"){
  
    F1b <- 1-F1
    F2b <- 1-F2
    lnF1b <- -log(F1b)
    lnF2b <- -log(F2b)
    u <- lnF1b^teta + lnF2b^teta
    ut <- u^(1/teta)
    z <- exp(-ut)*F2b^(-1)*ut*u^(-1)*lnF2b^(teta-1)
    z <- ifelse(z<0.99999999,z,0.99999999)
    z <- ifelse(z>0.00000001,z,0.00000001)
    p <- -lnF1b^(teta-1)*F1b^(-1)*u^(-1)*(1-teta-ut)
    h <- - F2b^(-1)*( lnF2b^(-1)*(F1b*lnF1b*p-ut) + 1 )
    E <- - F2b^(-1)*( 1 - h^(-1)*F2b^(-1)*(teta-1)*lnF2b^(-2)*( -1+lnF2b^teta*( (1-teta-ut)*u^(-1)+(teta+ut)*u^(-2)*lnF2b^teta ) ) )
    if(ver!=2){
       P <- -p*( F1b^(-1)*ph*( (teta-1)*lnF1b^(-1) + 1 + (teta-1)*lnF1b^(teta-1)*u^(-1)*(ut+teta)*(1-teta-ut)^(-1)) + eta1 )
       A <- - i1*(teta-1)*p*lnF2b^(teta-1)*F2b^(-1)*u^(-1)*(teta+ut)*(1-teta-ut)^(-1)*ph
    }
    if(ver!=1){
       C <- lnF1b^teta*log(lnF1b) + lnF2b^teta*log(lnF2b)
       b <- (teta-1)*( teta^(-1)*(C/u-log(u)/teta)*(1-ut) + log(lnF2b) - C/u )
       Ct <- lnF1b^teta*(log(lnF1b))^2 + lnF2b^teta*(log(lnF2b))^2
       B <- - i1*(teta-1)*F2b^(-1)*lnF2b^(-1)*( 1+lnF2b^teta*( (1-teta-ut)*u^(-1)*(log(lnF2b)-C/u)-u^(-1)-teta^(-1)*u^(1/teta-1)*(C/u-log(u)*teta^(-1)) ) )
       if(ver!=2)  h14 <- i1*(teta-1)*p*(
	      log(lnF1b)-C/u-(1-teta-ut)^(-1)*(1+ut*teta^(-1)*(C/u-log(u)*teta^(-1))) 
	      )*ph
       h44 <- i1*( b +(teta-1)^2*(
                          (1-teta-ut)*teta^(-1)*(Ct/u-C^2/u^2) - teta^(-2)*(C/u-log(u)/teta)*(2*(1-ut)+ut*(C/u-log(u)/teta)) 
		) )
    }

    z <- 1-z

  }

}





 if(ver==1) b <- B <- h14 <- h44 <- 0
 if(ver==2) p <- P <- A <- h14 <- 0


 list(z=z,p=p,h=h,b=b,P=P,A=A,h14=h14,E=E,B=B,h44=h44)



}