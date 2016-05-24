copgHsAT <- function(p1, p2, teta, BivD, Ln = FALSE){


c.copula.be2 <- c.copula2.be1be2 <- 1


########################################################################################
# Rotations
########################################################################################

if(BivD %in% c("C90","J90","G90") ) {
p1 <- 1 - p1 
teta <- -teta
}  

if(BivD %in% c("C180","J180","G180") ) {
p1 <- 1 - p1
p2 <- 1 - p2
}  

if(BivD %in% c("C270","J270","G270") ) {
p2 <- 1 - p2 
teta <- -teta 
}   
   
########################################################################################   
########################################################################################




if(Ln == FALSE){




if(BivD == "AMH"){
                    
c.copula.be2 <- p1 * (1 - p2 * teta * (1 - p1)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - 
                  teta * (1 - p1) * (1 - p2))

}

if(BivD == "FGM"){
                    
c.copula.be2 <- p1 * (1 + teta * (1 - 2 * p2) * (1 - p1))
}




if(BivD == "N"){
                    
c.copula.be2 <- pnorm( (qnorm(p1) - teta*qnorm(p2))/sqrt(1 - teta^2)   )  

}


if(BivD == "F"){

c.copula.be2 <- (exp(teta)* (-1 + exp(p1* teta)))/(-exp((p1 + p2)* teta) + 
                     exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))
    
}


if(BivD %in% c("C0","C90","C180","C270") ){

c.copula.be2 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p2^((-teta) - 1) * (-teta)))

}



if(BivD %in% c("G0","G90","G180","G270") ){

c.copula.be2 <- (exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/teta)))* (-log(p2))^(-1 + 
  teta)* ((-log(p2))^teta + (-log(p1))^teta)^(-1 + 1/teta))/p2

}


if(BivD %in% c("J0","J90","J180","J270") ){

  c.copula.be2 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p2)^(teta - 1) * teta - (1 - p1)^teta * 
    ((1 - p2)^(teta - 1) * teta)))

}



if(BivD %in% c("C90","J90","G90") )    c.copula.be2  <- 1 - c.copula.be2
if(BivD %in% c("C180","J180","G180") ) c.copula.be2  <- 1 - c.copula.be2

epsilon <- 0.0000001 
max.p   <- 0.9999999
c.copula.be2 <- ifelse(c.copula.be2 > max.p, max.p, c.copula.be2) 
c.copula.be2 <- ifelse(c.copula.be2 < epsilon,     epsilon, c.copula.be2)



}







if(Ln == TRUE){






if(BivD == "AMH"){
                    
c.copula2.be1be2 <- (1 - teta * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - 
    p2)/(1 - teta * (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - 
    teta * (1 - p1) * (1 - p2)))/(1 - teta * (1 - p1) * (1 - 
    p2))


}

if(BivD == "FGM"){
                    
c.copula2.be1be2 <- 1 + teta * (1 - 2 * p1) * (1 - 2 * p2)

}






if(BivD == "N"){
                    
c.copula2.be1be2 <- 1/sqrt(1 - teta^2)*exp(  - (teta^2*( qnorm(p1)^2 +  qnorm(p2)^2 ) - 2*teta*qnorm(p1)*qnorm(p2) ) / (2*(1 - teta^2)) ) 

}


if(BivD == "F"){

c.copula2.be1be2 <- (exp((1 + p1 + p2)* teta)* (-1 + exp(teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2
    
}


if(BivD %in% c("C0","C90","C180","C270") ){

c.copula2.be1be2 <- p1^(-1 - teta)* p2^(-1 - teta)* (-1 + p1^-teta + p2^-teta)^(-2 - 1/
  teta) *(1 + teta)
  
}



if(BivD %in% c("G0","G90","G180","G270") ){

c.copula2.be1be2 <- (exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/teta)))* (-log(p1))^(-1 + 
  teta) *(-1 + teta + ((-log(p1))^teta + (-log(p2))^teta)^(1/
   teta))* ((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* (-log(p2))^(-1 + teta))/(p1 *p2)

}


if(BivD %in% c("J0","J90","J180","J270") ){

c.copula2.be1be2 <- (1 - p1)^(-1 + 
  teta)* ((1 - p1)^teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
  teta) *(1 - p2)^(-1 + 
  teta)* (-(-1 + (1 - p1)^teta)* (-1 + (1 - p2)^teta) + teta)

}






#if(BivD %in% c("C270","J270","G270") ) {
#c.copula2.be1be2 <- c.copula2.be1be2
# rotations
#} 


epsilon <- 0.0000001 
max.p   <- 0.9999999
c.copula2.be1be2 <- ifelse(c.copula2.be1be2 < epsilon, epsilon, c.copula2.be1be2)



  
}






ifef <- function(dv){

epsilon <- 0.0000001 
dv <- ifelse(is.na(dv), epsilon, dv ) 
dv <- ifelse(dv == Inf ,  8.218407e+20, dv )
dv <- ifelse(dv == -Inf ,  -8.218407e+20, dv )
dv 

}

# safety check

c.copula.be2     <- ifef(c.copula.be2)
c.copula2.be1be2 <- ifef(c.copula2.be1be2)







list(c.copula.be2 = c.copula.be2, c.copula2.be1be2 = c.copula2.be1be2)     


}




     























