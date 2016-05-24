copgHs <- function(p1, p2, eta1 = NULL, eta2 = NULL, teta, teta.st, BivD){

########################################################################################

cjg <- c("C0","J0","G0","C90","J90","G90","C180","J180","G180","C270","J270","G270")

if(BivD %in% cjg) {

derteta.derteta.st <- der2teta.derteta.stteta.st <- exp(teta.st) 
   
}   


if(BivD %in% c("N","FGM","AMH") ) {

derteta.derteta.st <- 1/cosh(teta.st)^2
der2teta.derteta.stteta.st <- -(2 * (sinh(teta.st) * cosh(teta.st))/(cosh(teta.st)^2)^2)

       
} 


if(BivD %in% c("F") ) {

derteta.derteta.st <- 1
der2teta.derteta.stteta.st <- 0
       
} 



   
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



if(BivD=="N"){


c.copula.be1 <- pnorm( (qnorm(p2) - teta*qnorm(p1))/sqrt(1 - teta^2)   )                            
c.copula.be2 <- pnorm( (qnorm(p1) - teta*qnorm(p2))/sqrt(1 - teta^2)   )  

c.copula.thet <- dbinorm(qnorm(p1),qnorm(p2), cov12=teta)

c.copula2.be1 <- dnorm((qnorm(p2)-teta*qnorm(p1))/sqrt(1 - teta^2))  * sqrt(2*pi) *(-teta)/sqrt(1 - teta^2)/exp(-qnorm(p1)^2/2)     
c.copula2.be2 <- dnorm((qnorm(p1)-teta*qnorm(p2))/sqrt(1 - teta^2))  * sqrt(2*pi) *(-teta)/sqrt(1 - teta^2)/exp(-qnorm(p2)^2/2) 

c.copula2.be1be2 <- 1/sqrt(1 - teta^2)*exp(  - (teta^2*( qnorm(p1)^2 +  qnorm(p2)^2 ) - 2*teta*qnorm(p1)*qnorm(p2) ) / (2*(1 - teta^2)) ) 

c.copula2.be1t <- (-(dnorm((qnorm(p2) - tanh(teta.st) * qnorm(p1))/sqrt(1 - tanh(teta.st)^2)) * 
     (1/cosh(teta.st)^2 * qnorm(p1)/sqrt(1 - tanh(teta.st)^2) - (qnorm(p2) - 
         tanh(teta.st) * qnorm(p1)) * (0.5 * (2 * (1/cosh(teta.st)^2 * 
         tanh(teta.st)) * (1 - tanh(teta.st)^2)^-0.5))/sqrt(1 - 
         tanh(teta.st)^2)^2)))/derteta.derteta.st                                                       
                                                                                                                                                                      
c.copula2.be2t <- (-(dnorm((qnorm(p1) - tanh(teta.st) * qnorm(p2))/sqrt(1 - tanh(teta.st)^2)) *  # this as well
    (1/cosh(teta.st)^2 * qnorm(p2)/sqrt(1 - tanh(teta.st)^2) - (qnorm(p1) - 
        tanh(teta.st) * qnorm(p2)) * (0.5 * (2 * (1/cosh(teta.st)^2 * 
        tanh(teta.st)) * (1 - tanh(teta.st)^2)^-0.5))/sqrt(1 - 
        tanh(teta.st)^2)^2)))/derteta.derteta.st
        
bit1.th2ATE <- (0.5 * (pi * teta/(pi * sqrt(1 - teta^2))^2) - 0.5 * 
    ((teta * (qnorm(p1) * (qnorm(p1) - 2 * (teta * qnorm(p2))) + 
        qnorm(p2)^2)/(1 - teta^2) - qnorm(p1) * qnorm(p2))/(pi * 
        (1 - teta^2)))) * exp(-(0.5 * ((qnorm(p1) * (qnorm(p1) - 
    2 * (teta * qnorm(p2))) + qnorm(p2)^2)/(1 - teta^2))))/sqrt(1 - 
    teta^2)


}





if(BivD=="F"){

# 1 - exp(-teta) = -expm1(-teta)
# recall log1p as well if needed. 
# epcl <- -expm1(-teta) 

  c.copula.be1 <- (exp(teta)* (-1 + exp(p2* teta)))/(-exp((p1 + p2)* teta) + exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))

  c.copula.be2 <- (exp(teta)* (-1 + exp(p1* teta)))/(-exp((p1 + p2)* teta) + exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))

  c.copula.thet <-  (exp(teta)* (1/(-1 + exp(teta)) + (-1 - exp(p2* teta)* (-1 + p1) + p1 - 
     exp(p1* teta)* (-1 + p2) + p2)/(
    exp((p1 + p2)* teta) - exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta))))* teta + 
 log((exp(-(p1 + p2)* teta)* (-exp((p1 + p2)* teta) + 
     exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta))))/(-1 + exp(teta))))/teta^2 


c.copula2.be1 <-   (exp(teta + 
  p1* teta)* (-1 + exp(p2* teta))* (-exp(teta) + exp(
   p2* teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2
  
 
 c.copula2.be2 <-     (exp(teta + 
  p2* teta)* (-1 + exp(p1* teta))* (-exp(teta) + exp(
   p1* teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))^2


c.copula2.be1be2 <- (exp((1 + p1 + p2)* teta)* (-1 + exp(teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2


c.copula2.be1t <- (exp(teta + 
  p1 *teta)* (exp(2* p2* teta)* (-1 + p1) + exp(teta)* p1 - 
   exp(p2* teta)* (-1 + p1 + exp(teta)* (p1 - p2) + p2)))/(exp((p1 + 
     p2)* teta) - exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2

c.copula2.be2t <-  (exp(teta + 
  p2 *teta)* (exp(2* p1* teta)* (-1 + p2) + exp(teta)* p2 - 
   exp(p1* teta)* (-1 + p2 + exp(teta)* (p2 - p1) + p1)))/(exp((p2 + 
     p1)* teta) - exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))^2


bit1.th2 <- bit1.th2ATE <- 1/teta^2 * ((1/(1 - exp(-teta)) * (exp(-teta) - (exp(-teta * 
    p1) * p1 * (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * 
    (exp(-teta * p2) * p2))) - exp(-teta)/(1 - exp(-teta))^2 * 
    ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2))))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))) - 2 * teta/(teta^2)^2 * 
    log(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - exp(-teta * 
        p1)) * (1 - exp(-teta * p2)))) + (1/teta^2 * ((1/(1 - 
    exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
    exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    p2))) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - exp(-teta)) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(1 - exp(-teta)) * 
    ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2))))) - -1/teta * ((1/(1 - exp(-teta)) * (exp(-teta) + 
    (exp(-teta * p1) * p1 * (exp(-teta * p2) * p2) - exp(-teta * 
        p1) * p1 * p1 * (1 - exp(-teta * p2)) + (exp(-teta * 
        p1) * p1 * (exp(-teta * p2) * p2) - (1 - exp(-teta * 
        p1)) * (exp(-teta * p2) * p2 * p2)))) + exp(-teta)/(1 - 
    exp(-teta))^2 * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
    exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    p2))) + (exp(-teta)/(1 - exp(-teta))^2 * (exp(-teta) - (exp(-teta * 
    p1) * p1 * (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * 
    (exp(-teta * p2) * p2))) - (exp(-teta)/(1 - exp(-teta))^2 + 
    exp(-teta) * (2 * (exp(-teta) * (1 - exp(-teta))))/((1 - 
        exp(-teta))^2)^2) * ((1 - exp(-teta)) - (1 - exp(-teta * 
    p1)) * (1 - exp(-teta * p2)))))/(1/(1 - exp(-teta)) * ((1 - 
    exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) + 
    (1/(1 - exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * 
        (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * 
        p2) * p2))) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - exp(-teta)) - 
        (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) * (1/(1 - 
        exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
        exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * 
        p2) * p2))) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - exp(-teta)) - 
        (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(1 - 
        exp(-teta)) * ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * 
        (1 - exp(-teta * p2))))^2))
    
    
    
    
    
    
    


}





if(BivD %in% c("C0","C90","C180","C270")){



c.copula.be1 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p1^((-teta) - 1) * (-teta)))

c.copula.be2 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p2^((-teta) - 1) * (-teta)))
  
  
  c.copula.thet <- ((-1 + p1^-teta + p2^-teta)^(-1/
  teta) *((teta *(p2^teta *log(p1) + p1^teta* log(p2)))/(
   p2^teta - p1^teta* (-1 + p2^teta)) + 
   log(-1 + p1^-teta + p2^-teta)))/teta^2
  


c.copula2.be1 <- (p1^(-2 + teta)* p2^teta* (-1 + p1^-teta + p2^-teta)^(-1/
  teta)* (-1 + p2^teta)* (1 + teta))/(p2^teta - 
  p1^teta* (-1 + p2^teta))^2
 

 c.copula2.be2 <- (p2^(-2 + teta)* p1^teta* (-1 + p2^-teta + p1^-teta)^(-1/
  teta)* (-1 + p1^teta)* (1 + teta))/(p1^teta - 
  p2^teta* (-1 + p1^teta))^2


c.copula2.be1be2 <- p1^(-1 - teta)* p2^(-1 - teta)* (-1 + p1^-teta + p2^-teta)^(-2 - 1/
  teta) *(1 + teta)
   
  
c.copula2.be1t <- ((1 + 1/teta) * (log(p1)/p1^teta + log(p2)/p2^teta)/(1/p1^teta + 
    1/p2^teta - 1)^(1/teta + 2) + log(1/p1^teta + 1/p2^teta - 
    1)/(teta^2 * (1/p1^teta + 1/p2^teta - 1)^(1 + 1/teta)))/p1^(1 + 
    teta) + (1/p1^(1 + teta) - (1/p1^(1 + teta) + teta * log(p1)/p1^(1 + 
    teta)))/(teta * (1/p1^teta + 1/p2^teta - 1)^(1 + 1/teta)) 
  
 c.copula2.be2t <- ((1 + 1/teta) * (log(p1)/p1^teta + log(p2)/p2^teta)/(1/p1^teta + 
    1/p2^teta - 1)^(1/teta + 2) + log(1/p1^teta + 1/p2^teta - 
    1)/(teta^2 * (1/p1^teta + 1/p2^teta - 1)^(1 + 1/teta)))/p2^(1 + 
    teta) + (1/p2^(1 + teta) - (1/p2^(1 + teta) + teta * log(p2)/p2^(1 + 
    teta)))/(teta * (1/p1^teta + 1/p2^teta - 1)^(1 + 1/teta))
  
bit1.th2ATE <- (teta * (2 * (p1^teta * log(p2) * (p2^teta * (p1^teta - 
    1) - p1^teta) * (teta - log(1/p1^teta + 1/p2^teta - 1))) + 
    2 * (p2^teta * ((p2^teta * (p1^teta - 1) - p1^teta) * (teta - 
        log(1/p1^teta + 1/p2^teta - 1)) + p1^teta * teta * (1 + 
        teta) * log(p2)) * log(p1)) + teta * (p1^teta * log(p2)^2 * 
    (p1^teta + p2^teta * teta * (p1^teta - 1)) + p2^teta * log(p1)^2 * 
    (p1^teta * teta * (p2^teta - 1) + p2^teta))) - (2 * teta - 
    log(1/p1^teta + 1/p2^teta - 1)) * log(1/p1^teta + 1/p2^teta - 
    1) * (p1^teta + p2^teta - p1^teta * p2^teta)^2)/(teta^4 * 
    (1/p1^teta + 1/p2^teta - 1)^(1/teta) * (p2^teta - p1^teta * 
    (p2^teta - 1))^2)
   
    
}








if(BivD %in% c("G0","G90","G180","G270")){


  c.copula.be1 <- (exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/teta)))* (-log(p1))^(-1 + 
  teta)* ((-log(p1))^teta + (-log(p2))^teta)^(-1 + 1/teta))/p1


  c.copula.be2 <- (exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/teta)))* (-log(p2))^(-1 + 
  teta)* ((-log(p2))^teta + (-log(p1))^teta)^(-1 + 1/teta))/p2


  c.copula.thet <-   (1/(teta^2))*exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/
  teta)))* ((-log(p1))^teta + (-log(p2))^teta)^(-1 + 1/
  teta)* (-teta* (-log(p1))^
    teta* log(-log(p1)) + ((-log(p1))^teta + (-log(p2))^
      teta)* log((-log(p1))^teta + (-log(p2))^teta) - 
   teta *(-log(p2))^teta* log(-log(p2)))

  
c.copula2.be1 <-(1/(p1^2))*exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/
  teta)))* (-log(p1))^(-2 + 
  teta) *((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* ((-log(p1))^
    teta* (log(p1) + ((-log(p1))^teta + (-log(p2))^teta)^(1/
      teta)) + (1 - teta + log(p1))* (-log(p2))^teta)

                
c.copula2.be2 <- (1/(p2^2))*exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/
  teta)))* (-log(p2))^(-2 + 
  teta) *((-log(p2))^teta + (-log(p1))^teta)^(-2 + 1/
  teta)* ((-log(p2))^
    teta* (log(p2) + ((-log(p2))^teta + (-log(p1))^teta)^(1/
      teta)) + (1 - teta + log(p2))* (-log(p1))^teta)


c.copula2.be1be2 <- (exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/teta)))* (-log(p1))^(-1 + 
  teta) *(-1 + teta + ((-log(p1))^teta + (-log(p2))^teta)^(1/
   teta))* ((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* (-log(p2))^(-1 + teta))/(p1 *p2)


 c.copula2.be1t <- ((-log(p1))^(-1 + teta) * (((-log(p1))^teta * log(-log(p1)) + 
    (-log(p2))^teta * log(-log(p2))) * ((-log(p1))^teta + (-log(p2))^teta)^(1/teta - 
    2) * (1/teta - 1) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta - 
    1) * log((-log(p1))^teta + (-log(p2))^teta)/teta^2) + ((-log(p1))^(-1 + 
    teta) * log(-log(p1)) - (-log(p1))^(-1 + teta) * (((-log(p1))^teta * 
    log(-log(p1)) + (-log(p2))^teta * log(-log(p2))) * ((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta - 1) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * 
    log((-log(p1))^teta + (-log(p2))^teta)/teta)/teta) * ((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta - 1)) * exp(-((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta))/p1
  
 c.copula2.be2t <- ((-log(p2))^(-1 + teta) * (((-log(p1))^teta * log(-log(p1)) + 
    (-log(p2))^teta * log(-log(p2))) * ((-log(p1))^teta + (-log(p2))^teta)^(1/teta - 
    2) * (1/teta - 1) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta - 
    1) * log((-log(p1))^teta + (-log(p2))^teta)/teta^2) + ((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta - 1) * ((-log(p2))^(-1 + teta) * 
    log(-log(p2)) - (-log(p2))^(-1 + teta) * (((-log(p1))^teta * 
    log(-log(p1)) + (-log(p2))^teta * log(-log(p2))) * ((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta - 1) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * 
    log((-log(p1))^teta + (-log(p2))^teta)/teta)/teta)) * exp(-((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta))/p2


bit1.th2ATE <- ((-log(p1))^teta + (-log(p2))^teta)^(1/teta - 2) * 
    (((-log(p1))^teta + (-log(p2))^teta) * (((-log(p1))^teta + 
        (-log(p2))^teta) * (((-log(p1))^teta + (-log(p2))^teta)^(1/teta) - 
        1) * log((-log(p1))^teta + (-log(p2))^teta) - 2 * (teta * 
        ((-log(p1))^teta + (-log(p2))^teta * ((((-log(p1))^teta + 
            (-log(p2))^teta)^(1/teta) - 1) * log(-log(p2)) + 
            1)))) * log((-log(p1))^teta + (-log(p2))^teta) + 
        teta * (2 * ((-log(p1))^teta * log(-log(p1)) * (teta * 
            ((-log(p1))^teta + (-log(p2))^teta * ((((-log(p1))^teta + 
                (-log(p2))^teta)^(1/teta) + teta - 1) * log(-log(p2)) + 
                1)) - ((-log(p1))^teta + (-log(p2))^teta) * (((-log(p1))^teta + 
            (-log(p2))^teta)^(1/teta) - 1) * log((-log(p1))^teta + 
            (-log(p2))^teta))) + teta * ((-log(p1))^teta * ((-log(p1))^teta * 
            (((-log(p1))^teta + (-log(p2))^teta)^(1/teta) - 1) - 
            teta * (-log(p2))^teta) * log(-log(p1))^2 + (-log(p2))^teta * 
            ((-log(p1))^teta * (2 - teta * log(-log(p2))) + (-log(p2))^teta * 
                ((((-log(p1))^teta + (-log(p2))^teta)^(1/teta) - 
                  1) * log(-log(p2)) + 2)) * log(-log(p2))))) * 
    exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta))/teta^4


}






if(BivD %in% c("J0","J90","J180","J270")){

  c.copula.be1 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p1)^(teta - 1) * teta - (1 - p1)^(teta - 
    1) * teta * (1 - p2)^teta))


  c.copula.be2 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p2)^(teta - 1) * teta - (1 - p1)^teta * 
    ((1 - p2)^(teta - 1) * teta)))


  c.copula.thet <-    (((1 - p1)^
   teta - (-1 + (1 - p1)^teta) *(1 - p2)^
    teta)^(1/teta)* (log((1 - p1)^
     teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta) + (
   teta* ((1 - p1)^
       teta* (-1 + (1 - p2)^teta)* log(
        1 - p1) + (-1 + (1 - p1)^teta) *(1 - p2)^
       teta* log(1 - p2)))/((1 - p1)^
    teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)))/teta^2
        
  
c.copula2.be1 <- (1 - p1)^(-2 + 
  teta)* (-1 + (1 - p2)^teta)* ((1 - p1)^
   teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
  teta) *(1 - p2)^teta* (-1 + teta)

       
c.copula2.be2 <- (1 - p2)^(-2 + 
  teta)* (-1 + (1 - p1)^teta)* ((1 - p2)^
   teta - (-1 + (1 - p2)^teta)* (1 - p1)^teta)^(-2 + 1/
  teta) *(1 - p1)^teta* (-1 + teta)


c.copula2.be1be2 <- (1 - p1)^(-1 + 
  teta)* ((1 - p1)^teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
  teta) *(1 - p2)^(-1 + 
  teta)* (-(-1 + (1 - p1)^teta)* (-1 + (1 - p2)^teta) + teta)

      
 c.copula2.be1t <- ((((1 - p1)^teta - (1 - p1)^teta * (1 - p2)^teta) * log(1 - p1) + 
    ((1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta) * log(1 - 
        p2)) * ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)^(1/teta - 2) * (1/teta - 1) - ((1 - p1)^teta + 
    (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(1/teta - 
    1) * log((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)/teta^2) * ((1 - p1)^(teta - 1) - (1 - p1)^(teta - 
    1) * (1 - p2)^teta) + ((1 - p1)^(teta - 1) + (1 - p1)^(teta - 
    1) * (1 - p2)^teta + teta * ((1 - p1)^(teta - 1) * log(1 - 
    p1) - (1 - p1)^(teta - 1) * (1 - p2)^teta * log(1 - p2)) - 
    (((1 - p1)^(teta - 1) + teta * (1 - p1)^(teta - 1) * log(1 - 
        p1)) * (1 - p2)^teta + (1 - p1)^(teta - 1))) * ((1 - 
    p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(1/teta - 
    1)/teta
  
  
 c.copula2.be2t <- ((((1 - p1)^teta - (1 - p1)^teta * (1 - p2)^teta) * log(1 - p1) + 
    ((1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta) * log(1 - 
        p2)) * ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)^(1/teta - 2) * (1/teta - 1) - ((1 - p1)^teta + 
    (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(1/teta - 
    1) * log((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)/teta^2) * ((1 - p2)^(teta - 1) - (1 - p1)^teta * 
    (1 - p2)^(teta - 1)) + ((1 - p1)^teta * (1 - p2)^(teta - 
    1) + (1 - p2)^(teta - 1) + teta * ((1 - p2)^(teta - 1) * 
    log(1 - p2) - (1 - p1)^teta * (1 - p2)^(teta - 1) * log(1 - 
    p1)) - (((1 - p2)^(teta - 1) + teta * (1 - p2)^(teta - 1) * 
    log(1 - p2)) * (1 - p1)^teta + (1 - p2)^(teta - 1))) * ((1 - 
    p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(1/teta - 
    1)/teta

  
bit1.th2ATE <- ((1 - p1)^teta - ((1 - p1)^teta - 1) * (1 - p2)^teta)^(1/teta - 
    2) * ((2 * (teta * (((1 - p1)^teta - 1) * (1 - p2)^teta - 
    (1 - p1)^teta) * (((1 - p1)^teta - 1) * (1 - p2)^teta * (log(1 - 
    p2) - 1) + (1 - p1)^teta)) - ((1 - p1)^teta - ((1 - p1)^teta - 
    1) * (1 - p2)^teta)^2 * log((1 - p1)^teta - ((1 - p1)^teta - 
    1) * (1 - p2)^teta)) * log((1 - p1)^teta - ((1 - p1)^teta - 
    1) * (1 - p2)^teta) + teta * (2 * (((((1 - p1)^teta - 1) * 
    (1 - p2)^teta - (1 - p1)^teta) * ((1 - p2)^teta - 1) * log((1 - 
    p1)^teta - ((1 - p1)^teta - 1) * (1 - p2)^teta) + teta * 
    ((((1 - p1)^teta + teta - 1) * log(1 - p2) + 1 - 2 * (1 - 
        p1)^teta) * (1 - p2)^teta + (1 - p1)^teta - ((1 - p1)^teta - 
        1) * (1 - p2)^(2 * teta) * (log(1 - p2) - 1))) * (1 - 
    p1)^teta * log(1 - p1)) + teta * (((1 - p1)^teta - 1) * ((1 - 
    p1)^teta * (teta * log(1 - p2) - 2) - ((1 - p1)^teta - 1) * 
    (1 - p2)^teta * (log(1 - p2) - 2)) * (1 - p2)^teta * log(1 - 
    p2) - (((1 - p2)^teta - 1) * (1 - p1)^teta - teta * (1 - 
    p2)^teta) * ((1 - p2)^teta - 1) * (1 - p1)^teta * log(1 - 
    p1)^2)))/teta^4
 
}








if(BivD == "AMH"){

  c.copula.be1 <- p2 * (1 - p1 * teta * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - 
                  teta * (1 - p1) * (1 - p2))


  c.copula.be2 <- p1 * (1 - p2 * teta * (1 - p1)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - 
                  teta * (1 - p1) * (1 - p2))


  c.copula.thet <- p1 * p2 * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))^2
        
  
c.copula2.be1 <- -(p2 * teta * (1 - p2) * (2 - 2 * (p1 * teta * (1 - p2)/(1 - 
    teta * (1 - p1) * (1 - p2))))/(1 - teta * (1 - p1) * (1 - 
    p2))^2)

       
c.copula2.be2 <- -(p1 * teta * (1 - p1) * (2 - 2 * (p2 * teta * (1 - p1)/(1 - 
    teta * (1 - p1) * (1 - p2))))/(1 - teta * (1 - p1) * (1 - 
    p2))^2)



c.copula2.be1be2 <- (1 - teta * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - 
    p2)/(1 - teta * (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - 
    teta * (1 - p1) * (1 - p2)))/(1 - teta * (1 - p1) * (1 - 
    p2))

      
 c.copula2.be1t <- p2 * (1 - p1 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * 
    (1 - p1) * (1 - p2))))) * (1 - p2)/(1 - teta * (1 - p1) * 
    (1 - p2))^2

  
  
 c.copula2.be2t <- p1 * (1 - p1) * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - 
    teta * (1 - p1) * (1 - p2)))))/(1 - teta * (1 - p1) * (1 - 
    p2))^2

  
bit1.th2ATE <- 2 * (p1 * p2 * (1 - p1)^2 * (1 - p2)^2/(1 - teta * (1 - p1) * 
    (1 - p2))^3)
 
}







if(BivD == "FGM"){

  c.copula.be1 <- p2 * (1 + teta * (1 - 2 * p1) * (1 - p2))


  c.copula.be2 <- p1 * (1 + teta * (1 - 2 * p2) * (1 - p1))


  c.copula.thet <- p1 * p2 * (1 - p1) * (1 - p2)

        
  
c.copula2.be1 <- -(2 * (p2 * teta * (1 - p2)))


       
c.copula2.be2 <- -(2 * (p1 * teta * (1 - p1)))


c.copula2.be1be2 <- 1 + teta * (1 - 2 * p1) * (1 - 2 * p2)


      
 c.copula2.be1t <- p2 * (1 - 2 * p1) * (1 - p2)
  
  
 c.copula2.be2t <- p1 * (1 - 2 * p2) * (1 - p1)


  
bit1.th2ATE <- 0
 
}










#########################
# modular derivatives

c.copula.theta  <- c.copula.thet*derteta.derteta.st
c.copula2.be1th <- c.copula2.be1t*derteta.derteta.st
c.copula2.be2th <- c.copula2.be2t*derteta.derteta.st 
bit1.th2 <- bit1.th2ATE*derteta.derteta.st^2 + c.copula.thet*der2teta.derteta.stteta.st   

#########################





if(BivD %in% c("C90","J90","G90") ) {

#c.copula.be1     <- c.copula.be1
c.copula.be2     <- 1 - c.copula.be2
c.copula.theta   <- - c.copula.theta
c.copula2.be1    <- - c.copula2.be1
c.copula2.be2    <- - c.copula2.be2
#c.copula2.be1be2 <- c.copula2.be1be2
#c.copula2.be1th  <- c.copula2.be1th 
c.copula2.be2th  <- - c.copula2.be2th
bit1.th2ATE      <- - bit1.th2ATE  
bit1.th2         <- - bit1.th2 

}  

if(BivD %in% c("C180","J180","G180") ) {

c.copula.be1     <- 1 - c.copula.be1 
c.copula.be2     <- 1 - c.copula.be2
#c.copula.theta   <- c.copula.theta
#c.copula2.be1    <- c.copula2.be1 
#c.copula2.be2    <- c.copula2.be2
#c.copula2.be1be2 <- c.copula2.be1be2
c.copula2.be1th  <- - c.copula2.be1th
c.copula2.be2th  <- - c.copula2.be2th
#bit1.th2ATE      <- bit1.th2ATE
#bit1.th2         <- bit1.th2  


}  


if(BivD %in% c("C270","J270","G270") ) {

c.copula.be1     <- 1 - c.copula.be1
#c.copula.be2     <- c.copula.be2
c.copula.theta   <- - c.copula.theta
c.copula2.be1    <- - c.copula2.be1
c.copula2.be2    <- - c.copula2.be2
#c.copula2.be1be2 <- c.copula2.be1be2
c.copula2.be1th  <- - c.copula2.be1th
#c.copula2.be2th  <-   c.copula2.be2th
bit1.th2ATE      <- - bit1.th2ATE
bit1.th2         <- - bit1.th2

}   




epsilon <- 0.0000001 
max.p   <- 0.9999999

# the bits below are probs
  
c.copula.be2 <- ifelse(c.copula.be2 > max.p, max.p, c.copula.be2) 
c.copula.be2 <- ifelse(c.copula.be2 < epsilon,     epsilon, c.copula.be2)
c.copula.be1 <- ifelse(c.copula.be1 > max.p, max.p, c.copula.be1) 
c.copula.be1 <- ifelse(c.copula.be1 < epsilon,     epsilon, c.copula.be1)
c.copula2.be1be2 <- ifelse(c.copula2.be1be2 < epsilon, epsilon, c.copula2.be1be2)


# this below is a safety check really

ifef <- function(dv){

epsilon <- 0.0000001 
dv <- ifelse(is.na(dv), epsilon, dv ) 
dv <- ifelse(dv == Inf ,  8.218407e+20, dv )
dv <- ifelse(dv == -Inf ,  -8.218407e+20, dv )
dv


}


c.copula.be1     <- ifef(c.copula.be1    )  
c.copula.be2     <- ifef(c.copula.be2    ) 
c.copula.theta   <- ifef(c.copula.theta  ) 
c.copula2.be1    <- ifef(c.copula2.be1   ) 
c.copula2.be2    <- ifef(c.copula2.be2   ) 
c.copula2.be1be2 <- ifef(c.copula2.be1be2) 
c.copula2.be1th  <- ifef(c.copula2.be1th ) 
c.copula2.be2th  <- ifef(c.copula2.be2th ) 
bit1.th2ATE      <- ifef(bit1.th2ATE     ) 
bit1.th2         <- ifef(bit1.th2        ) 



list(
c.copula.be1     = c.copula.be1,    
c.copula.be2     = c.copula.be2,    
c.copula.theta   = c.copula.theta,  
c.copula2.be1    = c.copula2.be1,  
c.copula2.be2    = c.copula2.be2,   
c.copula2.be1be2 = c.copula2.be1be2,
c.copula2.be1th  = c.copula2.be1th, 
c.copula2.be2th  = c.copula2.be2th, 
bit1.th2ATE      = bit1.th2ATE,     
bit1.th2         = bit1.th2 )     


}

