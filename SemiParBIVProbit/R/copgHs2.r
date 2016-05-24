copgHs2 <- function(p1, p2, eta1 = NULL, eta2 = NULL, teta, teta.st, BivD){

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





if(BivD=="FGM"){

c.copula.be2 <- p1 * (1 + teta * (1 - 2 * p2) * (1 - p1))
}


if(BivD=="AMH"){

c.copula.be2 <- p1 * (1 - p2 * teta * (1 - p1)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - 
                  teta * (1 - p1) * (1 - p2))
 

}





if(BivD=="N"){

c.copula.be2 <- pnorm( (qnorm(p1) - teta*qnorm(p2))/sqrt(1 - teta^2)   ) 

}


if(BivD=="F"){
 
  c.copula.be2 <-   (exp(teta)* (-1 + exp(p1* teta)))/(-exp((p1 + p2)* teta) + 
                     exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))

}



if(BivD %in% c("C0","C90","C180","C270")){

c.copula.be2 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p2^((-teta) - 1) * (-teta)))
  
}



if(BivD %in% c("G0","G90","G180","G270")){

  c.copula.be2 <- (exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/teta)))* (-log(p2))^(-1 + 
  teta)* ((-log(p2))^teta + (-log(p1))^teta)^(-1 + 1/teta))/p2

}


if(BivD %in% c("J0","J90","J180","J270")){

  c.copula.be2 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p2)^(teta - 1) * teta - (1 - p1)^teta * 
    ((1 - p2)^(teta - 1) * teta)))

}



if(BivD %in% c("C90","J90","G90") ) {


c.copula.be2  <- 1 - c.copula.be2

}  

if(BivD %in% c("C180","J180","G180") ) {


c.copula.be2  <- 1 - c.copula.be2


}  


#if(BivD %in% c("C270","J270","G270") ) {
#
#
#c.copula.be2 <- c.copula.be2
#
#
#}   



ifef <- function(dv){

epsilon <- 0.0000001 
dv <- ifelse(is.na(dv), epsilon, dv ) 
dv <- ifelse(dv == Inf ,  8.218407e+20, dv )
dv <- ifelse(dv == -Inf ,  -8.218407e+20, dv )
dv

}

# safety check

c.copula.be2 <- ifef(c.copula.be2)




epsilon <- 0.0000001 
max.p   <- 0.9999999
  
c.copula.be2 <- ifelse(c.copula.be2 > max.p, max.p, c.copula.be2) 
c.copula.be2 <- ifelse(c.copula.be2 < epsilon,     epsilon, c.copula.be2)


list( c.copula.be2 = c.copula.be2 )     


}
