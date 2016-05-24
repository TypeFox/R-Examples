distrHsAT <- function(y2, eta2, sigma2, nu, margin2){


if(margin2 == "probit"){

  pdf2          <- dnorm(y2, 0, 1)
    p2          <- pnorm(y2, 0, 1)
     
}


if( margin2 == "logit" ){
 
  p2  <- plogis(y2)
  pdf2 <- dlogis(y2)
  
  
}



if( margin2 == "cloglog" ){
 
  p2  <- 1-exp(-exp(y2))
  pdf2 <- exp(-exp(y2)) * exp(y2)
  
}


if( margin2 == "cauchit" ){
 
  p2  <- 1 / pi * atan(y2) + 0.5
  pdf2 <- 1 / (pi * (1 + y2^2))
  
  
}




if(margin2 == "N"){

  pdf2          <- dnorm(y2, mean=eta2, sd = sqrt(sigma2))
    p2          <- pnorm(y2, mean=eta2, sd = sqrt(sigma2))
     
}


if(margin2 == "LN"){

  pdf2          <- dlnorm(y2, meanlog=eta2, sdlog = sqrt(sigma2))
    p2          <- plnorm(y2, meanlog=eta2, sdlog = sqrt(sigma2))

}


if(margin2 == "WEI"){

  pdf2          <- sqrt(sigma2)/exp(eta2)*(y2/exp(eta2))^(sqrt(sigma2)-1) * exp(-(y2/exp(eta2))^sqrt(sigma2))
                   
    p2          <-  1-exp(-(y2/exp(eta2))^sqrt(sigma2)) 

}


if(margin2 == "FISK"){

  pdf2          <- sqrt(sigma2)*y2^(sqrt(sigma2)-1) / (exp(eta2)^sqrt(sigma2)*(1+(y2/exp(eta2))^sqrt(sigma2))^2)
                   
    p2          <-  1/(1+(y2/exp(eta2))^-sqrt(sigma2))

}




if(margin2 == "iG"){


#sigma2 <- ifelse(sigma2 < 0.001, 0.001, sigma2)


  pdf2          <- exp(-0.5 * log(2 * pi) - log(sqrt(sigma2)) - (3/2) * log(y2) - 
                   ((y2 - exp(eta2))^2)/(2 * sigma2 * (exp(eta2)^2) * y2))
                   
    p2          <-  pnorm(((y2/exp(eta2)) - 1)/(sqrt(sigma2) * sqrt(y2))) + 
                    exp(2/(exp(eta2)*sigma2))* pnorm(-((y2/exp(eta2)) + 1)/(sqrt(sigma2) * sqrt(y2)))
                
}



if(margin2 == "LO"){

  pdf2          <- dlogis(y2,eta2,sqrt(sigma2)) # exp(-(y2-eta2)/sqrt(sigma2))/(sqrt(sigma2)*(1+exp(-(y2-eta2)/sqrt(sigma2)))^2)
    p2          <- plogis(y2,eta2,sqrt(sigma2)) #1/(1+exp(-(y2-eta2)/sqrt(sigma2)))
                
}


if(margin2 == "rGU"){

  pdf2          <- 1/sqrt(sigma2)*exp(-((y2-eta2)/sqrt(sigma2)+exp(-((y2-eta2)/sqrt(sigma2)))))
    p2          <- exp(-(exp(-(y2-eta2)/sqrt(sigma2))))
                
}



if(margin2 == "GU"){


  
  bit <- (exp((y2 - eta2)/sqrt(sigma2)) * (1/sqrt(sigma2)))

  bit <- ifelse(bit=="Inf", 1e+290, bit)
  
  pdf2          <- exp(-exp((y2 - eta2)/sqrt(sigma2))) * bit
  
    p2          <- 1 - exp(-exp((y2 - eta2)/sqrt(sigma2)))
    

}


if(margin2 %in% c("GA")){

sigma2    <- ifelse(sigma2 < 0.006, 0.006, sigma2) # related to gamma function

 pdf2          <-  dgamma(y2, shape = 1/sigma2, scale = exp(eta2) * sigma2)
                   
    p2          <-  pgamma(y2, shape = 1/sigma2, scale = exp(eta2) * sigma2)
    
    }
    
if(margin2 %in% c("GAi")){

sigma2 <- ifelse(sigma2 < 0.006, 0.006, sigma2) # related to gamma function
eta2   <- ifelse(eta2 < 0.0000001, 0.0000001, eta2)

 pdf2  <-  dgamma(y2, shape = 1/sigma2, scale = eta2 * sigma2)       
    p2 <-  pgamma(y2, shape = 1/sigma2, scale = eta2 * sigma2)
    
    }    
    


if(margin2 == "DAGUM"){

pdf2 <- sqrt(sigma2)*nu/y2*( ((y2/exp(eta2))^(sqrt(sigma2)*nu)) /  ( (y2/exp(eta2))^sqrt(sigma2) + 1 )^(nu+1) )            

    p2  <- ( 1 + (y2/exp(eta2))^-sqrt(sigma2) )^-nu 


}

if(margin2 == "SM"){

pdf2 <- sqrt(sigma2)*nu*y2^(sqrt(sigma2)-1)*(exp(eta2)^sqrt(sigma2)*(1+(y2/exp(eta2))^sqrt(sigma2))^(nu+1) )^-1            

    p2  <- 1 - (1+(y2/exp(eta2))^sqrt(sigma2))^-nu 


}



if(margin2 == "BE"){

pdf2 <- dbeta(y2, shape1 = plogis(eta2) * (1 - sigma2)/(sigma2), shape2 = (1-plogis(eta2))*(1 - sigma2)/(sigma2))          

    p2  <- pbeta(y2, shape1 = plogis(eta2) * (1 - sigma2)/(sigma2), shape2 = (1-plogis(eta2))*(1 - sigma2)/(sigma2))


}










epsilon <- 0.0000001 
max.p   <- 0.9999999

  #pdf2 <- ifelse(pdf2 < epsilon, epsilon, pdf2 )

  p2 <- ifelse(p2 > max.p, max.p, p2) 
  p2 <- ifelse(p2 < epsilon, epsilon, p2) 


ifef <- function(dv){

epsilon <- 0.0000001 
dv <- ifelse(is.na(dv), epsilon, dv ) 
#dv <- ifelse(dv == Inf ,  8.218407e+307, dv )
#dv <- ifelse(dv == -Inf ,  8.218407e+307, dv )
dv <- ifelse(dv == Inf ,  8.218407e+20, dv )
dv <- ifelse(dv == -Inf ,  -8.218407e+20, dv )
dv

}

# for safety

pdf2 = ifef(pdf2)
p2   = ifef(p2)





list(pdf2 = pdf2, p2 = p2)     


}




     























