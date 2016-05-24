distrHs <- function(y2, eta2, sigma2, sigma2.st, nu, nu.st, margin2, naive = FALSE){


p2 <- derp2.dersigma.st <- derp2.dereta2 <- der2p2.dereta2eta2 <- der2p2.dersigma2.st2 <- der2p2.dereta2dersigma2.st <- 1

der2pdf2.dereta2dernu.st    = 1
der2pdf2.sigma2.st2dernu.st = 1
derpdf2.dernu.st            = 1
der2pdf2.dernu.st2          = 1
derp2.nu.st                 = 1
der2p2.dernu.st2            = 1
der2p2.dereta2dernu.st      = 1
der2p2.dersigma2.stdernu.st = 1


cont2par <- c("WEI","iG","LO","rGU","GU","GA","GAi","BE","FISK") # "N" escluded for a reason
cont3par <- c("DAGUM", "SM")

# library(Deriv); library(numDeriv)
# expr <- expression(  )
# Simplify( D(D(expr, "mu2"),"sigma2") )
# func0 <- function(mu2){   }
# grad(func0 , mu2)
############################################################################
# remember that eta2 will have to disappear if we change default link on mu2
# this only applies to cases in which mu2 must be positive
# otherwise things are fine
############################################################################

#######################################################################



if(margin2 %in% c("N","LN")){

  pdf2          <- dnorm(y2, mean=eta2, sd = sqrt(sigma2))

derpdf2.dereta2 <- (1/(sqrt(2 * pi * sigma2))) * (exp(-0.5 * (y2 - eta2)^2/sigma2) * (0.5 * (2 * (y2 - eta2))/sigma2))    

derpdf2.sigma2 <-  -((0.5 - 0.5 * ((y2 - eta2)^2/sigma2)) * exp(-(0.5 * 
    ((y2 - eta2)^2/sigma2)))/(sigma2 * sqrt(2 * (pi * sigma2))))
                   
dersigma2.dersigma2.st <- exp(sigma2.st)   

derpdf2.dersigma2.st <- derpdf2.sigma2 * dersigma2.dersigma2.st   # left here for a reason

der2pdf2.dereta2 <- ((y2 - eta2)^2/sigma2 - 1) * exp(-(0.5 * ((y2 - eta2)^2/sigma2)))/(sigma2 * 
    sqrt(2 * (pi * sigma2)))
                     
                     
der2pdf2.dersigma2.st2 <- (((0.25 * ((y2 - eta2)^2/exp(sigma2.st)) - 1.5) * 
    (y2 - eta2)^2/exp(sigma2.st)^2 + 3 * (pi^2 * exp(sigma2.st)/(2 * 
    (pi * exp(sigma2.st)))^2)) * exp(sigma2.st) + 0.5 * ((y2 - 
    eta2)^2/exp(sigma2.st)) - 0.5) * exp(-(0.5 * ((y2 - eta2)^2/exp(sigma2.st))))/sqrt(2 * 
    (pi * exp(sigma2.st)))
                          

der2pdf2.dereta2dersigma2.st <-  -((1.5 - 0.5 * ((y2 - eta2)^2/exp(sigma2.st))) * exp(-(0.5 * 
    ((y2 - eta2)^2/exp(sigma2.st)))) * (y2 - eta2)/(exp(sigma2.st) * 
    sqrt(2 * (pi * exp(sigma2.st)))))
  
  
  
if(naive == FALSE){  
  
    p2          <- pnorm(y2, mean=eta2, sd = sqrt(sigma2))

derp2.dereta2    <- -pdf2
                    
derp2.dersigma.st <- 0.5*(y2-eta2)*derp2.dereta2     


der2p2.dereta2eta2 <- -((1/(sqrt(2 * pi * sigma2))) * (exp(-0.5 * (y2 - eta2)^2/sigma2) * 
                        (0.5 * (2 * (y2 - eta2))/sigma2)))


der2p2.dersigma2.st2 <-  0.5 * ((0.5 - 0.5 * ((y2 - eta2)^2/exp(sigma2.st))) * 
    exp(-(0.5 * ((y2 - eta2)^2/exp(sigma2.st)))) * (y2 - eta2)/sqrt(2 * 
    (pi * exp(sigma2.st))))

der2p2.dereta2dersigma2.st <-  (0.5 - 0.5 * ((y2 - eta2)^2/exp(sigma2.st))) * exp(-(0.5 * 
    ((y2 - eta2)^2/exp(sigma2.st))))/sqrt(2 * (pi * exp(sigma2.st))) 
    
                }


}



###################

if(margin2 == "DAGUM"){


mu2 <- dermu2.dereta2 <- der2mu2.dereta2eta2 <- exp(eta2) 
dersigma2.dersigma2.st <- exp(sigma2.st)   
dernu.dernu.st <- exp(nu.st)
 

pdf2 <- sqrt(sigma2)*nu/y2*( ((y2/mu2)^(sqrt(sigma2)*nu))/( (y2/mu2)^sqrt(sigma2) + 1 )^(nu+1) )            
  
derpdf2.dermu2 <-  -(nu * (nu * (y2/mu2)^(nu * sqrt(sigma2) - 1)/((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) - ((y2/mu2)^sqrt(sigma2) + 1)^(nu - 2 * (1 + 
    nu)) * (1 + nu) * (y2/mu2)^((1 + nu) * sqrt(sigma2) - 1)) * 
    sigma2^1/mu2^2)
    
derpdf2.sigma2 <- nu * ((0.5 * (nu * (y2/mu2)^(nu * sqrt(sigma2))/((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu)) - 0.5 * (((y2/mu2)^sqrt(sigma2) + 1)^(nu - 2 * 
    (1 + nu)) * (1 + nu) * (y2/mu2)^((1 + nu) * sqrt(sigma2)))) * 
    (log(y2) - log(mu2)) + 0.5 * ((y2/mu2)^(nu * sqrt(sigma2))/(((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) * sqrt(sigma2))))/y2    
    


derpdf2.nu <-  ((y2/mu2)^(nu * sqrt(sigma2))/((y2/mu2)^sqrt(sigma2) + 1)^(1 + 
    nu) + nu * ((log(y2) - log(mu2)) * sqrt(sigma2) * (y2/mu2)^(nu * 
    sqrt(sigma2))/((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) - ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu - 2 * (1 + nu)) * log1p((y2/mu2)^sqrt(sigma2)) * 
    (y2/mu2)^(nu * sqrt(sigma2)))) * sqrt(sigma2)/y2



der2pdf2.dermu2 <- nu * (nu * ((2 * (y2/mu2)^(nu * sqrt(sigma2) - 1) + 
    y2 * (nu * sqrt(sigma2) - 1) * (y2/mu2)^(nu * sqrt(sigma2) - 
        2)/mu2) * sqrt(sigma2)/((y2/mu2)^sqrt(sigma2) + 1)^(1 + 
    nu) - sigma2 * y2 * ((y2/mu2)^sqrt(sigma2) + 1)^(nu - 2 * 
    (1 + nu)) * (1 + nu) * (y2/mu2)^((1 + nu) * sqrt(sigma2) - 
    2)/mu2) - (((((y2/mu2)^sqrt(sigma2) + 1)^nu * (2 * (y2/mu2)^(sqrt(sigma2) - 
    1) + y2 * (sqrt(sigma2) - 1) * (y2/mu2)^(sqrt(sigma2) - 2)/mu2) * 
    sqrt(sigma2) + nu * sigma2 * y2 * ((y2/mu2)^sqrt(sigma2) + 
    1)^(nu - 1) * (y2/mu2)^(2 * (sqrt(sigma2) - 1))/mu2) * (y2/mu2)^(nu * 
    sqrt(sigma2)) + nu * sigma2 * y2 * ((y2/mu2)^sqrt(sigma2) + 
    1)^nu * (y2/mu2)^((1 + nu) * sqrt(sigma2) - 2)/mu2)/((y2/mu2)^sqrt(sigma2) + 
    1)^(2 * (1 + nu)) - 2 * (sigma2 * y2 * ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + 3 * nu - 4 * (1 + nu)) * (1 + nu) * (y2/mu2)^((2 + 
    nu) * sqrt(sigma2) - 2)/mu2)) * (1 + nu)) * sqrt(sigma2)/mu2^3


der2pdf2.dersigma22 <-  nu * (((nu * ((0.25 * (nu * (log(y2) - log(mu2)) * (y2/mu2)^(nu * 
    sqrt(sigma2))) - 0.25 * ((y2/mu2)^(nu * sqrt(sigma2))/sqrt(sigma2)))/((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) - 0.25 * (((y2/mu2)^sqrt(sigma2) + 1)^(nu - 2 * 
    (1 + nu)) * (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^((1 + 
    nu) * sqrt(sigma2)))) - (((((y2/mu2)^sqrt(sigma2) + 1)^nu * 
    (0.25 * ((log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2)) - 
        0.25 * ((y2/mu2)^sqrt(sigma2)/sqrt(sigma2))) + 0.25 * 
    (nu * ((y2/mu2)^sqrt(sigma2) + 1)^(nu - 1) * (log(y2) - log(mu2)) * 
        (y2/mu2)^(2 * sqrt(sigma2)))) * (y2/mu2)^(nu * sqrt(sigma2)) + 
    0.25 * (nu * ((y2/mu2)^sqrt(sigma2) + 1)^nu * (log(y2) - 
        log(mu2)) * (y2/mu2)^((1 + nu) * sqrt(sigma2))))/((y2/mu2)^sqrt(sigma2) + 
    1)^(2 * (1 + nu)) - 0.5 * (((y2/mu2)^sqrt(sigma2) + 1)^(1 + 
    3 * nu - 4 * (1 + nu)) * (1 + nu) * (log(y2) - log(mu2)) * 
    (y2/mu2)^((2 + nu) * sqrt(sigma2)))) * (1 + nu)) * sqrt(sigma2) + 
    0.5 * (nu * (y2/mu2)^(nu * sqrt(sigma2))/((y2/mu2)^sqrt(sigma2) + 
        1)^(1 + nu)) - 0.5 * (((y2/mu2)^sqrt(sigma2) + 1)^(nu - 
    2 * (1 + nu)) * (1 + nu) * (y2/mu2)^((1 + nu) * sqrt(sigma2)))) * 
    (log(y2) - log(mu2)) - 0.25 * ((y2/mu2)^(nu * sqrt(sigma2))/(((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) * sqrt(sigma2))))/(sigma2 * y2)




    der2pdf2.dernu2 <- (2 * ((log(y2) - log(mu2)) * sqrt(sigma2) * (y2/mu2)^(nu * sqrt(sigma2))/((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu)) + nu * ((log(y2) - log(mu2)) * (sigma2 * (log(y2) - 
    log(mu2)) * (y2/mu2)^(nu * sqrt(sigma2))/((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) - ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu - 2 * (1 + 
    nu)) * log1p((y2/mu2)^sqrt(sigma2)) * sqrt(sigma2) * (y2/mu2)^(nu * 
    sqrt(sigma2))) - ((((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) * 
    (log(y2) - log(mu2)) * sqrt(sigma2) * (y2/mu2)^(nu * sqrt(sigma2)) + 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) * log1p((y2/mu2)^sqrt(sigma2)) * 
        (y2/mu2)^(nu * sqrt(sigma2)))/((y2/mu2)^sqrt(sigma2) + 
    1)^(2 * (1 + nu)) - 2 * (((y2/mu2)^sqrt(sigma2) + 1)^(1 + 
    nu - 2 * (1 + nu)) * log1p((y2/mu2)^sqrt(sigma2)) * (y2/mu2)^(nu * 
    sqrt(sigma2)))) * log1p((y2/mu2)^sqrt(sigma2))) - 2 * (((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu - 2 * (1 + nu)) * log1p((y2/mu2)^sqrt(sigma2)) * 
    (y2/mu2)^(nu * sqrt(sigma2)))) * sqrt(sigma2)/y2
    
    
    der2pdf2.dersigma2dernu <-  ((0.5 * ((y2/mu2)^(nu * sqrt(sigma2))/((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu)) + 0.5 * (nu * ((log(y2) - log(mu2)) * sqrt(sigma2) * 
    (y2/mu2)^(nu * sqrt(sigma2))/((y2/mu2)^sqrt(sigma2) + 1)^(1 + 
    nu) - ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu - 2 * (1 + nu)) * 
    log1p((y2/mu2)^sqrt(sigma2)) * (y2/mu2)^(nu * sqrt(sigma2)))))/sqrt(sigma2) + 
    (log(y2) - log(mu2)) * (nu * (((((y2/mu2)^sqrt(sigma2) + 
        1)^(nu - 2 * (1 + nu)) * (1 + nu) * (y2/mu2)^((1 + nu) * 
        sqrt(sigma2)) - 0.5 * (nu * ((y2/mu2)^sqrt(sigma2) + 
        1)^(1 + nu - 2 * (1 + nu)) * (y2/mu2)^(nu * sqrt(sigma2)))) * 
        log1p((y2/mu2)^sqrt(sigma2))/sqrt(sigma2) + (0.5 * ((y2/mu2)^(nu * 
        sqrt(sigma2))/sqrt(sigma2)) + 0.5 * (nu * (log(y2) - 
        log(mu2)) * (y2/mu2)^(nu * sqrt(sigma2))))/((y2/mu2)^sqrt(sigma2) + 
        1)^(1 + nu) - ((0.5 * (((y2/mu2)^sqrt(sigma2) + 1)^nu * 
        (1 + nu) * log1p((y2/mu2)^sqrt(sigma2)) * (y2/mu2)^sqrt(sigma2)) + 
        0.5 * (((y2/mu2)^sqrt(sigma2) + 1)^nu * (y2/mu2)^sqrt(sigma2))) * 
        (y2/mu2)^(nu * sqrt(sigma2))/sqrt(sigma2) + 0.5 * (((y2/mu2)^sqrt(sigma2) + 
        1)^nu * (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^((1 + 
        nu) * sqrt(sigma2))))/((y2/mu2)^sqrt(sigma2) + 1)^(2 * 
        (1 + nu))) * sqrt(sigma2) + 0.5 * ((y2/mu2)^(nu * sqrt(sigma2))/((y2/mu2)^sqrt(sigma2) + 
        1)^(1 + nu))) - 0.5 * (((y2/mu2)^sqrt(sigma2) + 1)^(nu - 
        2 * (1 + nu)) * (1 + nu) * (y2/mu2)^((1 + nu) * sqrt(sigma2)))))/y2
 
 
 
 der2pdf2.mu2dernu <- -((nu * ((2 * (((y2/mu2)^sqrt(sigma2) + 1)^(nu - 2 * 
    (1 + nu)) * (1 + nu) * (y2/mu2)^((1 + nu) * sqrt(sigma2) - 
    1)) - nu * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu - 2 * (1 + 
    nu)) * (y2/mu2)^(nu * sqrt(sigma2) - 1)) * log1p((y2/mu2)^sqrt(sigma2)) * 
    sqrt(sigma2) + (nu * sigma2 * (log(y2) - log(mu2)) * (y2/mu2)^(nu * 
    sqrt(sigma2) - 1) + sqrt(sigma2) * (y2/mu2)^(nu * sqrt(sigma2) - 
    1))/((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) - ((((y2/mu2)^sqrt(sigma2) + 
    1)^nu * (1 + nu) * log1p((y2/mu2)^sqrt(sigma2)) * (y2/mu2)^(sqrt(sigma2) - 
    1) + ((y2/mu2)^sqrt(sigma2) + 1)^nu * (y2/mu2)^(sqrt(sigma2) - 
    1)) * sqrt(sigma2) * (y2/mu2)^(nu * sqrt(sigma2)) + sigma2 * 
    ((y2/mu2)^sqrt(sigma2) + 1)^nu * (1 + nu) * (log(y2) - log(mu2)) * 
    (y2/mu2)^((1 + nu) * sqrt(sigma2) - 1))/((y2/mu2)^sqrt(sigma2) + 
    1)^(2 * (1 + nu))) * sqrt(sigma2) + sigma2 * (nu * (y2/mu2)^(nu * 
    sqrt(sigma2) - 1)/((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) - 
    ((y2/mu2)^sqrt(sigma2) + 1)^(nu - 2 * (1 + nu)) * (1 + nu) * 
        (y2/mu2)^((1 + nu) * sqrt(sigma2) - 1)))/mu2^2)
    
    
 
 der2pdf2.mu2dersigma2 <- -(nu * ((nu * ((0.5 * ((y2/mu2)^(nu * sqrt(sigma2) - 
    1)/sqrt(sigma2)) + 0.5 * (nu * (log(y2) - log(mu2)) * (y2/mu2)^(nu * 
    sqrt(sigma2) - 1)))/((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) - 
    0.5 * (((y2/mu2)^sqrt(sigma2) + 1)^(nu - 2 * (1 + nu)) * 
        (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^((1 + nu) * 
        sqrt(sigma2) - 1))) - (((((y2/mu2)^sqrt(sigma2) + 1)^nu * 
    (0.5 * ((log(y2) - log(mu2)) * (y2/mu2)^(sqrt(sigma2) - 1)) + 
        0.5 * ((y2/mu2)^(sqrt(sigma2) - 1)/sqrt(sigma2))) + 0.5 * 
    (nu * ((y2/mu2)^sqrt(sigma2) + 1)^(nu - 1) * (log(y2) - log(mu2)) * 
        (y2/mu2)^(2 * sqrt(sigma2) - 1))) * (y2/mu2)^(nu * sqrt(sigma2)) + 
    0.5 * (nu * ((y2/mu2)^sqrt(sigma2) + 1)^nu * (log(y2) - log(mu2)) * 
        (y2/mu2)^((1 + nu) * sqrt(sigma2) - 1)))/((y2/mu2)^sqrt(sigma2) + 
    1)^(2 * (1 + nu)) - ((y2/mu2)^sqrt(sigma2) + 1)^(1 + 3 * 
    nu - 4 * (1 + nu)) * (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^((2 + 
    nu) * sqrt(sigma2) - 1)) * (1 + nu)) * sqrt(sigma2) + 0.5 * 
    (nu * (y2/mu2)^(nu * sqrt(sigma2) - 1)/((y2/mu2)^sqrt(sigma2) + 
        1)^(1 + nu) - ((y2/mu2)^sqrt(sigma2) + 1)^(nu - 2 * (1 + 
        nu)) * (1 + nu) * (y2/mu2)^((1 + nu) * sqrt(sigma2) - 
        1)))/mu2^2)
        
        
        
if(naive == FALSE){   
 
 
    p2  <- ( 1 + (y2/mu2)^-sqrt(sigma2) )^-nu 


derp2.dermu2 <-  -((1 + (y2/mu2)^-sqrt(sigma2))^-(nu + 1) * (nu * ((y2/mu2)^-(sqrt(sigma2) + 
    1) * (sqrt(sigma2) * (y2/mu2^2)))))
               
derp2.dersigma2 <- 0.5 * (nu * (log(y2) - log(mu2))/((1 + 1/(y2/mu2)^sqrt(sigma2))^(1 + 
    nu) * sqrt(sigma2) * (y2/mu2)^sqrt(sigma2)))        
    

derp2.dernu <- -(log1p(1/(y2/mu2)^sqrt(sigma2))/(1 + 1/(y2/mu2)^sqrt(sigma2))^nu)



der2p2.dermu22 <-  -(nu * y2 * (sqrt(sigma2) * (y2 * (1 + sqrt(sigma2))/(mu2 * 
    (y2/mu2)^(2 + sqrt(sigma2))) - 2/(y2/mu2)^(1 + sqrt(sigma2)))/(1 + 
    1/(y2/mu2)^sqrt(sigma2))^(1 + nu) - sigma2 * y2 * (1 + nu)/(mu2 * 
    (1 + 1/(y2/mu2)^sqrt(sigma2))^(2 + nu) * (y2/mu2)^(2 * (1 + 
    sqrt(sigma2)))))/mu2^3)



der2p2.dersigma22 <-  nu * (0.25 * ((1 + nu) * (log(y2) - log(mu2))/((1 + 1/(y2/mu2)^sqrt(sigma2))^(2 + 
    nu) * (y2/mu2)^(2 * sqrt(sigma2)))) - (0.25 * ((log(y2) - 
    log(mu2))/(y2/mu2)^sqrt(sigma2)) + 0.25/(sqrt(sigma2) * (y2/mu2)^sqrt(sigma2)))/(1 + 
    1/(y2/mu2)^sqrt(sigma2))^(1 + nu)) * (log(y2) - log(mu2))/sigma2
   
   
   
der2p2.dernu2 <- log1p(1/(y2/mu2)^sqrt(sigma2))^2/(1 + 1/(y2/mu2)^sqrt(sigma2))^nu
    
    
        
der2p2.dersigma2dernu <- (0.5/((1 + 1/(y2/mu2)^sqrt(sigma2))^(1 + nu) * (y2/mu2)^sqrt(sigma2)) - 
    0.5 * (nu * log1p(1/(y2/mu2)^sqrt(sigma2))/((1 + 1/(y2/mu2)^sqrt(sigma2))^(1 + 
        nu) * (y2/mu2)^sqrt(sigma2)))) * (log(y2) - log(mu2))/sqrt(sigma2)
    
    

der2p2.dermu2dernu <-  -(y2 * (1/((1 + 1/(y2/mu2)^sqrt(sigma2))^(1 + nu) * 
    (y2/mu2)^(1 + sqrt(sigma2))) - nu * log1p(1/(y2/mu2)^sqrt(sigma2))/((1 + 
    1/(y2/mu2)^sqrt(sigma2))^(1 + nu) * (y2/mu2)^(1 + sqrt(sigma2)))) * 
    sqrt(sigma2)/mu2^2)


der2p2.derdermu2sigma2 <-  -(nu * y2 * ((0.5/(sqrt(sigma2) * (y2/mu2)^(1 + sqrt(sigma2))) - 
    0.5 * ((log(y2) - log(mu2))/(y2/mu2)^(1 + sqrt(sigma2))))/(1 + 
    1/(y2/mu2)^sqrt(sigma2))^(1 + nu) + 0.5 * ((1 + nu) * (log(y2) - 
    log(mu2))/((1 + 1/(y2/mu2)^sqrt(sigma2))^(2 + nu) * (y2/mu2)^(1 + 
    2 * sqrt(sigma2)))))/mu2^2)
    

    
                                        }


}


####



if(margin2 == "SM"){


mu2 <- dermu2.dereta2 <- der2mu2.dereta2eta2 <- exp(eta2) 
dersigma2.dersigma2.st <- exp(sigma2.st)   
dernu.dernu.st <- exp(nu.st)
 

pdf2 <- sqrt(sigma2)*nu*y2^(sqrt(sigma2)-1)*(mu2^sqrt(sigma2)*(1+(y2/mu2)^sqrt(sigma2))^(nu+1) )^-1


derpdf2.dermu2 <- -(nu * y2^(sqrt(sigma2) - 1) * (mu2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) - mu2^(sqrt(sigma2) - 2) * y2 * ((y2/mu2)^sqrt(sigma2) + 
    1)^nu * (1 + nu) * (y2/mu2)^(sqrt(sigma2) - 1)) * sigma2^1/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))^2)  

   
derpdf2.sigma2 <- nu * ((0.5 * (y2^(sqrt(sigma2) - 1) * log(y2)) + 0.5 * (y2^(sqrt(sigma2) - 
    1)/sqrt(sigma2)))/(mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu)) - y2^(sqrt(sigma2) - 1) * (0.5 * (mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) * log(mu2)) + 0.5 * 
    (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 1)^nu * (1 + 
        nu) * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2)))/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))^2) 


derpdf2.nu <- sqrt(sigma2) * (y2^(sqrt(sigma2) - 1)/(mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu)) - mu2^sqrt(sigma2) * nu * y2^(sqrt(sigma2) - 
    1) * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) * log1p((y2/mu2)^sqrt(sigma2))/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))^2) 


der2pdf2.dermu2 <- -(nu * sigma2 * (y2^(sqrt(sigma2) - 1) * (mu2^(sqrt(sigma2) - 
    2) * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) * (sqrt(sigma2) - 
    1) - y2 * ((mu2^(sqrt(sigma2) - 3) * ((y2/mu2)^sqrt(sigma2) + 
    1)^nu * (sqrt(sigma2) - 2) - mu2^(sqrt(sigma2) - 4) * nu * 
    y2 * ((y2/mu2)^sqrt(sigma2) + 1)^(nu - 1) * sqrt(sigma2) * 
    (y2/mu2)^(sqrt(sigma2) - 1)) * (y2/mu2)^(sqrt(sigma2) - 1) + 
    mu2^(sqrt(sigma2) - 3) * ((y2/mu2)^sqrt(sigma2) + 1)^nu * 
        sqrt(sigma2) * (y2/mu2)^(sqrt(sigma2) - 1) - mu2^(sqrt(sigma2) - 
    4) * y2 * ((y2/mu2)^sqrt(sigma2) + 1)^nu * (sqrt(sigma2) - 
    1) * (y2/mu2)^(sqrt(sigma2) - 2)) * (1 + nu)) - 2 * (mu2^sqrt(sigma2) * 
    y2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + 
    nu) * (mu2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) - mu2^(sqrt(sigma2) - 2) * y2 * ((y2/mu2)^sqrt(sigma2) + 
    1)^nu * (1 + nu) * (y2/mu2)^(sqrt(sigma2) - 1))^2 * sqrt(sigma2)/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))^2))/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))^2) 


der2pdf2.dersigma22 <- nu * ((0.25 * (y2^(sqrt(sigma2) - 1) * log(y2)^2/sqrt(sigma2)) + 
    0.5 * ((0.5 * (y2^(sqrt(sigma2) - 1) * log(y2)) - 0.5 * (y2^(sqrt(sigma2) - 
        1)/sqrt(sigma2)))/sigma2))/(mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu)) - ((0.5 * (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) * log(mu2)) + 0.5 * (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^nu * (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2))) * 
    (0.5 * (y2^(sqrt(sigma2) - 1)/sqrt(sigma2)) + y2^(sqrt(sigma2) - 
        1) * log(y2) - 2 * (mu2^sqrt(sigma2) * y2^(sqrt(sigma2) - 
        1) * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) * (0.5 * (mu2^sqrt(sigma2) * 
        ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) * log(mu2)) + 0.5 * 
        (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 1)^nu * 
            (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2)))/(mu2^sqrt(sigma2) * 
        ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))^2)) + y2^(sqrt(sigma2) - 
    1) * (0.5 * (((0.5 * (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^nu * log(mu2)) + 0.5 * (mu2^sqrt(sigma2) * nu * ((y2/mu2)^sqrt(sigma2) + 
    1)^(nu - 1) * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2))) * 
    (y2/mu2)^sqrt(sigma2) + 0.5 * (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^nu * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2))) * 
    (1 + nu) * (log(y2) - log(mu2))) + 0.5 * ((0.5 * (mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) * log(mu2)) + 0.5 * 
    (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 1)^nu * (1 + 
        nu) * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2))) * 
    log(mu2))))/((mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu))^2 * sqrt(sigma2))) 


der2pdf2.dernu2 <- -(log1p((y2/mu2)^sqrt(sigma2)) * (mu2^sqrt(sigma2) * y2^(sqrt(sigma2) - 
    1) * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) + mu2^sqrt(sigma2) * 
    y2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + 
    nu) + nu * y2^(sqrt(sigma2) - 1) * log1p((y2/mu2)^sqrt(sigma2)) * 
    (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) - 
        2 * (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 1)^(2 * 
            (1 + nu) - (1 + nu))))) * sqrt(sigma2)/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))^2) 
    
       
der2pdf2.dersigma2dernu <- (0.5 * (y2^(sqrt(sigma2) - 1) * log(y2)) + 0.5 * (y2^(sqrt(sigma2) - 
    1)/sqrt(sigma2)))/(mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu)) - (nu * (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) * (0.5 * (y2^(sqrt(sigma2) - 1) * log(y2)) + 
    0.5 * (y2^(sqrt(sigma2) - 1)/sqrt(sigma2))) * log1p((y2/mu2)^sqrt(sigma2)) + 
    y2^(sqrt(sigma2) - 1) * ((0.5 * (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
        1)^(1 + nu) * log(mu2)) - 2 * (0.5 * (mu2^sqrt(sigma2) * 
        ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) * log(mu2)) + 0.5 * 
        (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 1)^nu * 
            (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2)))) * 
        log1p((y2/mu2)^sqrt(sigma2)) + mu2^sqrt(sigma2) * (0.5 * 
        (((y2/mu2)^sqrt(sigma2) + 1)^nu * (1 + nu) * log1p((y2/mu2)^sqrt(sigma2)) * 
            (y2/mu2)^sqrt(sigma2)) + 0.5 * (((y2/mu2)^sqrt(sigma2) + 
        1)^nu * (y2/mu2)^sqrt(sigma2))) * (log(y2) - log(mu2)))) + 
    y2^(sqrt(sigma2) - 1) * (0.5 * (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
        1)^(1 + nu) * log(mu2)) + 0.5 * (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
        1)^nu * (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2))))/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))^2 
 
 
der2pdf2.mu2dernu <- -(sigma2 * (nu * y2^(sqrt(sigma2) - 1) * (log1p((y2/mu2)^sqrt(sigma2)) * 
    (mu2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + 
        nu) - 2 * (mu2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 
        1)^(1 + nu) - mu2^(sqrt(sigma2) - 2) * y2 * ((y2/mu2)^sqrt(sigma2) + 
        1)^nu * (1 + nu) * (y2/mu2)^(sqrt(sigma2) - 1))) - mu2^(sqrt(sigma2) - 
    2) * y2 * (((y2/mu2)^sqrt(sigma2) + 1)^nu * (1 + nu) * log1p((y2/mu2)^sqrt(sigma2)) * 
    (y2/mu2)^(sqrt(sigma2) - 1) + ((y2/mu2)^sqrt(sigma2) + 1)^nu * 
    (y2/mu2)^(sqrt(sigma2) - 1))) + y2^(sqrt(sigma2) - 1) * (mu2^(sqrt(sigma2) - 
    1) * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) - mu2^(sqrt(sigma2) - 
    2) * y2 * ((y2/mu2)^sqrt(sigma2) + 1)^nu * (1 + nu) * (y2/mu2)^(sqrt(sigma2) - 
    1)))/(mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + 
    nu))^2)
    
    
der2pdf2.mu2dersigma2 <- -(nu * ((0.5 * (y2^(sqrt(sigma2) - 1) * log(y2)) + 0.5 * (y2^(sqrt(sigma2) - 
    1)/sqrt(sigma2))) * (mu2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) - mu2^(sqrt(sigma2) - 2) * y2 * ((y2/mu2)^sqrt(sigma2) + 
    1)^nu * (1 + nu) * (y2/mu2)^(sqrt(sigma2) - 1)) + y2^(sqrt(sigma2) - 
    1) * (((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) * (0.5 * (mu2^(sqrt(sigma2) - 
    1) * log(mu2)) + 0.5 * (mu2^(sqrt(sigma2) - 1)/sqrt(sigma2))) + 
    (0.5 * (mu2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 
        1)^nu * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2)) - 
        y2 * (0.5 * (mu2^(sqrt(sigma2) - 2) * ((y2/mu2)^sqrt(sigma2) + 
            1)^nu * log(mu2) * (y2/mu2)^(sqrt(sigma2) - 1)) + 
            mu2^(sqrt(sigma2) - 2) * (((y2/mu2)^sqrt(sigma2) + 
                1)^nu * (0.5 * ((log(y2) - log(mu2)) * (y2/mu2)^(sqrt(sigma2) - 
                1)) + 0.5 * ((y2/mu2)^(sqrt(sigma2) - 1)/sqrt(sigma2))) + 
                0.5 * (nu * ((y2/mu2)^sqrt(sigma2) + 1)^(nu - 
                  1) * (log(y2) - log(mu2)) * (y2/mu2)^(2 * sqrt(sigma2) - 
                  1))))) * (1 + nu) - 2 * ((0.5 * (mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) * log(mu2)) + 0.5 * 
    (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 1)^nu * (1 + 
        nu) * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2))) * 
    (mu2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 1)^(1 + 
        nu) - mu2^(sqrt(sigma2) - 2) * y2 * ((y2/mu2)^sqrt(sigma2) + 
        1)^nu * (1 + nu) * (y2/mu2)^(sqrt(sigma2) - 1))/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))))) * sqrt(sigma2)/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))^2) 
        
        
        
if(naive == FALSE){   
 
 
p2  <-  1 - (1+(y2/mu2)^sqrt(sigma2))^-nu

derp2.dermu2 <- -(nu * y2 * sqrt(sigma2) * (y2/mu2)^(sqrt(sigma2) - 1)/(mu2^2 * 
    ((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))) 

                                    
derp2.dersigma2 <- 0.5 * (nu * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2)/(((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) * sqrt(sigma2)))            
    

derp2.dernu <- log1p((y2/mu2)^sqrt(sigma2))/((y2/mu2)^sqrt(sigma2) + 1)^nu 


der2p2.dermu22 <- -(nu * y2 * (sigma2 * y2 * (1 + nu) * (y2/mu2)^(2 * (sqrt(sigma2) - 
    1))/(mu2 * ((y2/mu2)^sqrt(sigma2) + 1)^(2 + nu)) - (2 * (y2/mu2)^(sqrt(sigma2) - 
    1) + y2 * (sqrt(sigma2) - 1) * (y2/mu2)^(sqrt(sigma2) - 2)/mu2) * 
    sqrt(sigma2)/((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu))/mu2^3)
 


der2p2.dersigma22 <- nu * ((0.25 * ((log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2)) - 
    0.25 * ((y2/mu2)^sqrt(sigma2)/sqrt(sigma2)))/((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) - 0.25 * ((1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^(2 * 
    sqrt(sigma2))/((y2/mu2)^sqrt(sigma2) + 1)^(2 + nu))) * (log(y2) - 
    log(mu2))/sigma2 
   
   
der2p2.dernu2 <- -(log1p((y2/mu2)^sqrt(sigma2))^2/((y2/mu2)^sqrt(sigma2) + 1)^nu)

  
der2p2.dersigma2dernu <- (0.5 * ((y2/mu2)^sqrt(sigma2)/((y2/mu2)^sqrt(sigma2) + 1)^(1 + 
    nu)) - 0.5 * (nu * log1p((y2/mu2)^sqrt(sigma2)) * (y2/mu2)^sqrt(sigma2)/((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu))) * (log(y2) - log(mu2))/sqrt(sigma2)

    
der2p2.dermu2dernu <- y2 * (nu * log1p((y2/mu2)^sqrt(sigma2)) * (y2/mu2)^(sqrt(sigma2) - 
    1)/((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu) - (y2/mu2)^(sqrt(sigma2) - 
    1)/((y2/mu2)^sqrt(sigma2) + 1)^(1 + nu)) * sqrt(sigma2)/mu2^2
 


der2p2.derdermu2sigma2 <- -(nu * y2 * ((0.5 * ((log(y2) - log(mu2)) * (y2/mu2)^(sqrt(sigma2) - 
    1)) + 0.5 * ((y2/mu2)^(sqrt(sigma2) - 1)/sqrt(sigma2)))/((y2/mu2)^sqrt(sigma2) + 
    1)^(1 + nu) - 0.5 * ((1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^(2 * 
    sqrt(sigma2) - 1)/((y2/mu2)^sqrt(sigma2) + 1)^(2 + nu)))/mu2^2)

 
                   }


}



###


if(margin2 == "FISK"){


mu2 <- dermu2.dereta2 <- der2mu2.dereta2eta2 <- exp(eta2) 
dersigma2.dersigma2.st <- exp(sigma2.st)  

pdf2 <- sqrt(sigma2)*y2^(sqrt(sigma2)-1) / (mu2^sqrt(sigma2)*(1+(y2/mu2)^sqrt(sigma2))^2)
  
derpdf2.dermu2 <- -(sigma2 * y2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 1) * 
    (mu2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 1) - 2 * 
        (mu2^(sqrt(sigma2) - 2) * y2 * (y2/mu2)^(sqrt(sigma2) - 
            1)))/(mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^2)^2)

                
derpdf2.sigma2 <- (0.5 * (y2^(sqrt(sigma2) - 1) * log(y2)) + 0.5 * (y2^(sqrt(sigma2) - 
    1)/sqrt(sigma2)))/(mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^2) - y2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 
    1) * (0.5 * (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1) * log(mu2)) + mu2^sqrt(sigma2) * (log(y2) - log(mu2)) * 
    (y2/mu2)^sqrt(sigma2))/(mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^2)^2


der2pdf2.dermu2 <- -(sqrt(sigma2) * (y2^(sqrt(sigma2) - 1) * (((y2/mu2)^sqrt(sigma2) + 
    1) * (mu2^(sqrt(sigma2) - 2) * ((y2/mu2)^sqrt(sigma2) + 1) * 
    sqrt(sigma2) * (sqrt(sigma2) - 1) - 2 * (mu2^(sqrt(sigma2) - 
    3) * sigma2 * y2 * (y2/mu2)^(sqrt(sigma2) - 1))) - y2 * (2 * 
    (mu2^(sqrt(sigma2) - 3) * sigma2 * ((y2/mu2)^sqrt(sigma2) + 
        1) * (y2/mu2)^(sqrt(sigma2) - 1)) - 2 * (mu2^(sqrt(sigma2) - 
    3) * (((y2/mu2)^sqrt(sigma2) + 1) * (2 * (y2/mu2)^(sqrt(sigma2) - 
    1) + y2 * (sqrt(sigma2) - 1) * (y2/mu2)^(sqrt(sigma2) - 2)/mu2) * 
    sqrt(sigma2) + sigma2 * y2 * (y2/mu2)^(2 * (sqrt(sigma2) - 
    1))/mu2)))) - 2 * (mu2^sqrt(sigma2) * sigma2 * y2^(sqrt(sigma2) - 
    1) * ((y2/mu2)^sqrt(sigma2) + 1)^4 * (mu2^(sqrt(sigma2) - 
    1) * ((y2/mu2)^sqrt(sigma2) + 1) - 2 * (mu2^(sqrt(sigma2) - 
    2) * y2 * (y2/mu2)^(sqrt(sigma2) - 1)))^2/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^2)^2))/(mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^2)^2)

            
der2pdf2.dersigma22 <- (((0.25 * (y2^(sqrt(sigma2) - 1) * log(y2)) - 0.25 * (y2^(sqrt(sigma2) - 
    1)/sqrt(sigma2))) * sqrt(sigma2) + 0.25 * y2^(sqrt(sigma2) - 
    1) + 0.25 * y2^(sqrt(sigma2) - 1)) * log(y2) - 0.25 * (y2^(sqrt(sigma2) - 
    1)/sqrt(sigma2)))/(mu2^sqrt(sigma2) * sigma2 * ((y2/mu2)^sqrt(sigma2) + 
    1)^2) - (((y2/mu2)^sqrt(sigma2) + 1) * (0.5 * (mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1) * log(mu2)) + mu2^sqrt(sigma2) * 
    (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2)) * (y2^(sqrt(sigma2) - 
    1) * log(y2) + y2^(sqrt(sigma2) - 1)/sqrt(sigma2) - 2 * (mu2^sqrt(sigma2) * 
    y2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 1)^3 * (0.5 * 
    (mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 1) * log(mu2)) + 
    mu2^sqrt(sigma2) * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2))/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^2)^2))/sqrt(sigma2) + y2^(sqrt(sigma2) - 
    1) * ((((y2/mu2)^sqrt(sigma2) + 1) * (0.25 * (mu2^sqrt(sigma2) * 
    log(mu2)) - 0.25 * (mu2^sqrt(sigma2)/sqrt(sigma2))) + 0.5 * 
    (mu2^sqrt(sigma2) * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2))) * 
    ((y2/mu2)^sqrt(sigma2) + 1) * log(mu2) + (0.5 * (mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1) * log(mu2) * (y2/mu2)^sqrt(sigma2)) + 
    2 * (mu2^sqrt(sigma2) * (((y2/mu2)^sqrt(sigma2) + 1) * (0.25 * 
        ((log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2)) - 0.25 * 
        ((y2/mu2)^sqrt(sigma2)/sqrt(sigma2))) + 0.25 * ((log(y2) - 
        log(mu2)) * (y2/mu2)^(2 * sqrt(sigma2)))))) * (log(y2) - 
    log(mu2))) * sqrt(sigma2)/sigma2)/(mu2^sqrt(sigma2) * ((y2/mu2)^sqrt(sigma2) + 
    1)^2)^2

         
der2pdf2.mu2dersigma2 <- -(((((y2/mu2)^sqrt(sigma2) + 1) * (0.5 * (y2^(sqrt(sigma2) - 
    1) * log(y2)) + 0.5 * (y2^(sqrt(sigma2) - 1)/sqrt(sigma2))) * 
    (mu2^(sqrt(sigma2) - 1) * ((y2/mu2)^sqrt(sigma2) + 1) - 2 * 
        (mu2^(sqrt(sigma2) - 2) * y2 * (y2/mu2)^(sqrt(sigma2) - 
            1))) + y2^(sqrt(sigma2) - 1) * ((((y2/mu2)^sqrt(sigma2) + 
    1) * (0.5 * (mu2^(sqrt(sigma2) - 1) * log(mu2)) + 0.5 * (mu2^(sqrt(sigma2) - 
    1)/sqrt(sigma2))) + mu2^(sqrt(sigma2) - 1) * (log(y2) - log(mu2)) * 
    (y2/mu2)^sqrt(sigma2)) * ((y2/mu2)^sqrt(sigma2) + 1) - y2 * 
    (2 * (mu2^(sqrt(sigma2) - 2) * (((y2/mu2)^sqrt(sigma2) + 
        1) * (0.5 * ((log(y2) - log(mu2)) * (y2/mu2)^(sqrt(sigma2) - 
        1)) + 0.5 * ((y2/mu2)^(sqrt(sigma2) - 1)/sqrt(sigma2))) + 
        0.5 * ((log(y2) - log(mu2)) * (y2/mu2)^(2 * sqrt(sigma2) - 
            1)))) + mu2^(sqrt(sigma2) - 2) * ((y2/mu2)^sqrt(sigma2) + 
        1) * log(mu2) * (y2/mu2)^(sqrt(sigma2) - 1)))) * sqrt(sigma2) - 
    2 * (mu2^sqrt(sigma2) * sigma2 * y2^(sqrt(sigma2) - 1) * 
        ((y2/mu2)^sqrt(sigma2) + 1)^4 * (0.5 * (mu2^sqrt(sigma2) * 
        ((y2/mu2)^sqrt(sigma2) + 1) * log(mu2)) + mu2^sqrt(sigma2) * 
        (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2)) * (mu2^(sqrt(sigma2) - 
        1) * ((y2/mu2)^sqrt(sigma2) + 1) - 2 * (mu2^(sqrt(sigma2) - 
        2) * y2 * (y2/mu2)^(sqrt(sigma2) - 1)))/((mu2^sqrt(sigma2) * 
        ((y2/mu2)^sqrt(sigma2) + 1)^2)^2 * sqrt(sigma2))))/(mu2^sqrt(sigma2) * 
    ((y2/mu2)^sqrt(sigma2) + 1)^2)^2)


 
  
 
if(naive == FALSE){   
 
    p2  <- 1/(1+(y2/mu2)^-sqrt(sigma2))
    
    derp2.dermu2 <- -(y2 * sqrt(sigma2)/(mu2^2 * (1 + 1/(y2/mu2)^sqrt(sigma2))^2 * 
    (y2/mu2)^(1 + sqrt(sigma2))))


                 
derp2.dersigma2 <-  0.5 * ((log(y2) - log(mu2))/((1 + 1/(y2/mu2)^sqrt(sigma2))^2 * 
    sqrt(sigma2) * (y2/mu2)^sqrt(sigma2)))
             
    

der2p2.dermu22 <- -(y2 * (sqrt(sigma2) * (y2 * (1 + sqrt(sigma2))/(mu2 * (y2/mu2)^(2 + 
    sqrt(sigma2))) - 2/(y2/mu2)^(1 + sqrt(sigma2))) - 2 * (sigma2 * 
    y2/(mu2 * (1 + 1/(y2/mu2)^sqrt(sigma2)) * (y2/mu2)^(2 * (1 + 
    sqrt(sigma2))))))/(mu2^3 * (1 + 1/(y2/mu2)^sqrt(sigma2))^2))


der2p2.dersigma22 <- ((0.5/((1 + 1/(y2/mu2)^sqrt(sigma2)) * (y2/mu2)^(2 * sqrt(sigma2))) - 
    0.25/(y2/mu2)^sqrt(sigma2)) * (log(y2) - log(mu2)) - 0.25/(sqrt(sigma2) * 
    (y2/mu2)^sqrt(sigma2))) * (log(y2) - log(mu2))/(sigma2 * 
    (1 + 1/(y2/mu2)^sqrt(sigma2))^2)



der2p2.derdermu2sigma2 <- -(y2 * ((1/((1 + 1/(y2/mu2)^sqrt(sigma2)) * (y2/mu2)^(1 + 2 * 
    sqrt(sigma2))) - 0.5/(y2/mu2)^(1 + sqrt(sigma2))) * (log(y2) - 
    log(mu2)) + 0.5/(sqrt(sigma2) * (y2/mu2)^(1 + sqrt(sigma2))))/(mu2^2 * 
    (1 + 1/(y2/mu2)^sqrt(sigma2))^2))

   
                                             
                                             
                                             
                   }



}





####



if(margin2 == "WEI"){


mu2 <- dermu2.dereta2 <- der2mu2.dereta2eta2 <- exp(eta2) 
dersigma2.dersigma2.st <- exp(sigma2.st)  

pdf2 <- sqrt(sigma2)/mu2*(y2/mu2)^(sqrt(sigma2)-1) * exp(-(y2/mu2)^sqrt(sigma2))  
  
derpdf2.dermu2 <-   exp(-(y2/mu2)^sqrt(sigma2)) * (sigma2 * y2 * (y2/mu2)^(2 * 
    (sqrt(sigma2) - 1))/mu2 - ((y2/mu2)^(sqrt(sigma2) - 1) + 
    y2 * (sqrt(sigma2) - 1) * (y2/mu2)^(sqrt(sigma2) - 2)/mu2) * 
    sqrt(sigma2))/mu2^2
  
                  
derpdf2.sigma2 <- ((0.5 * (y2/mu2)^(sqrt(sigma2) - 1) - 0.5 * (y2/mu2)^(2 * sqrt(sigma2) - 
    1)) * (log(y2) - log(mu2)) + 0.5 * ((y2/mu2)^(sqrt(sigma2) - 
    1)/sqrt(sigma2))) * exp(-(y2/mu2)^sqrt(sigma2))/mu2


der2pdf2.dermu2 <- ((2 * (y2/mu2)^(sqrt(sigma2) - 1) + y2 * (4 * (y2/mu2)^(sqrt(sigma2) - 
    2) + y2 * (sqrt(sigma2) - 2) * (y2/mu2)^(sqrt(sigma2) - 3)/mu2) * 
    (sqrt(sigma2) - 1)/mu2) * sqrt(sigma2) + y2 * ((sigma2 * 
    y2 * (y2/mu2)^(2 * (sqrt(sigma2) - 1))/mu2 - (2 * (y2/mu2)^(sqrt(sigma2) - 
    1) + y2 * (sqrt(sigma2) - 1) * (y2/mu2)^(sqrt(sigma2) - 2)/mu2) * 
    sqrt(sigma2)) * sqrt(sigma2) * (y2/mu2)^(sqrt(sigma2) - 1) - 
    2 * (sigma2 * ((y2/mu2)^(sqrt(sigma2) - 1) + y2 * (sqrt(sigma2) - 
        1) * (y2/mu2)^(sqrt(sigma2) - 2)/mu2) * (y2/mu2)^(sqrt(sigma2) - 
        1)))/mu2) * exp(-(y2/mu2)^sqrt(sigma2))/mu2^3
            
            
 
der2pdf2.dersigma22 <- ((((0.25 * ((log(y2) - log(mu2)) * (y2/mu2)^(sqrt(sigma2) - 1)) - 
    0.25 * ((y2/mu2)^(sqrt(sigma2) - 1)/sqrt(sigma2))) * sqrt(sigma2) + 
    0.25 * (y2/mu2)^(sqrt(sigma2) - 1) + 0.25 * (y2/mu2)^(sqrt(sigma2) - 
    1)) * (log(y2) - log(mu2)) - 0.25 * ((y2/mu2)^(sqrt(sigma2) - 
    1)/sqrt(sigma2)))/sigma2 - (((0.25 * (y2/mu2)^sqrt(sigma2) - 
    0.25 * (y2/mu2)^(2 * sqrt(sigma2))) * (log(y2) - log(mu2)) - 
    0.25 * ((y2/mu2)^sqrt(sigma2)/sqrt(sigma2))) * sqrt(sigma2) * 
    (y2/mu2)^(sqrt(sigma2) - 1)/sigma2 + (0.5 * ((log(y2) - log(mu2)) * 
    (y2/mu2)^(sqrt(sigma2) - 1)) + 0.5 * ((y2/mu2)^(sqrt(sigma2) - 
    1)/sqrt(sigma2))) * (y2/mu2)^sqrt(sigma2)/sqrt(sigma2)) * 
    (log(y2) - log(mu2))) * exp(-(y2/mu2)^sqrt(sigma2))/mu2
       
       
       
der2pdf2.mu2dersigma2 <- exp(-(y2/mu2)^sqrt(sigma2)) * (y2 * ((((0.5 * (y2/mu2)^(sqrt(sigma2) - 
    1) - 0.5 * (y2/mu2)^(2 * sqrt(sigma2) - 1)) * (log(y2) - 
    log(mu2)) + 0.5 * ((y2/mu2)^(sqrt(sigma2) - 1)/sqrt(sigma2))) * 
    (y2/mu2)^(sqrt(sigma2) - 1) + (0.5 * ((log(y2) - log(mu2)) * 
    (y2/mu2)^(sqrt(sigma2) - 1)) + 0.5 * ((y2/mu2)^(sqrt(sigma2) - 
    1)/sqrt(sigma2))) * (y2/mu2)^(sqrt(sigma2) - 1)) * sqrt(sigma2) - 
    ((0.5 * ((log(y2) - log(mu2)) * (y2/mu2)^(sqrt(sigma2) - 
        2)) + 0.5 * ((y2/mu2)^(sqrt(sigma2) - 2)/sqrt(sigma2))) * 
        (sqrt(sigma2) - 1) + 0.5 * (y2/mu2)^(sqrt(sigma2) - 2)))/mu2 - 
    ((0.5 * (y2/mu2)^(sqrt(sigma2) - 1) - 0.5 * (((y2/mu2)^(sqrt(sigma2) - 
        1) + y2 * (sqrt(sigma2) - 1) * (y2/mu2)^(sqrt(sigma2) - 
        2)/mu2) * (y2/mu2)^sqrt(sigma2))) * (log(y2) - log(mu2)) + 
        0.5 * ((y2/mu2)^(sqrt(sigma2) - 1)/sqrt(sigma2))))/mu2^2

 
  
 
if(naive == FALSE){   
 
    p2  <-  1-exp(-(y2/mu2)^sqrt(sigma2)) 
    
    derp2.dermu2 <- -(exp(-(y2/mu2)^sqrt(sigma2)) * ((y2/mu2)^(sqrt(sigma2) - 1) * 
    (sqrt(sigma2) * (y2/mu2^2))))

                 
derp2.dersigma2 <- 0.5 * (exp(-(y2/mu2)^sqrt(sigma2)) * (log(y2) - log(mu2)) * (y2/mu2)^sqrt(sigma2)/sqrt(sigma2))         
    

der2p2.dermu22 <-  -(y2 * exp(-(y2/mu2)^sqrt(sigma2)) * (sigma2 * y2 * 
    (y2/mu2)^(2 * (sqrt(sigma2) - 1))/mu2 - (2 * (y2/mu2)^(sqrt(sigma2) - 
    1) + y2 * (sqrt(sigma2) - 1) * (y2/mu2)^(sqrt(sigma2) - 2)/mu2) * 
    sqrt(sigma2))/mu2^3)



der2p2.dersigma22 <-  ((0.25 * (y2/mu2)^sqrt(sigma2) - 0.25 * (y2/mu2)^(2 * sqrt(sigma2))) * 
    (log(y2) - log(mu2)) - 0.25 * ((y2/mu2)^sqrt(sigma2)/sqrt(sigma2))) * 
    exp(-(y2/mu2)^sqrt(sigma2)) * (log(y2) - log(mu2))/sigma2




der2p2.derdermu2sigma2 <-   -(y2 * ((0.5 * (y2/mu2)^(sqrt(sigma2) - 1) - 0.5 * 
    (y2/mu2)^(2 * sqrt(sigma2) - 1)) * (log(y2) - log(mu2)) + 
    0.5 * ((y2/mu2)^(sqrt(sigma2) - 1)/sqrt(sigma2))) * exp(-(y2/mu2)^sqrt(sigma2))/mu2^2)
   
                                             }


}





###########

#a <- mu * (1 - sigma^2)/(sigma^2)
#b <- a * (1 - mu)/mu


if(margin2 == "BE"){



#########################

mu2 <- plogis(eta2)

dermu2.dereta2 <- (1 - exp(eta2)/(1 + exp(eta2))) * exp(eta2)/(1 + exp(eta2))

der2mu2.dereta2eta2 <-  (1 - (3 - 2 * (exp(eta2)/(1 + exp(eta2)))) * exp(eta2)/(1 + exp(eta2))) * 
                        exp(eta2)/(1 + exp(eta2))
                        
sigma2 <- plogis(sigma2.st)  

dersigma2.dersigma2.st <- (1 - exp(sigma2.st)/(1 + exp(sigma2.st))) * exp(sigma2.st)/(1 + exp(sigma2.st))

#########################







pdf2 <- dbeta(y2, shape1 = mu2 * (1 - sigma2)/(sigma2), shape2 = (1-mu2)*(1 - sigma2)/(sigma2))
  
  
derpdf2.dermu2 <- ((-1 + sigma2)* (1 - y2)^(-1 + ((-1 + mu2)* (-1 + sigma2))/sigma2)*
    y2^(-1 + mu2 * (-1 + 1/sigma2))* (log(1 - y2) - log(y2) + 
     psigamma(mu2* (-1 + 1/sigma2)) - 
     psigamma(((-1 + mu2)* (-1 + sigma2))/sigma2)))/(sigma2* beta(
    mu2 *(-1 + 1/sigma2), ((-1 + mu2)* (-1 + sigma2))/sigma2))
  
            
            
derpdf2.sigma2 <- -((1 - y2)^(-1 + ((-1 + mu2)* (-1 + sigma2))/sigma2)*
     y2^(-1 + mu2 * (-1 + 1/sigma2)) *(log(1 - y2) - mu2* log(1 - y2) + 
      mu2 * log(y2) + psigamma(-1 + 1/sigma2) - 
      mu2* psigamma(mu2 *(-1 + 1/sigma2)) - 
      psigamma(((-1 + mu2) *(-1 + sigma2))/sigma2) + 
      mu2* psigamma(((-1 + mu2)* (-1 + sigma2))/
        sigma2)))/(sigma2^2 *beta(
     mu2* (-1 + 1/sigma2), ((-1 + mu2)* (-1 + sigma2))/sigma2))



der2pdf2.dermu2 <- ((-1 + sigma2)^2* (1 - y2)^(-1 + ((-1 + mu2)* (-1 + sigma2))/sigma2)*
    y2^(-1 + 
    mu2* (-1 + 1/sigma2)) *(log(1 - y2)^2 - 2* log(1 - y2)* log(y2) + 
     log(y2)^2 + psigamma(mu2* (-1 + 1/sigma2))^2 + 
     2 *psigamma( 
       mu2* (-1 + 1/sigma2))* (log(1 - y2) - log(y2) - 
        psigamma(((-1 + mu2) *(-1 + sigma2))/sigma2)) - 
     2 *(log(1 - y2) - log(y2))* psigamma(
       ((-1 + mu2)* (-1 + sigma2))/sigma2) + 
     psigamma(((-1 + mu2)* (-1 + sigma2))/sigma2)^2 - 
     psigamma(mu2 *(-1 + 1/sigma2),1) - 
     psigamma( ((-1 + mu2)* (-1 + sigma2))/sigma2,1)))/(sigma2^2* beta(
    mu2* (-1 + 1/sigma2), ((-1 + mu2) *(-1 + sigma2))/sigma2))
            
         
         
 
der2pdf2.dersigma22 <- -((1 - y2)^(-1 + ((-1 + mu2) *(-1 + sigma2))/sigma2)*
     y2^(-1 + 
     mu2* (-1 + 1/sigma2))* (-2 *sigma2* log(1 - y2) + 
      2* mu2* sigma2* log(1 - y2) - log(1 - y2)^2 + 2* mu2* log(1 - y2)^2 -
       mu2^2* log(1 - y2)^2 - 2* mu2* sigma2* log(y2) - 
      2 *mu2* log(1 - y2)* log(y2) + 2* mu2^2 *log(1 - y2)* log(y2) - 
      mu2^2* log(y2)^2 - psigamma(-1 + 1/sigma2)^2 - 
      mu2^2 *psigamma(mu2 *(-1 + 1/sigma2))^2 + 
      2 *sigma2* psigamma(((-1 + mu2) *(-1 + sigma2))/sigma2) - 
      2 *mu2* sigma2* psigamma( ((-1 + mu2) *(-1 + sigma2))/sigma2) + 
      2* log(1 - y2)* psigamma(((-1 + mu2) *(-1 + sigma2))/sigma2) - 
      4 *mu2* log(1 - y2)* psigamma(((-1 + mu2)* (-1 + sigma2))/
        sigma2) + 
      2* mu2^2* log(1 - y2)* psigamma(((-1 + mu2)* (-1 + sigma2))/
        sigma2) + 
      2 *mu2 *log(y2)* psigamma(((-1 + mu2)* (-1 + sigma2))/sigma2) - 
      2 *mu2^2* log(y2)* psigamma( ((-1 + mu2)* (-1 + sigma2))/
        sigma2) - psigamma( ((-1 + mu2)* (-1 + sigma2))/sigma2)^2 + 
      2* mu2* psigamma(((-1 + mu2)* (-1 + sigma2))/sigma2)^2 - 
      mu2^2 *psigamma( ((-1 + mu2)* (-1 + sigma2))/sigma2)^2 + 
      2* mu2* psigamma( 
        mu2* (-1 + 1/sigma2))* (sigma2 - (-1 + mu2)* log(1 - y2) + 
         mu2 *log(y2) + (-1 + mu2)* psigamma(
           ((-1 + mu2)* (-1 + sigma2))/sigma2)) - 
      2* psigamma(-1 + 1/sigma2)* (sigma2 + log(1 - y2) - mu2* log(1 - y2) + 
         mu2 *log(y2) - 
         mu2* psigamma( mu2* (-1 + 1/sigma2)) + (-1 + mu2)* psigamma(((-1 + mu2)* (-1 + sigma2))/sigma2)) - 
      psigamma(-1 + 1/sigma2,1) + 
      mu2^2* psigamma( mu2* (-1 + 1/sigma2),1) + 
      psigamma( ((-1 + mu2)* (-1 + sigma2))/sigma2,1) - 
      2* mu2* psigamma(((-1 + mu2)* (-1 + sigma2))/sigma2,1) + 
      mu2^2 *psigamma(((-1 + mu2)* (-1 + sigma2))/
        sigma2,1)))/(sigma2^4* beta(
     mu2 *(-1 + 1/sigma2), ((-1 + mu2)* (-1 + sigma2))/sigma2))
       
       
       
       
der2pdf2.mu2dersigma2 <- ((1 - y2)^(-1 + ((-1 + mu2) *(-1 + sigma2))/sigma2)*
    y2^(-1 + 
    mu2* (-1 + 1/
       sigma2))* (-(-1 + sigma2)* sigma2* (log(1 - y2) - log(y2) + 
        psigamma(mu2* (-1 + 1/sigma2)) - 
        psigamma(((-1 + mu2)* (-1 + sigma2))/sigma2)) + 
     sigma2^2* (log(1 - y2) - log(y2) + 
        psigamma(mu2 *(-1 + 1/sigma2)) - 
        psigamma(((-1 + mu2)* (-1 + sigma2))/sigma2)) + (-1 + 
        mu2)* (-1 + sigma2)* log(
       1 - y2)* (log(1 - y2) - log(y2) + 
        psigamma(mu2* (-1 + 1/sigma2)) - 
        psigamma(((-1 + mu2)* (-1 + sigma2))/sigma2)) + 
     mu2* (-1 + sigma2)* log(
       y2)* (-log(1 - y2) + log(y2) - 
        psigamma(mu2 *(-1 + 1/sigma2)) + 
        psigamma(((-1 + mu2)* (-1 + sigma2))/sigma2)) + (-1 + 
        sigma2)* (log(1 - y2) - log(y2) + 
        psigamma(mu2* (-1 + 1/sigma2)) - 
        psigamma(((-1 + mu2)* (-1 + sigma2))/
         sigma2))* (-psigamma(-1 + 1/sigma2) + 
        mu2* psigamma(mu2 *(-1 + 1/sigma2)) - (-1 + mu2)* psigamma(
          ((-1 + mu2)* (-1 + sigma2))/sigma2)) - (-1 + 
        sigma2)* (mu2* psigamma( 
          mu2 *(-1 + 1/sigma2),1) + (-1 + mu2)* psigamma(
          ((-1 + mu2) *(-1 + sigma2))/sigma2,1))))/(sigma2^3* beta(
    mu2 *(-1 + 1/sigma2), ((-1 + mu2)* (-1 + sigma2))/sigma2))

 
  
  
 
if(naive == FALSE){   


#set.seed(1)
#
#mu2 <- runif(2)
#sigma2 <- runif(2)
#y2 <- runif(2)
#
#mu2    <- mu2[2]
#sigma2 <- sigma2[2] 
#y2     <- y2[2]
#
#
#
#  func0 <- function(mu2){  pbeta(y2, shape1 = mu2 * (1 - sigma2)/sigma2, shape2 = (1-mu2)*(1 - sigma2)/sigma2) }
#  grad(func0 , mu2)
#  numDeriv::hessian(func0 , mu2)
#  #numericNHessian(func0, mu2)
#
#
#  func0 <- function(sigma2){  pbeta(y2, shape1 = mu2 * (1 - sigma2)/sigma2, shape2 = (1-mu2)*(1 - sigma2)/sigma2) }
#  grad(func0 , sigma2)
#  numDeriv::hessian(func0 , sigma2)
#  #numericNHessian(func0, sigma2)


p2  <-  pbeta(y2, shape1 = mu2 * (1 - sigma2)/sigma2, shape2 = (1-mu2)*(1 - sigma2)/sigma2)

funcD <- function(para) pbeta(y2, shape1 = para * (1 - sigma2)/sigma2, shape2 = (1-para)*(1 - sigma2)/sigma2)
 
nde <- numgh(funcD, mu2) 
 
derp2.dermu2   <- nde$fi
der2p2.dermu22 <- nde$se

funcD <- function(para) pbeta(y2, shape1 = mu2 * (1 - para)/para, shape2 = (1-mu2)*(1 - para)/para )
 

nde <- numgh(funcD, sigma2) 
 
derp2.dersigma2   <- nde$fi
der2p2.dersigma22 <- nde$se

funcD1 <- function(pms1, pms2) pbeta(y2, shape1 = pms1 * (1 - pms2)/pms2, shape2 = (1-pms1)*(1 - pms2)/pms2)

der2p2.derdermu2sigma2 <- numch(funcD1, mu2, sigma2)


#funcD1 <- function(pms) pbeta(y2, shape1 = pms[1] * (1 - pms[2])/pms[2], shape2 = (1-pms[1])*(1 - pms[2])/pms[2])
#numericNHessian(funcD1, c(mu2,sigma2))



                                             }


}




















####

if(margin2 == "iG"){


#
#sigma2    <- ifelse(sigma2 < 0.0001234098, 0.0001234098, sigma2)
#sigma2.st <- log(sigma2) 
#



mu2 <- dermu2.dereta2 <- der2mu2.dereta2eta2 <- exp(eta2) 
dersigma2.dersigma2.st <- exp(sigma2.st)                

                
pdf2          <- exp(-0.5 * log(2 * pi) - log(sqrt(sigma2)) - (3/2) * log(y2) - 
                   ((y2 - mu2)^2)/(2 * sigma2 * (mu2^2) * y2))
  
derpdf2.dermu2 <- (1/(mu2^2 * sigma2 * y2) + 4 * (mu2 * sigma2 * y2 * 
    (y2 - mu2)/(2 * (mu2^2 * sigma2 * y2))^2)) * exp(-((y2 - 
    mu2)^2/(2 * (mu2^2 * sigma2 * y2)) + (0.5 * (0.693147180559945 + 
    log(pi)) + 0.5 * log(sigma2) + 1.5 * log(y2)))) * (y2 - mu2)
         
           
   derpdf2.sigma2 <-  -((0.5/sigma2 - 2 * (mu2^2 * y2 * (y2 - mu2)^2/(2 * (mu2^2 * 
    sigma2 * y2))^2)) * exp(-((y2 - mu2)^2/(2 * (mu2^2 * sigma2 * 
    y2)) + (0.5 * (0.693147180559945 + log(pi)) + 0.5 * log(sigma2) + 
    1.5 * log(y2)))))


der2pdf2.dermu2 <-  (((1/(mu2^2 * sigma2 * y2) + 4 * (mu2 * sigma2 * y2 * 
    (y2 - mu2)/(2 * (mu2^2 * sigma2 * y2))^2))^2 * (y2 - mu2) + 
    sigma2 * y2 * (4 * (y2 - mu2) - mu2 * (16 + 64 * (mu2^3 * 
        sigma2^2 * y2^2 * (y2 - mu2)/(2 * (mu2^2 * sigma2 * y2))^2)))/(2 * 
        (mu2^2 * sigma2 * y2))^2) * (y2 - mu2) - 1/(mu2^2 * sigma2 * 
    y2)) * exp(-((y2 - mu2)^2/(2 * (mu2^2 * sigma2 * y2)) + (0.5 * 
    (0.693147180559945 + log(pi)) + 0.5 * log(sigma2) + 1.5 * 
    log(y2))))

   
der2pdf2.dersigma22 <- -((16 * (mu2^6 * sigma2 * y2^3 * (y2 - mu2)^2/(2 * (mu2^2 * sigma2 * 
    y2))^4) - ((0.5/sigma2 - 2 * (mu2^2 * y2 * (y2 - mu2)^2/(2 * 
    (mu2^2 * sigma2 * y2))^2))^2 + 0.5/(sigma2^1.5 * sqrt(sigma2)))) * 
    exp(-((y2 - mu2)^2/(2 * (mu2^2 * sigma2 * y2)) + (0.5 * (0.693147180559945 + 
        log(pi)) + 0.5 * log(sigma2) + 1.5 * log(y2)))))
 

der2pdf2.mu2dersigma2 <-  exp(-((y2 - mu2)^2/(2 * (mu2^2 * sigma2 * y2)) + (0.5 * 
    (0.693147180559945 + log(pi)) + 0.5 * log(sigma2) + 1.5 * 
    log(y2)))) * (mu2 * y2 * ((4 - 32 * (mu2^4 * sigma2^2 * y2^2/(2 * 
    (mu2^2 * sigma2 * y2))^2)) * (y2 - mu2) - 4 * mu2)/(2 * (mu2^2 * 
    sigma2 * y2))^2 - (0.5/sigma2 - 2 * (mu2^2 * y2 * (y2 - mu2)^2/(2 * 
    (mu2^2 * sigma2 * y2))^2)) * (1/(mu2^2 * sigma2 * y2) + 4 * 
    (mu2 * sigma2 * y2 * (y2 - mu2)/(2 * (mu2^2 * sigma2 * y2))^2))) * 
    (y2 - mu2)


                  
if(naive == FALSE){                   
          

p2          <-  pnorm(((y2/mu2) - 1)/(sqrt(sigma2) * sqrt(y2))) + 
                    exp(2/(mu2*sigma2))* pnorm(-((y2/mu2) + 1)/(sqrt(sigma2) * sqrt(y2)))
                                   
derp2.dermu2 <- exp(2/(mu2 * sigma2)) * (y2 * dnorm(-((1 + y2/mu2)/(sqrt(sigma2) * 
    sqrt(y2))))/(mu2^2 * sqrt(sigma2) * sqrt(y2)) - 2 * (sigma2 * 
    pnorm(-((1 + y2/mu2)/(sqrt(sigma2) * sqrt(y2))))/(mu2 * sigma2)^2)) - 
    y2 * dnorm((y2/mu2 - 1)/(sqrt(sigma2) * sqrt(y2)))/(mu2^2 * 
        sqrt(sigma2) * sqrt(y2))

      
derp2.dersigma2 <- (0.5 * ((1 + y2/mu2) * dnorm(-((1 + y2/mu2)/(sqrt(sigma2) * sqrt(y2)))) * 
    sqrt(y2)/(sqrt(sigma2) * (sqrt(sigma2) * sqrt(y2))^2)) - 
    2 * (mu2 * pnorm(-((1 + y2/mu2)/(sqrt(sigma2) * sqrt(y2))))/(mu2 * 
        sigma2)^2)) * exp(2/(mu2 * sigma2)) - 0.5 * (dnorm((y2/mu2 - 
    1)/(sqrt(sigma2) * sqrt(y2))) * sqrt(y2) * (y2/mu2 - 1)/(sqrt(sigma2) * 
    (sqrt(sigma2) * sqrt(y2))^2))          
    

der2p2.dermu22 <-   -(exp(2/(mu2 * sigma2)) * (sigma2 * (2 * (y2 * dnorm(-((1 + 
    y2/mu2)/(sqrt(sigma2) * sqrt(y2))))/(mu2^2 * sqrt(sigma2) * 
    sqrt(y2))) - sigma2 * (4 + 4 * (mu2 * sigma2)) * pnorm(-((1 + 
    y2/mu2)/(sqrt(sigma2) * sqrt(y2))))/(mu2 * sigma2)^2)/(mu2 * 
    sigma2)^2 + y2 * ((2 - (1 + y2/mu2)/(mu2 * sigma2))/mu2 + 
    2 * (sigma2/(mu2 * sigma2)^2)) * dnorm(-((1 + y2/mu2)/(sqrt(sigma2) * 
    sqrt(y2))))/(mu2^2 * sqrt(sigma2) * sqrt(y2))) + y2 * ((y2/mu2 - 
    1)/(mu2 * sigma2) - 2) * dnorm((y2/mu2 - 1)/(sqrt(sigma2) * 
    sqrt(y2)))/(mu2^3 * sqrt(sigma2) * sqrt(y2)))

  
der2p2.dersigma22 <-  (((0.25 * ((1 + y2/mu2)^2/(sqrt(sigma2) * sqrt(y2))^2) - 0.25)/sigma2 - 
    (0.5 * (y2/(sqrt(sigma2) * sqrt(y2))^2) + mu2/(mu2 * sigma2)^2)) * 
    (1 + y2/mu2) * dnorm(-((1 + y2/mu2)/(sqrt(sigma2) * sqrt(y2)))) * 
    sqrt(y2)/(sqrt(sigma2) * (sqrt(sigma2) * sqrt(y2))^2) - mu2 * 
    ((1 + y2/mu2) * dnorm(-((1 + y2/mu2)/(sqrt(sigma2) * sqrt(y2)))) * 
        sqrt(y2)/(sqrt(sigma2) * (sqrt(sigma2) * sqrt(y2))^2) - 
        mu2 * (4 + 4 * (mu2 * sigma2)) * pnorm(-((1 + y2/mu2)/(sqrt(sigma2) * 
            sqrt(y2))))/(mu2 * sigma2)^2)/(mu2 * sigma2)^2) * 
    exp(2/(mu2 * sigma2)) - ((0.25 * ((y2/mu2 - 1)^2/(sqrt(sigma2) * 
    sqrt(y2))^2) - 0.25)/sigma2 - 0.5 * (y2/(sqrt(sigma2) * sqrt(y2))^2)) * 
    dnorm((y2/mu2 - 1)/(sqrt(sigma2) * sqrt(y2))) * sqrt(y2) * 
    (y2/mu2 - 1)/(sqrt(sigma2) * (sqrt(sigma2) * sqrt(y2))^2)



der2p2.derdermu2sigma2 <-  -((((2 - mu2 * sigma2 * (4 + 4 * (mu2 * sigma2))/(mu2 * 
    sigma2)^2) * pnorm(-((1 + y2/mu2)/(sqrt(sigma2) * sqrt(y2)))) + 
    sigma2 * (1 + y2/mu2) * dnorm(-((1 + y2/mu2)/(sqrt(sigma2) * 
        sqrt(y2)))) * sqrt(y2)/(sqrt(sigma2) * (sqrt(sigma2) * 
        sqrt(y2))^2))/(mu2 * sigma2)^2 + y2 * ((0.5 * sqrt(y2) - 
    0.5 * ((1 + y2/mu2)^2/(sigma2 * sqrt(y2))))/(mu2 * (sqrt(sigma2) * 
    sqrt(y2))^2) + 2/((mu2 * sigma2)^2 * sqrt(y2))) * dnorm(-((1 + 
    y2/mu2)/(sqrt(sigma2) * sqrt(y2))))/(mu2 * sqrt(sigma2))) * 
    exp(2/(mu2 * sigma2)) + y2 * (0.5 * ((y2/mu2 - 1)^2/(sigma2 * 
    sqrt(y2))) - 0.5 * sqrt(y2)) * dnorm((y2/mu2 - 1)/(sqrt(sigma2) * 
    sqrt(y2)))/(mu2^2 * sqrt(sigma2) * (sqrt(sigma2) * sqrt(y2))^2))


                                                }





}


####

if(margin2 == "LO"){

mu2 <- eta2
dermu2.dereta2 <- 1
der2mu2.dereta2eta2 <- 0 
dersigma2.dersigma2.st <- exp(sigma2.st)  


  pdf2          <- dlogis(y2,eta2,sqrt(sigma2)) 
  
derpdf2.dermu2 <-     (exp((eta2 + y2)/sqrt(sigma2))* (-exp((eta2/sqrt(sigma2))) + exp(y2/sqrt(
                        sigma2))))/((exp(eta2/sqrt(sigma2)) + exp(y2/sqrt(sigma2)))^3 *sigma2)  
  
derpdf2.sigma2 <-  ((0.5 * eta2 - (0.5 * sqrt(sigma2) + 0.5 * y2)) * 
    exp(eta2/sqrt(sigma2)) + (0.5 * y2 - (0.5 * eta2 + 0.5 * 
    sqrt(sigma2))) * exp(y2/sqrt(sigma2))) * exp((eta2 + y2)/sqrt(sigma2))/(sigma2^2 * 
    (exp(eta2/sqrt(sigma2)) + exp(y2/sqrt(sigma2)))^3)
  
 
 der2pdf2.dermu2 <-  exp((eta2 + y2)/sqrt(sigma2)) * (exp(2 * (eta2/sqrt(sigma2))) + 
    exp(2 * (y2/sqrt(sigma2))) - 4 * exp((eta2 + y2)/sqrt(sigma2)))/(sigma2^1.5 * 
    (exp(eta2/sqrt(sigma2)) + exp(y2/sqrt(sigma2)))^4)
  
  
der2pdf2.dersigma22 <-  ((0.75 * sigma2 + eta2 * (0.25 * eta2 - (0.5 * y2 + 
    1.25 * sqrt(sigma2))) + y2 * (0.25 * y2 + 1.25 * sqrt(sigma2))) * 
    exp(2 * (eta2/sqrt(sigma2))) + (0.75 * sigma2 + eta2 * (0.25 * 
    eta2 + 1.25 * sqrt(sigma2) - 0.5 * y2) + y2 * (0.25 * y2 - 
    1.25 * sqrt(sigma2))) * exp(2 * (y2/sqrt(sigma2))) + (1.5 * 
    sigma2 + eta2 * (2 * y2 - eta2) - y2^2) * exp((eta2 + y2)/sqrt(sigma2))) * 
    exp((eta2 + y2)/sqrt(sigma2))/(sigma2^3.5 * (exp(eta2/sqrt(sigma2)) + 
    exp(y2/sqrt(sigma2)))^4)
  
  
  
  
der2pdf2.mu2dersigma2 <-  ((0.5 * y2 - (0.5 * eta2 + sqrt(sigma2))) * exp(2 * 
    (y2/sqrt(sigma2))) + (0.5 * y2 + sqrt(sigma2) - 0.5 * eta2) * 
    exp(2 * (eta2/sqrt(sigma2))) + (2 * eta2 - 2 * y2) * exp((eta2 + 
    y2)/sqrt(sigma2))) * exp((eta2 + y2)/sqrt(sigma2))/(sigma2^2.5 * 
    (exp(eta2/sqrt(sigma2)) + exp(y2/sqrt(sigma2)))^4)
  
  
  
if(naive == FALSE){   
  
    p2          <- plogis(y2,eta2,sqrt(sigma2)) 
             

derp2.dermu2    <- -(exp(-(y2 - eta2)/sqrt(sigma2)) * (1/sqrt(sigma2))/(1 + exp(-(y2 - 
                      eta2)/sqrt(sigma2)))^2)
                    
         
derp2.dersigma2 <-  -(exp(-(y2 - eta2)/sqrt(sigma2)) * ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)/(1 + 
                      exp(-(y2 - eta2)/sqrt(sigma2)))^2)                     
    


der2p2.dermu22 <-  (exp((2 * eta2 + y2)/sqrt(sigma2)) - exp((2 * y2 + 
    eta2)/sqrt(sigma2)))/(sigma2 * (exp(eta2/sqrt(sigma2)) + 
    exp(y2/sqrt(sigma2)))^3)


der2p2.dersigma22 <-    (0.5 * ((eta2 - y2) * exp(2 * ((eta2 - y2)/sqrt(sigma2)))/(1 + 
    exp((eta2 - y2)/sqrt(sigma2)))^3) - 0.25 * ((3 * sqrt(sigma2) + 
    eta2 - y2) * exp((eta2 + y2)/sqrt(sigma2))/(exp(eta2/sqrt(sigma2)) + 
    exp(y2/sqrt(sigma2)))^2)) * (eta2 - y2)/sigma2^3

der2p2.derdermu2sigma2 <-  ((0.5 * eta2 + 0.5 * sqrt(sigma2) - 0.5 * y2) * exp(y2/sqrt(sigma2)) + 
    (0.5 * sqrt(sigma2) + 0.5 * y2 - 0.5 * eta2) * exp(eta2/sqrt(sigma2))) * 
    exp((eta2 + y2)/sqrt(sigma2))/(sigma2^2 * (exp(eta2/sqrt(sigma2)) + 
    exp(y2/sqrt(sigma2)))^3)
            

                                         }



}


##########################


if(margin2 == "rGU"){

mu2 <- eta2
dermu2.dereta2 <- 1
der2mu2.dereta2eta2 <- 0 
dersigma2.dersigma2.st <- exp(sigma2.st)    


  pdf2          <- 1/sqrt(sigma2)*exp(-((y2-eta2)/sqrt(sigma2)+exp(-((y2-eta2)/sqrt(sigma2)))))
  
derpdf2.dermu2 <-   -(exp(-((y2 - eta2)/sqrt(sigma2) + exp(-((y2 - eta2)/sqrt(sigma2))))) * 
    (exp(-((y2 - eta2)/sqrt(sigma2))) - 1)/sigma2^1)
  
derpdf2.sigma2 <-  ((0.5 * eta2 - 0.5 * y2) * exp(eta2/sqrt(sigma2)) + 
    (0.5 * y2 - (0.5 * eta2 + 0.5 * sqrt(sigma2))) * exp(y2/sqrt(sigma2))) * 
    exp((eta2 - 2 * y2)/sqrt(sigma2) - exp((eta2 - y2)/sqrt(sigma2)))/sigma2^2


  
 der2pdf2.dermu2 <-   -(exp(-((y2 - eta2)/sqrt(sigma2) + exp(-((y2 - eta2)/sqrt(sigma2))))) * 
    (exp(-((y2 - eta2)/sqrt(sigma2))) - (exp(-((y2 - eta2)/sqrt(sigma2))) - 
        1)^2)/(sigma2 * sqrt(sigma2)))
  
der2pdf2.dersigma22 <-  ((0.25 * y2^2 + eta2 * (0.25 * eta2 - 0.5 * y2)) * 
    exp(2 * (eta2/sqrt(sigma2))) + (0.75 * sigma2 + eta2 * (0.25 * 
    eta2 + 1.25 * sqrt(sigma2) - 0.5 * y2) + y2 * (0.25 * y2 - 
    1.25 * sqrt(sigma2))) * exp(2 * (y2/sqrt(sigma2))) + (eta2 * 
    (1.5 * y2 - (0.75 * eta2 + 1.25 * sqrt(sigma2))) + y2 * (1.25 * 
    sqrt(sigma2) - 0.75 * y2)) * exp((eta2 + y2)/sqrt(sigma2))) * 
    exp((eta2 - 3 * y2)/sqrt(sigma2) - exp((eta2 - y2)/sqrt(sigma2)))/sigma2^3.5
  
  

der2pdf2.mu2dersigma2 <-  ((0.5 * y2 - (0.5 * eta2 + sqrt(sigma2))) * exp(2 * 
    (y2/sqrt(sigma2))) + (0.5 * y2 - 0.5 * eta2) * exp(2 * (eta2/sqrt(sigma2))) + 
    (1.5 * eta2 + sqrt(sigma2) - 1.5 * y2) * exp((eta2 + y2)/sqrt(sigma2))) * 
    exp((eta2 - 3 * y2)/sqrt(sigma2) - exp((eta2 - y2)/sqrt(sigma2)))/sigma2^2.5
  
  
 
  
if(naive == FALSE){   
  
  
  
    p2          <- exp(-(exp(-(y2-eta2)/sqrt(sigma2))))
                

derp2.dermu2    <- -(exp(-(exp(-(y2 - eta2)/sqrt(sigma2)))) * (exp(-(y2 - eta2)/sqrt(sigma2)) * 
                      (1/sqrt(sigma2))))
                      
                          
derp2.dersigma2 <-   -(exp(-(exp(-(y2 - eta2)/sqrt(sigma2)))) * (exp(-(y2 - eta2)/sqrt(sigma2)) * 
                     ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)))                      
    

der2p2.dermu22 <-  -((1 - exp(-((y2 - eta2)/sqrt(sigma2)))) * exp(-((y2 - 
    eta2)/sqrt(sigma2))) * exp(-exp(-((y2 - eta2)/sqrt(sigma2))))/sigma2)


der2p2.dersigma22 <-  (((0.25 * y2 - 0.25 * eta2)/sqrt(sigma2) - 0.75) * 
    exp(y2/sqrt(sigma2)) + 0.25 * ((eta2 - y2) * exp(eta2/sqrt(sigma2))/sqrt(sigma2))) * 
    (eta2 - y2) * exp((eta2 - 2 * y2)/sqrt(sigma2) - exp((eta2 - 
    y2)/sqrt(sigma2)))/sigma2^2.5

der2p2.derdermu2sigma2 <- -(exp(-(exp(-(y2 - eta2)/sqrt(sigma2)))) * (exp(-(y2 - eta2)/sqrt(sigma2)) * 
                            ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2) * (1/sqrt(sigma2)) - 
                             exp(-(y2 - eta2)/sqrt(sigma2)) * (0.5 * sigma2^-0.5/sqrt(sigma2)^2)) - 
                             exp(-(exp(-(y2 - eta2)/sqrt(sigma2)))) * (exp(-(y2 - eta2)/sqrt(sigma2)) * 
                            ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)) * 
                            (exp(-(y2 - eta2)/sqrt(sigma2)) * (1/sqrt(sigma2))))
            

                                   }

}


######################################


if(margin2 == "GU"){

mu2 <- eta2
dermu2.dereta2 <- 1
der2mu2.dereta2eta2 <- 0 
dersigma2.dersigma2.st <- exp(sigma2.st)  

  pdf2          <- exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * (1/sqrt(sigma2)))
                   
 derpdf2.dermu2 <-    exp(-exp((y2 - eta2)/sqrt(sigma2))) * exp((y2 - eta2)/sqrt(sigma2)) * 
    (exp((y2 - eta2)/sqrt(sigma2)) - 1)/sigma2            
                   
derpdf2.sigma2 <-  ((0.5 * eta2 - (0.5 * sqrt(sigma2) + 0.5 * y2)) * 
    exp(eta2/sqrt(sigma2)) + (0.5 * y2 - 0.5 * eta2) * exp(y2/sqrt(sigma2))) * 
    exp((y2 - 2 * eta2)/sqrt(sigma2) - exp((y2 - eta2)/sqrt(sigma2)))/sigma2^2

          
der2pdf2.dermu2 <-  exp((y2 - 3 * eta2)/sqrt(sigma2) - exp((y2 - eta2)/sqrt(sigma2))) * 
    (exp(2 * (eta2/sqrt(sigma2))) + exp(2 * (y2/sqrt(sigma2))) - 
        3 * exp((eta2 + y2)/sqrt(sigma2)))/sigma2^1.5             
                   
der2pdf2.dersigma22 <-  ((0.25 * y2^2 + eta2 * (0.25 * eta2 - 0.5 * y2)) * 
    exp(2 * (y2/sqrt(sigma2))) + (0.75 * sigma2 + eta2 * (0.25 * 
    eta2 - (0.5 * y2 + 1.25 * sqrt(sigma2))) + y2 * (0.25 * y2 + 
    1.25 * sqrt(sigma2))) * exp(2 * (eta2/sqrt(sigma2))) + (eta2 * 
    (1.25 * sqrt(sigma2) + 1.5 * y2 - 0.75 * eta2) - y2 * (0.75 * 
    y2 + 1.25 * sqrt(sigma2))) * exp((eta2 + y2)/sqrt(sigma2))) * 
    exp((y2 - 3 * eta2)/sqrt(sigma2) - exp((y2 - eta2)/sqrt(sigma2)))/sigma2^3.5               
                   
                   
                   
der2pdf2.mu2dersigma2 <- ((0.5 * y2 - 0.5 * eta2) * exp(2 * (y2/sqrt(sigma2))) + 
    (0.5 * y2 + sqrt(sigma2) - 0.5 * eta2) * exp(2 * (eta2/sqrt(sigma2))) + 
    (1.5 * eta2 - (1.5 * y2 + sqrt(sigma2))) * exp((eta2 + y2)/sqrt(sigma2))) * 
    exp((y2 - 3 * eta2)/sqrt(sigma2) - exp((y2 - eta2)/sqrt(sigma2)))/sigma2^2.5        
                   
            
            
            
if(naive == FALSE){                    
                   
    p2          <- 1 - exp(-exp((y2 - eta2)/sqrt(sigma2)))
    
 

derp2.dermu2    <- -(exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * 
                     (1/sqrt(sigma2))))
                  


                          
derp2.dersigma2 <-   -(exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * 
                      ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)))                       
    


der2p2.dermu22 <-   -(exp(-exp((y2 - eta2)/sqrt(sigma2))) * exp((y2 - 
    eta2)/sqrt(sigma2)) * (exp((y2 - eta2)/sqrt(sigma2)) - 1)/sigma2)


 der2p2.dersigma22 <-  (0.25 * ((eta2 - (3 * sqrt(sigma2) + y2)) * (eta2 - 
    y2) * exp(eta2/sqrt(sigma2))) - 0.25 * (exp(y2/sqrt(sigma2)) * 
    (y2 - eta2)^2)) * exp((y2 - 2 * eta2)/sqrt(sigma2) - exp((y2 - 
    eta2)/sqrt(sigma2)))/sigma2^3

der2p2.derdermu2sigma2 <-  -(exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * 
                            ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)) * (exp((y2 - 
                              eta2)/sqrt(sigma2)) * (1/sqrt(sigma2))) - exp(-exp((y2 - 
                              eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * (0.5 * 
                              sigma2^-0.5/sqrt(sigma2)^2) + exp((y2 - eta2)/sqrt(sigma2)) * 
                            ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2) * (1/sqrt(sigma2))))
            

                                           }

}



if(margin2 %in% c("GA","GAi")){

#
sigma2    <- ifelse(sigma2 < 0.006, 0.006, sigma2) # related to gamma function
sigma2.st <- log(sigma2) 
#


if(margin2 == "GAi"){eta2 <- ifelse(eta2 < 0.0000001, 0.0000001, eta2); mu2 <- eta2; dermu2.dereta2 <- 1; der2mu2.dereta2eta2 <- 0}
if(margin2 == "GA") { mu2 <- dermu2.dereta2 <- der2mu2.dereta2eta2 <- exp(eta2) }


dersigma2.dersigma2.st <- exp(sigma2.st)   

pdf2 <- dgamma(y2, shape = 1/sigma2, scale = mu2 * sigma2)

#expr <- expression(1/(sigma2*mu2)^(1/sigma2)*( y2^(1/sigma2-1)*exp(-y2/(sigma2*mu2)) )  /gamma(1/sigma2))
#Simplify( D(expr, "sigma2") )
#Simplify( D(D(expr, "sigma2"),"sigma2") )


derpdf2.dermu2 <-  -(sigma2 * exp((log(y2) - (log(mu2) + log(sigma2) + 
    y2/mu2))/sigma2 - (lgamma(1/sigma2) + log(y2))) * (mu2 - 
    y2)/(mu2 * sigma2)^2)
    
 
derpdf2.sigma2 <-  ((mu2 * y2^(1/sigma2)/(mu2 * sigma2)^2 - y2^(1/sigma2 - 1) * 
    log(y2)/sigma2^2)/(mu2 * sigma2)^(1/sigma2) + (y2^(1/sigma2 - 
    1) * digamma(1/sigma2)/(sigma2 * (mu2 * sigma2)^(1/sigma2)) - 
    y2^(1/sigma2 - 1) * (mu2 * (mu2 * sigma2)^(1/sigma2 - 1) - 
        (log(mu2) + log(sigma2)) * (mu2 * sigma2)^(1/sigma2)/sigma2)/(mu2 * 
        sigma2)^(2/sigma2))/sigma2) * exp(-(y2/(mu2 * sigma2)))/gamma(1/sigma2) 
 
 der2pdf2.dermu2 <-  sigma2 * exp((log(y2) - (log(mu2) + log(sigma2) + 
    y2/mu2))/sigma2 - (lgamma(1/sigma2) + log(y2))) * (sigma2 * 
    ((mu2 - y2)^2 + mu2 * sigma2 * (2 * mu2 - 2 * y2)) - y2^2/(y2/(mu2 * 
    sigma2))^2)/(mu2 * sigma2)^4
            
      
der2pdf2.dersigma22 <-  (((((2 * (y2^(1/sigma2 - 1) * digamma(1/sigma2)/(sigma2 * (mu2 * 
    sigma2)^(1/sigma2))) - 2 * (y2^(1/sigma2 - 1) * (mu2 * (mu2 * 
    sigma2)^(1/sigma2 - 1) - (log(mu2) + log(sigma2)) * (mu2 * 
    sigma2)^(1/sigma2)/sigma2)/(mu2 * sigma2)^(2/sigma2)))/sigma2 + 
    2 * ((mu2 * y2^(1/sigma2)/(mu2 * sigma2)^2 - y2^(1/sigma2 - 
        1) * log(y2)/sigma2^2)/(mu2 * sigma2)^(1/sigma2))) * 
    digamma(1/sigma2) - y2^(1/sigma2 - 1) * ((digamma(1/sigma2)^2 + 
    trigamma(1/sigma2))/sigma2 + 2 * digamma(1/sigma2))/(sigma2 * 
    (mu2 * sigma2)^(1/sigma2)))/sigma2 - (2 * ((mu2 * (mu2 * 
    sigma2)^(1/sigma2 - 1) - (log(mu2) + log(sigma2)) * (mu2 * 
    sigma2)^(1/sigma2)/sigma2) * (mu2 * y2^(1/sigma2)/(mu2 * 
    sigma2)^2 - y2^(1/sigma2 - 1) * log(y2)/sigma2^2)/(mu2 * 
    sigma2)^(2/sigma2)) + y2^(1/sigma2 - 1) * ((mu2 * (mu2 * 
    (1/sigma2 - 1) * (mu2 * sigma2)^(1/sigma2 - 2) - ((log(mu2) + 
    log(sigma2)) * (mu2 * sigma2)^(1/sigma2 - 1)/sigma2 + (mu2 * 
    sigma2)^(1/sigma2 - 1))/sigma2) - ((1 - 2 * (log(mu2) + log(sigma2))) * 
    (mu2 * sigma2)^(1/sigma2) + (log(mu2) + log(sigma2)) * (mu2 * 
    (mu2 * sigma2)^(1/sigma2 - 1) - (log(mu2) + log(sigma2)) * 
    (mu2 * sigma2)^(1/sigma2)/sigma2))/sigma2^2)/(mu2 * sigma2)^(2/sigma2) - 
    2 * ((mu2 * (mu2 * sigma2)^(1/sigma2 - 1) - (log(mu2) + log(sigma2)) * 
        (mu2 * sigma2)^(1/sigma2)/sigma2)^2/(sigma2 * (mu2 * 
        sigma2)^(3/sigma2))))))/sigma2 + (mu2 * (mu2 * y2^(1/sigma2) * 
    (y2 - 2 * (mu2 * sigma2))/(mu2 * sigma2)^2 - y2^(1/sigma2) * 
    log(y2)/sigma2^2)/(mu2 * sigma2)^2 - log(y2) * (mu2 * y2^(1/sigma2)/(mu2 * 
    sigma2)^2 - (2 * y2^(1/sigma2 - 1) + y2^(1/sigma2 - 1) * 
    log(y2)/sigma2)/sigma2)/sigma2^2)/(mu2 * sigma2)^(1/sigma2)) * 
    exp(-(y2/(mu2 * sigma2)))/gamma(1/sigma2)
                    
 der2pdf2.mu2dersigma2 <- -(exp((1/sigma2) * log(y2/(mu2 * sigma2)) - y2/(mu2 * sigma2) - 
    log(y2) - lgamma(1/sigma2)) * ((1/sigma2) * ((y2/(mu2 * sigma2)^2 - 
    y2 * sigma2 * (2 * (mu2 * (mu2 * sigma2)))/((mu2 * sigma2)^2)^2)/(y2/(mu2 * 
    sigma2)) + y2 * sigma2/(mu2 * sigma2)^2 * (y2 * mu2/(mu2 * 
    sigma2)^2)/(y2/(mu2 * sigma2))^2) - 1/sigma2^2 * (y2 * sigma2/(mu2 * 
    sigma2)^2/(y2/(mu2 * sigma2))) - (y2/(mu2 * sigma2)^2 - y2 * 
    sigma2 * (2 * (mu2 * (mu2 * sigma2)))/((mu2 * sigma2)^2)^2)) - 
    exp((1/sigma2) * log(y2/(mu2 * sigma2)) - y2/(mu2 * sigma2) - 
        log(y2) - lgamma(1/sigma2)) * ((1/sigma2) * (y2 * mu2/(mu2 * 
        sigma2)^2/(y2/(mu2 * sigma2))) + 1/sigma2^2 * log(y2/(mu2 * 
        sigma2)) - y2 * mu2/(mu2 * sigma2)^2 - 1/sigma2^2 * digamma(1/sigma2)) * 
        ((1/sigma2) * (y2 * sigma2/(mu2 * sigma2)^2/(y2/(mu2 * 
            sigma2))) - y2 * sigma2/(mu2 * sigma2)^2))

    
if(naive == FALSE){      
              
    p2  <-  pgamma(y2, shape = 1/sigma2, scale = mu2 * sigma2)
                                  

   derp2.dermu2 <- -((exp(-(y2/(mu2 *sigma2)))* y2 *(y2/(mu2* sigma2))^(-1 + 1/sigma2))/(
                  mu2^2 *sigma2*gamma(1/sigma2)))                                           #done in Mathematica
   
  
der2p2.dermu22 <- y2 * (2 * (y2/(mu2 * sigma2))^(1/sigma2 - 1) + y2 * 
    ((1/sigma2 - 1) * (y2/(mu2 * sigma2))^(1/sigma2 - 2) - (y2/(mu2 * 
        sigma2))^(1/sigma2 - 1))/(mu2 * sigma2)) * exp(-(y2/(mu2 * 
    sigma2)))/(mu2^3 * sigma2 * gamma(1/sigma2))
   
      
der2p2.derdermu2sigma2 <- y2 * (((1 + log(y2) - (log(mu2) + log(sigma2)))/sigma2 - 
    1) * (y2/(mu2 * sigma2))^(1/sigma2 - 1) + (y2/(mu2 * sigma2))^(1/sigma2 - 
    1) - (psigamma(1/sigma2, 0) * (y2/(mu2 * sigma2))^(1/sigma2 - 
    1) + y2 * (y2/(mu2 * sigma2))^(1/sigma2 - 1)/mu2)/sigma2) * 
    exp(-(y2/(mu2 * sigma2)))/(mu2^2 * sigma2^2 * gamma(1/sigma2))
   

funcD <- function(para) pgamma(y2, shape = 1/para, scale = mu2 * para)
 
nde <- numgh(funcD, sigma2) 
 
derp2.dersigma2   <- nde$fi
der2p2.dersigma22 <- nde$se

                              }



}


##########################################################




if(margin2 %in% c(cont2par,cont3par)){ #### the stuff below does not apply to N case and there is a reason ####







derpdf2.dereta2              <- derpdf2.dermu2*dermu2.dereta2       
der2pdf2.dereta2             <- der2pdf2.dermu2* dermu2.dereta2^2 + derpdf2.dermu2*der2mu2.dereta2eta2        
der2pdf2.dereta2dersigma2    <- der2pdf2.mu2dersigma2* dermu2.dereta2
  
derpdf2.dersigma2.st         <- derpdf2.sigma2 * dersigma2.dersigma2.st   
der2pdf2.dersigma2.st2       <- der2pdf2.dersigma22 * dersigma2.dersigma2.st^2 + derpdf2.sigma2  * dersigma2.dersigma2.st     
der2pdf2.dereta2dersigma2.st <- der2pdf2.dereta2dersigma2 *  dersigma2.dersigma2.st




if(margin2 %in% cont3par){


der2pdf2.dereta2dernu        <- der2pdf2.mu2dernu* dermu2.dereta2

der2pdf2.dereta2dernu.st     <- der2pdf2.dereta2dernu * dernu.dernu.st
der2pdf2.sigma2.st2dernu.st  <- der2pdf2.dersigma2dernu * dersigma2.dersigma2.st * dernu.dernu.st                              
derpdf2.dernu.st             <- derpdf2.nu * dernu.dernu.st 
der2pdf2.dernu.st2           <- der2pdf2.dernu2 * dernu.dernu.st^2 +  derpdf2.nu  * dernu.dernu.st 

                         }
                         
                         
                         
                         
                        

if(naive == FALSE){  



derp2.dereta2                <- derp2.dermu2*dermu2.dereta2
der2p2.dereta2eta2           <- der2p2.dermu22*dermu2.dereta2^2+derp2.dermu2*der2mu2.dereta2eta2      
der2p2.dereta2dersigma2      <- der2p2.derdermu2sigma2* dermu2.dereta2    

derp2.dersigma.st            <- derp2.dersigma2 *  dersigma2.dersigma2.st 
der2p2.dersigma2.st2         <- der2p2.dersigma22 * dersigma2.dersigma2.st^2 + derp2.dersigma2 * dersigma2.dersigma2.st
der2p2.dereta2dersigma2.st   <- der2p2.dereta2dersigma2 *  dersigma2.dersigma2.st  



if(margin2 %in% cont3par){

der2p2.dereta2dernu          <- der2p2.dermu2dernu* dermu2.dereta2 
                             
derp2.nu.st                  <- derp2.dernu *  dernu.dernu.st 
der2p2.dernu.st2             <- der2p2.dernu2 * dernu.dernu.st^2 + derp2.dernu * dernu.dernu.st
der2p2.dereta2dernu.st       <- der2p2.dereta2dernu * dernu.dernu.st 
der2p2.dersigma2.stdernu.st  <- der2p2.dersigma2dernu * dersigma2.dersigma2.st * dernu.dernu.st

                         }


                 }
                 
                 
                 

}###############



epsilon <- 0.0000001 
max.p   <- 0.9999999

  pdf2 <- ifelse(pdf2 < epsilon, epsilon, pdf2 )

  p2   <- ifelse(p2 > max.p, max.p, p2) 
  p2   <- ifelse(p2 < epsilon, epsilon, p2) 



ifef <- function(dv){

epsilon <- 0.0000001 
dv <- ifelse(is.na(dv), epsilon, dv ) 
dv <- ifelse(dv == Inf ,  8.218407e+20, dv )
dv <- ifelse(dv == -Inf ,  -8.218407e+20, dv )
dv

}


# for safety

pdf2                         = ifef(pdf2)
p2                           = ifef(p2) 
derpdf2.dereta2              = ifef(derpdf2.dereta2) 
derpdf2.dersigma2.st         = ifef(derpdf2.dersigma2.st) 
derp2.dersigma.st            = ifef(derp2.dersigma.st)
derp2.dereta2                = ifef(derp2.dereta2)
der2p2.dereta2eta2           = ifef(der2p2.dereta2eta2) 
der2pdf2.dereta2             = ifef(der2pdf2.dereta2)
der2p2.dersigma2.st2         = ifef(der2p2.dersigma2.st2) 
der2pdf2.dersigma2.st2       = ifef(der2pdf2.dersigma2.st2)
der2p2.dereta2dersigma2.st   = ifef(der2p2.dereta2dersigma2.st)  
der2pdf2.dereta2dersigma2.st = ifef(der2pdf2.dereta2dersigma2.st)
der2pdf2.dereta2dernu.st     = ifef(der2pdf2.dereta2dernu.st)   
der2pdf2.sigma2.st2dernu.st  = ifef(der2pdf2.sigma2.st2dernu.st)
derpdf2.dernu.st             = ifef(derpdf2.dernu.st)           
der2pdf2.dernu.st2           = ifef(der2pdf2.dernu.st2)         
derp2.nu.st                  = ifef(derp2.nu.st)                
der2p2.dernu.st2             = ifef(der2p2.dernu.st2)           
der2p2.dereta2dernu.st       = ifef(der2p2.dereta2dernu.st)     
der2p2.dersigma2.stdernu.st  = ifef(der2p2.dersigma2.stdernu.st )




list(pdf2                         = pdf2,
     p2                           = p2, 
     derpdf2.dereta2              = derpdf2.dereta2, 
     derpdf2.dersigma2.st         = derpdf2.dersigma2.st, 
     derp2.dersigma.st            = derp2.dersigma.st,
     derp2.dereta2                = derp2.dereta2,
     der2p2.dereta2eta2           = der2p2.dereta2eta2, 
     der2pdf2.dereta2             = der2pdf2.dereta2,
     der2p2.dersigma2.st2         = der2p2.dersigma2.st2, 
     der2pdf2.dersigma2.st2       = der2pdf2.dersigma2.st2,
     der2p2.dereta2dersigma2.st   = der2p2.dereta2dersigma2.st,            
     der2pdf2.dereta2dersigma2.st = der2pdf2.dereta2dersigma2.st,
     der2pdf2.dereta2dernu.st     = der2pdf2.dereta2dernu.st,   
     der2pdf2.sigma2.st2dernu.st  = der2pdf2.sigma2.st2dernu.st,
     derpdf2.dernu.st             = derpdf2.dernu.st,           
     der2pdf2.dernu.st2           = der2pdf2.dernu.st2,         
     derp2.nu.st                  = derp2.nu.st,                
     der2p2.dernu.st2             = der2p2.dernu.st2,           
     der2p2.dereta2dernu.st       = der2p2.dereta2dernu.st,     
     der2p2.dersigma2.stdernu.st  = der2p2.dersigma2.stdernu.st )     


}




    