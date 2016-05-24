jc.probs <- function(x, y1, y2, newdata, type = "bivariate", cond = 0, intervals = FALSE, n.sim = 100, prob.lev = 0.05){

if(x$triv == TRUE ) stop("This function is currently not suitable for trivariate probit models.")


if(missing(y1)) stop("You must provide a value for y1.")
if(missing(y2)) stop("You must provide a value for y2.")

epsilon <- 0.0000001 
  
cont2par  <- x$VC$m2 
cont3par  <- x$VC$m3 

nu1 <- nu2 <- nu <- 1
CIp12 <- NULL
bin.link <- x$VC$bl  

if(!(type %in% c("bivariate","independence"))) stop("Error in parameter type value. It should be one of: bivariate, independence.")
if(!(cond %in% c(0,1,2))) stop("Error in parameter cond value. It should be one of: 0, 1, 2.")

if( type %in% c("independence") && x$VC$gamlssfit == FALSE) stop("You need to re-fit the model and set gamlssfit = TRUE to obtain probabilities under independence.")

if( x$VC$Cont == "NO" && !(x$margins[2] %in% bin.link) && !(y1 %in% c(0,1)) ) stop("The value for y1 must be either 0 or 1.")

if( x$VC$Cont == "NO" && x$margins[2] %in% bin.link){ if( !(y1 %in% c(0,1)) || !(y2 %in% c(0,1))   ) stop("The value for y1 and/or y2 must be either 0 or 1.") }



######################################################################################################
######################################################################################################

if(x$VC$Cont == "YES"){ ###


if(type == "bivariate"){ ##



if(!missing(newdata)){

nu1 <- nu2 <- NA

eta1 <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata)
eta2 <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata)

if( !is.null(x$X3) ){

sigma21 <- esp.tr(predict.SemiParBIVProbit(x, eq = 3, newdata = newdata), x$margins[1])$vrb
sigma22 <- esp.tr(predict.SemiParBIVProbit(x, eq = 4, newdata = newdata), x$margins[2])$vrb

if(x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){ eq.nu1 <- 5;  eq.nu2 <- 6;  eq.th <- 7}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){ eq.nu1 <- NA; eq.nu2 <- NA; eq.th <- 5}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 5;  eq.th <- 6}
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){ eq.nu1 <- 5;  eq.nu2 <- NA; eq.th <- 6}

if(x$margins[1] %in% cont3par) nu1 <- esp.tr(predict.SemiParBIVProbit(x, eq = eq.nu1, newdata = newdata), x$margins[1])$vrb
if(x$margins[2] %in% cont3par) nu2 <- esp.tr(predict.SemiParBIVProbit(x, eq = eq.nu2, newdata = newdata), x$margins[2])$vrb

theta <- teta.tr(x$VC, predict.SemiParBIVProbit(x, eq = eq.th, newdata = newdata))$teta

}

if( is.null(x$X3) ){

sigma21 <- x$sigma21
sigma22 <- x$sigma22

nu1 <- x$nu1 
nu2 <- x$nu2

theta <- x$theta 


}


}



if(missing(newdata)){

eta1 <- x$eta1
eta2 <- x$eta2

sigma21 <- x$sigma21
sigma22 <- x$sigma22

nu1 <- x$nu1 
nu2 <- x$nu2

theta <- x$theta 

}


p1  <- distrHsAT(y1, eta1, sigma21, nu1, x$margins[1])$p2 # I need to evaluate this at the y point of interest
p2  <- distrHsAT(y2, eta2, sigma22, nu2, x$margins[2])$p2 # so it can not be taken from the fitted model
p12 <- BiCDF(p1, p2, x$nC, theta)

if(cond == 1) p12 <- p12/p1
if(cond == 2) p12 <- p12/p2




if(intervals == TRUE){

bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)  
lf <- length(x$coefficients)


#############  
# etas
#############  

# try with 1 number

if(!missing(newdata)){ X1 <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix") 
                       X2 <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix") }
                       
if( missing(newdata)){ X1 <- x$X1 
                       X2 <- x$X2 }                       


eta1s <- eta.tr( X1%*%t(bs[,1:x$X1.d2])                     , x$VC$margins[1]) 
eta2s <- eta.tr( X2%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]) , x$VC$margins[2])

#############  
# thetas
#############  

if(  is.null(x$X3) ) epds <- bs[,lf]
  
  	if( !is.null(x$X3) ){ 
  	
  	
  	
  if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont2par){
  
  
       if(!missing(newdata)) X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X5 <- x$X5 
 
       epds <- X5%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)])
 
  }
 
  if((x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par) || (x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par) ){
  
       if(!missing(newdata)) X6 <- predict.SemiParBIVProbit(x, eq = 6, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X6 <- x$X6   
  
       epds <- X6%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2)])
  
  }
  
  if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par){
  
       if(!missing(newdata)) X7 <- predict.SemiParBIVProbit(x, eq = 7, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X7 <- x$X7    
  
       epds <- X7%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+x$X7.d2)])
  
  }
  
  
                            }
  	                        


est.RHOb <- teta.tr(x$VC, epds)$teta
   
#############  
# sigmas
#############  

      if( is.null(x$X3) ) { # could remove the conditions below and make it more general but need to define
                            # XX.d2 and change the multiplication operator. For now we will leave it like this   
      
if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont2par ){ ps1 <- lf - 2; ps2 <- lf - 1 }
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par ){ ps1 <- lf - 4; ps2 <- lf - 3 }
if((x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par) || (x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par) ){ ps1 <- lf - 3; ps2 <- lf - 2 }
      
                                sigma2.1.star <- bs[, ps1] 
                                sigma2.2.star <- bs[, ps2] 
                                
                                }
  
  
      if( !is.null(x$X3) ) {
      
      
       if(!missing(newdata)){ X3 <- predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")
                              X4 <- predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix") }
       if( missing(newdata)){ X3 <- x$X3; X4 <- x$X4 }        
      
      

       sigma2.1.star <- X3%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]) 
       sigma2.2.star <- X4%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)]) 

                                }  


    sigma21 <- esp.tr(sigma2.1.star, x$VC$margins[1])$vrb   
    sigma22 <- esp.tr(sigma2.2.star, x$VC$margins[2])$vrb   
    
#############  
# NUs
#############    
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3) )  {    pn1 <- lf - 2 
                            pn2 <- lf - 1
                            nu1.st <- bs[, pn1]    # t(as.matrix(bs[, pn1]))
                            nu2.st <- bs[, pn2]  } # t(as.matrix(bs[, pn2]))  } 
  
  if( !is.null(x$X3) ) {  
  
  
       if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")
                              X6 <- predict.SemiParBIVProbit(x, eq = 6, newdata = newdata, type = "lpmatrix") }
       if( missing(newdata)){ X5 <- x$X5; X6 <- x$X6 }   
 
       nu1.st <- X5%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)]) 
       nu2.st <- X6%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2)])
  
                       }   
   
nu1 <- esp.tr(nu1.st, x$VC$margins[1])$vrb   
nu2 <- esp.tr(nu2.st, x$VC$margins[2])$vrb   
  
} 


if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par ){  
  
  if( is.null(x$X3) )  {  pn2 <- lf - 1; nu2.st <- bs[, pn2]  } # nu2.st <- t(as.matrix(bs[, pn2]))
  
  
       if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")}
       if( missing(newdata)){ X5 <- x$X5}     
  
  if( !is.null(x$X3) ) nu2.st <- X5%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)]) 
  
  nu2 <- esp.tr(nu2.st, x$VC$margins[2])$vrb   

} 
  
  
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par ){  
    
  if( is.null(x$X3) )  {  pn1 <- lf - 1; nu1.st <- bs[, pn1]  } # nu1.st <- t(as.matrix(bs[, pn1])) 
  
       if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")}
       if( missing(newdata)){ X5 <- x$X5}    
  
  if( !is.null(x$X3) ) nu1.st <- X5%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)]) 
  
  nu1 <- esp.tr(nu1.st, x$VC$margins[1])$vrb  
   
}  


####################################################


if( is.null(x$X3) ){

est.RHOb <- matrix(rep(est.RHOb, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
sigma21  <- matrix(rep(sigma21, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
sigma22  <- matrix(rep(sigma22, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
nu1      <- matrix(rep(nu1, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
nu2      <- matrix(rep(nu2, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)

                   }

p1s  <- distrHsAT(y1, eta1s, sigma21, nu1, x$margins[1])$p2
p2s  <- distrHsAT(y2, eta2s, sigma22, nu2, x$margins[2])$p2


if(x$VC$BivD == "N") p12s <- matrix(BiCDF(p1s, p2s, x$nC, est.RHOb, test = FALSE), dim(p1s)[1], n.sim) else p12s <- BiCDF(p1s, p2s, x$nC, est.RHOb, test = FALSE)

if(cond == 1) p12s <- p12s/p1s
if(cond == 2) p12s <- p12s/p2s





} # interv






                        }## biv




if(type == "independence"){



if(!missing(newdata)){

nu1 <- nu2 <- NA

eta1 <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%x$gamlss1$fit$argument[1:x$X1.d2]
eta2 <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%x$gamlss2$fit$argument[1:x$X2.d2]

if( !is.null(x$X3) ){ 

sigma21 <- esp.tr(predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")%*%x$gamlss1$fit$argument[(x$X1.d2+1):(x$X1.d2+x$X3.d2)], x$margins[1])$vrb
sigma22 <- esp.tr(predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix")%*%x$gamlss2$fit$argument[(x$X2.d2+1):(x$X2.d2+x$X4.d2)], x$margins[2])$vrb

if(x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){ eq.nu1 <- 5;  eq.nu2 <- 6  ; indn1 <- (x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2); indn2 <- (x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X6.d2)}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){ eq.nu1 <- NA; eq.nu2 <- NA ; indn1 <- NA; indn2 <- NA}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 5  ; indn1 <- NA; indn2 <- (x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X5.d2)}
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){ eq.nu1 <- 5;  eq.nu2 <- NA ; indn1 <- (x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2); indn2 <- NA}

if(x$margins[1] %in% cont3par) nu1 <- esp.tr(predict.SemiParBIVProbit(x, eq = eq.nu1, newdata = newdata, type = "lpmatrix")%*%x$gamlss1$fit$argument[indn1], x$margins[1])$vrb
if(x$margins[2] %in% cont3par) nu2 <- esp.tr(predict.SemiParBIVProbit(x, eq = eq.nu2, newdata = newdata, type = "lpmatrix")%*%x$gamlss2$fit$argument[indn2], x$margins[2])$vrb

}

if( is.null(x$X3) ){ 

sigma21 <- x$gamlss1$fit$sigma2
sigma22 <- x$gamlss2$fit$sigma2

nu1 <- x$gamlss1$fit$nu
nu2 <- x$gamlss2$fit$nu

}

}



if(missing(newdata)){


eta1 <- x$gamlss1$fit$eta2
eta2 <- x$gamlss2$fit$eta2

sigma21 <- x$gamlss1$fit$sigma2
sigma22 <- x$gamlss2$fit$sigma2

nu1 <- x$gamlss1$fit$nu
nu2 <- x$gamlss2$fit$nu

}




p1  <- distrHsAT(y1, eta1, sigma21, nu1, x$margins[1])$p2
p2  <- distrHsAT(y2, eta2, sigma22, nu2, x$margins[2])$p2
p12 <- p1*p2

if(cond == 1) p12 <- p2
if(cond == 2) p12 <- p1



if(intervals == TRUE){

bs1 <- rMVN(n.sim, mean = x$gamlss1$fit$argument, sigma=x$gamlss1$magpp$Vb)
bs2 <- rMVN(n.sim, mean = x$gamlss2$fit$argument, sigma=x$gamlss2$magpp$Vb)


#############  
# etas
############# 

if(!missing(newdata)){ X1 <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix")
                       X2 <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix")}
if( missing(newdata)){ X1 <- x$X1; X2 <- x$X2}  

eta1s <- eta.tr( X1%*%t(bs1[,1:x$X1.d2]), x$VC$margins[1]) 
eta2s <- eta.tr( X2%*%t(bs2[,1:x$X2.d2]), x$VC$margins[2]) 

#############  
# sigmas
#############  

      if( is.null(x$X3) ) {
      
       sigma2.1.star <- bs1[, x$X1.d2 + 1] 
       sigma2.2.star <- bs2[, x$X2.d2 + 1] 
                            
                          }
  
  
      if( !is.null(x$X3) ) {

if(!missing(newdata)){ X3 <- predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")
                       X4 <- predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix")}
if( missing(newdata)){ X3 <- x$X3; X4 <- x$X4}  

       sigma2.1.star <- X3%*%t(bs1[,(x$X1.d2+1):(x$X1.d2+x$X3.d2)]) 
       sigma2.2.star <- X4%*%t(bs2[,(x$X2.d2+1):(x$X2.d2+x$X4.d2)]) 

                           }  

    sigma21 <- esp.tr(sigma2.1.star, x$VC$margins[1])$vrb   
    sigma22 <- esp.tr(sigma2.2.star, x$VC$margins[2])$vrb  

#############  
# NUs
#############    
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3) )  {     
       nu1.st <- bs1[, x$X1.d2 + 2]  # t(as.matrix(bs1[, x$X1.d2 + 2]))
       nu2.st <- bs2[, x$X2.d2 + 2]  # t(as.matrix(bs2[, x$X2.d2 + 2]))   
                        } 
  
  if( !is.null(x$X3) ) {  
  
  
if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")
                       X6 <- predict.SemiParBIVProbit(x, eq = 6, newdata = newdata, type = "lpmatrix")}
if( missing(newdata)){ X5 <- x$X5; X6 <- x$X6}    
  
       nu1.st <- X5%*%t(bs1[,(x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2)]) 
       nu2.st <- X6%*%t(bs2[,(x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X6.d2)])
                       }   
   
nu1 <- esp.tr(nu1.st, x$VC$margins[1])$vrb   
nu2 <- esp.tr(nu2.st, x$VC$margins[2])$vrb   
  
} 


if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par ){  
  
  if( is.null(x$X3) )  nu2.st <- bs2[, x$X2.d2 + 2] # t(as.matrix(bs2[, x$X2.d2 + 2]))  
  
  
  if( !is.null(x$X3) ){
  
  
    if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")}
    if( missing(newdata)){ X5 <- x$X5}  
  
    nu2.st <- X5%*%t(bs2[,(x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X5.d2)]) 
  
  
                      }
  
  
                       nu2    <- esp.tr(nu2.st, x$VC$margins[2])$vrb   
} 
  
  
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par ){  
    
  if( is.null(x$X3) )  nu1.st <- bs1[, x$X1.d2 + 2] # t(as.matrix(bs1[, x$X1.d2 + 2]))   


  if( !is.null(x$X3) ){
    
    if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")}
    if( missing(newdata)){ X5 <- x$X5}    
  
    nu1.st <- X5%*%t(bs1[,(x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2)])                  
  
                      }
  
  nu1    <- esp.tr(nu1.st, x$VC$margins[1])$vrb    
}  


if( is.null(x$X3) ){

sigma21  <- matrix(rep(sigma21, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
sigma22  <- matrix(rep(sigma22, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
nu1      <- matrix(rep(nu1, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
nu2      <- matrix(rep(nu2, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)

                   }


p1s  <- distrHsAT(y1, eta1s, sigma21, nu1, x$margins[1])$p2
p2s  <- distrHsAT(y2, eta2s, sigma22, nu2, x$margins[2])$p2
p12s <- p1s*p2s

if(cond == 1) p12s <- p2s
if(cond == 2) p12s <- p1s



} # intervals


} # independence





} ### cont yes



######################################################################################################
######################################################################################################












if( x$VC$Cont == "NO" && !(x$margins[2] %in% bin.link) ){ 


if(type == "bivariate"){


if(!missing(newdata)){

nu <- NA

p1   <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "response")
eta2 <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata)


if( !is.null(x$X3) ){

sigma2 <- esp.tr(predict.SemiParBIVProbit(x, eq = 3, newdata = newdata), x$margins[2])$vrb

if(x$margins[2] %in% cont2par) eq.th <- 4
if(x$margins[2] %in% cont3par){ eq.nu <- 4; eq.th <- 5}

if(x$margins[2] %in% cont3par) nu <- esp.tr(predict.SemiParBIVProbit(x, eq = eq.nu, newdata = newdata), x$margins[2])$vrb

theta <- teta.tr(x$VC, predict.SemiParBIVProbit(x, eq = eq.th, newdata = newdata))$teta

}


if( is.null(x$X3) ){

sigma2 <- x$sigma2
nu     <- x$nu 
theta  <- x$theta 

}


}



if(missing(newdata)){

p1   <- x$p1
eta2 <- x$eta2

sigma2 <- x$sigma2
nu     <- x$nu 
theta  <- x$theta 

}


if(y1 == 0 || y1 == 1){              ### joint prob of y1=0 and y2=x ###

p0  <- 1 - p1         
p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2])$p2
p12 <- BiCDF(p0, p2, x$nC, theta)

if(cond == 1) p12 <- p12/p0
if(cond == 2) p12 <- p12/p2

                      }



if(y1 == 1){                         ### joint prob of y1=1 and y2=x ###

p12  <- p2 - p12

if(cond == 1) p12 <- p12/p1
if(cond == 2) p12 <- p12/p2

           }
           





           
           
if(intervals == TRUE){

bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)  

#############  
# etas
############# 

if(!missing(newdata)){ X1  <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix") 
                       X2s <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix") }
                       
if( missing(newdata)){ X1 <- x$X1 
                       if(x$VC$ccss == "yes") X2s <- x$X2s else X2s <- x$X2  
                       }  


p1s   <- probm( X1%*%t(bs[,1:x$X1.d2]), x$VC$margins[1])$pr 
eta2s <- eta.tr( X2s%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]) , x$VC$margins[2])

#############  
# thetas
#############  

if( is.null(x$X3) ) epds <- bs[,length(x$coefficients)]
  
if( !is.null(x$X3) ){ 
  if(x$VC$margins[2] %in% cont2par){   
  
  
if(!missing(newdata)){ X4s <- predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}    
  
                epds <- X4s%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)])
  
                                   }
                                   
  if(x$VC$margins[2] %in% cont3par){
  
if(!missing(newdata)){ X5s <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X5s <- x$X5s else X5s <- x$X5}    
  
                epds <- X5s%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)])
                
                                   }
  
  
  
  }

est.RHOb <- teta.tr(x$VC, epds)$teta
   
#############  
# sigmas
#############  

      if( is.null(x$X3) )   sigma2.star <- bs[, x$X1.d2 + x$X2.d2 + 1] 
      if( !is.null(x$X3) ) {
      
if(!missing(newdata)){ X3s <- predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X3s <- x$X3s else X3s <- x$X3}        
                
            sigma2.star <- X3s%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]) 
  
                           }
  
sigma2 <- esp.tr(sigma2.star, x$VC$margins[2])$vrb   
    
#############  
# NUs
#############    
  
if( x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3)  ) nu.st <- bs[, x$X1.d2 + x$X2.d2 + 2] # t(as.matrix(bs[,  x$X1.d2 + x$X2.d2 + 2]))
  if( !is.null(x$X3) ){
  
if(!missing(newdata)){ X4s <- predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix")}                    
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}    
  
              nu.st <- X4s%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)]) 
  
                      }
  
 nu <- esp.tr(nu.st, x$VC$margins[2])$vrb   
  
} 


#################


if( is.null(x$X3) ){

est.RHOb <- matrix(rep(est.RHOb, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
sigma2   <- matrix(rep(sigma2, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
nu       <- matrix(rep(nu, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)

                   }


if(y1 == 0 || y1 == 1){             
p0s  <- 1 - p1s         
p2s  <- distrHsAT(y2, eta2s, sigma2, nu, x$margins[2])$p2

if(x$VC$BivD == "N") p12s <- matrix(BiCDF(p0s, p2s, x$nC, est.RHOb, test = FALSE), dim(p1s)[1], n.sim) else p12s <- BiCDF(p0s, p2s, x$nC, est.RHOb, test = FALSE)

if(cond == 1) p12s <- p12s/p0s
if(cond == 2) p12s <- p12s/p2s

                      }

if(y1 == 1){                         
  
p12s  <- p2s - p12s

if(cond == 1) p12s <- p12s/p1s
if(cond == 2) p12s <- p12s/p2s

           }
           
           

} # int




     
} # biv








if(type == "independence"){



if(!missing(newdata)){

nu     <- NA
p1     <- probm(predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1])$pr
eta2   <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%x$gamlss$fit$argument[1:x$X2.d2]


if( !is.null(x$X3) ){

sigma2 <- esp.tr(predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")%*%x$gamlss$fit$argument[(x$X2.d2+1):(x$X2.d2+x$X3.d2)], x$margins[2])$vrb
if(x$margins[2] %in% cont3par) nu <- esp.tr(predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix")%*%x$gamlss$fit$argument[(x$X2.d2 + x$X3.d2 + 1):(x$X2.d2 + x$X3.d2 + x$X4.d2)], x$margins[2])$vrb

}

if( is.null(x$X3) ){

sigma2 <- x$gamlss$fit$sigma2
nu     <- x$gamlss$fit$nu 


}

}



if(missing(newdata)){

p1   <- probm(predict.SemiParBIVProbit(x, eq = 1, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1])$pr
eta2 <- x$gamlss$fit$eta2

sigma2 <- x$gamlss$fit$sigma2
nu     <- x$gamlss$fit$nu 

}



if(y1 == 0 || y1 == 1){  

p0  <- 1 - p1                     
p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2])$p2
p12 <- p0*p2

if(cond == 1) p12 <- p2
if(cond == 2) p12 <- p0

                      }

if(y1 == 1){                                          

p12 <- p1*p2 

if(cond == 1) p12 <- p2
if(cond == 2) p12 <- p1

           }
           
           
           




      
if(intervals == TRUE){


bs1 <- rMVN(n.sim, mean = x$gam1$coefficients, sigma=x$gam1$Vp)
bs2 <- rMVN(n.sim, mean = x$gamlss$fit$argument, sigma=x$gamlss$magpp$Vb)


#############  
# etas
#############  

if(!missing(newdata)){ X1  <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix") 
                       X2s <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix") }
                       
if( missing(newdata)){ X1 <- x$X1 
                       if(x$VC$ccss == "yes") X2s <- x$X2s else X2s <- x$X2  
                       } 
  
p1s   <- probm(   X1%*%t(bs1[,1:x$X1.d2]) , x$VC$margins[1])$pr   
eta2s <- eta.tr( X2s%*%t(bs2[,1:x$X2.d2]) , x$VC$margins[2])

#############  
# sigmas
#############  

      if( is.null(x$X3) )   sigma2.star <- bs2[, x$X2.d2 + 1] 
      
      if( !is.null(x$X3) ){
      
if(!missing(newdata)){ X3s <- predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")}                     
if( missing(newdata)){ if(x$VC$ccss == "yes") X3s <- x$X3s else X3s <- x$X3}        
      
                      sigma2.star <- X3s%*%t(bs2[,(x$X2.d2+1):(x$X2.d2+x$X3.d2)]) 

                          }

sigma2 <- esp.tr(sigma2.star, x$VC$margins[2])$vrb   
    
#############  
# NUs
#############    
  
if( x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3)  ) nu.st <- bs2[,  x$X2.d2 + 2] # t(as.matrix(bs2[,  x$X2.d2 + 2]))
  
  if( !is.null(x$X3) ){
  
if(!missing(newdata)){ X4s <- predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix")}                      
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}     
  
             nu.st <- X4s%*%t(bs2[,(x$X2.d2 + x$X3.d2 + 1):(x$X2.d2 + x$X3.d2 + x$X4.d2)]) 
                      }
                      
 nu <- esp.tr(nu.st, x$VC$margins[2])$vrb   
  
} 


#################


if( is.null(x$X3) ){

sigma2   <- matrix(rep(sigma2, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
nu       <- matrix(rep(nu, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)

                   }



if(y1 == 0 || y1 == 1){                              
p0s  <- 1 - p1s                     
p2s  <- distrHsAT(y2, eta2s, sigma2, nu, x$margins[2])$p2
p12s <- p0s*p2s

if(cond == 1) p12s <- p2s
if(cond == 2) p12s <- p0s

                      }

if(y1 == 1){                                          

p12s <- p1s*p2s 

if(cond == 1) p12s <- p2s
if(cond == 2) p12s <- p1s

           }

           
           

}# int
         
                          } # indep


} ## cont no




##################################################################################################
##################################################################################################













if(x$VC$Cont == "NO" && x$margins[2] %in% bin.link){ 


if(type == "bivariate"){


if(!missing(newdata)){

p1 <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "response")
p2 <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "response")

if(!is.null(x$X3)) theta <- teta.tr(x$VC, predict.SemiParBIVProbit(x, eq = 3, newdata = newdata))$teta
if(is.null(x$X3))  theta  <- x$theta 

}



if(missing(newdata)){

p1 <- x$p1
p2 <- x$p2

theta  <- x$theta 

}



p12 <- BiCDF(p1, p2, x$nC, theta) 



if(y1 == 1 && y2 == 1){ 

  p12 <- p12 
  if(cond == 1) p12 <- p12/p1
  if(cond == 2) p12 <- p12/p2

                       }

if(y1 == 1 && y2 == 0){ 

  p12 <- pmax(p1 - p12, epsilon)      
  if(cond == 1) p12 <- p12/p1
  if(cond == 2) p12 <- p12/(1-p2)

                        }

if(y1 == 0 && y2 == 1){ 

  p12 <- pmax(p2 - p12, epsilon)      
  if(cond == 1) p12 <- p12/(1-p1)
  if(cond == 2) p12 <- p12/p2

                      }

if(y1 == 0 && y2 == 0){ 

  p12 <- pmax(1 - p12 - pmax(p1 - p12, epsilon) - pmax(p2 - p12, epsilon) , epsilon)      
  if(cond == 1) p12 <- p12/(1-p1)
  if(cond == 2) p12 <- p12/(1-p2)

                      }




if(intervals == TRUE){

bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)  

#############  
# etas
############# 


if(!missing(newdata)){ X1  <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix") 
                       X2s <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix") }
                       
if( missing(newdata)){ X1 <- x$X1 
                       if(x$Model == "BSS") X2s <- x$X2s else X2s <- x$X2  
                       } 


p1s <- probm( X1%*%t(bs[,1:x$X1.d2]), x$VC$margins[1])$pr 
p2s <- probm(X2s%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]) , x$VC$margins[2])$pr

#############  
# thetas
#############  

if(x$Model != "BPO0"){

if(is.null(x$X3))  epds <- bs[,length(x$coefficients)]


if(!is.null(x$X3)){

if(!missing(newdata)){ X3s <- predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix") }                  
if( missing(newdata)){ if(x$Model == "BSS") X3s <- x$X3s else X3s <- x$X3 } 

epds <- X3s%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)])

                   }

            
est.RHOb <- teta.tr(x$VC, epds)$teta

}

if(x$Model == "BPO0") est.RHOb <- rep(0, n.sim )
  
if( is.null(x$X3) ) est.RHOb <- matrix(rep(est.RHOb, each = dim(p1s)[1]), ncol = n.sim, byrow=FALSE)


#############

if(x$VC$BivD == "N") p12s <- matrix(BiCDF(p1s, p2s, x$nC, est.RHOb, test = FALSE), dim(p1s)[1], n.sim) else p12s <- BiCDF(p1s, p2s, x$nC, est.RHOb, test = FALSE)

if(y1 == 1 && y2 == 1){ 

  if(cond == 1) p12s <- p12s/p1s
  if(cond == 2) p12s <- p12s/p2s

                      }

if(y1 == 1 && y2 == 0){ 

  p12s <- pmax(p1s - p12s, epsilon)     
  if(cond == 1) p12s <- p12s/p1s
  if(cond == 2) p12s <- p12s/(1-p2s)

                      }

if(y1 == 0 && y2 == 1){ 

  p12s <- pmax(p2s - p12s, epsilon)      
  if(cond == 1) p12s <- p12s/(1-p1s)
  if(cond == 2) p12s <- p12s/p2s

                      }

if(y1 == 0 && y2 == 0){ 

  p12s <- pmax(1 - p12s - pmax(p1s - p12s, epsilon) - pmax(p2s - p12s, epsilon) , epsilon)      
  if(cond == 1) p12s <- p12s/(1-p1s)
  if(cond == 2) p12s <- p12s/(1-p2s)

                       }

          

}




     
}







if(type == "independence"){



if(!missing(newdata)){

p1 <- probm( predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1])$pr
p2 <- probm( predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%x$gam2$coefficients, x$VC$margins[2])$pr

}



if(missing(newdata)){

p1 <- probm( predict.SemiParBIVProbit(x, eq = 1, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1])$pr
p2 <- probm( predict.SemiParBIVProbit(x, eq = 2, type = "lpmatrix")%*%x$gam2$coefficients, x$VC$margins[2])$pr


}



if(y1 == 1 && y2 == 1){ 

  p12 <- p1*p2     
  if(cond == 1) p12 <- p2
  if(cond == 2) p12 <- p1

                      }

if(y1 == 1 && y2 == 0){ 

  p12 <- p1*(1-p2)     
  if(cond == 1) p12 <- 1-p2
  if(cond == 2) p12 <- p1

                      }

if(y1 == 0 && y2 == 1){ 

  p12 <- (1-p1)*p2     
  if(cond == 1) p12 <- p2
  if(cond == 2) p12 <- 1-p1

                      } 

if(y1 == 0 && y2 == 0){ 

  p12 <- (1-p1)*(1-p2)      
  if(cond == 1) p12 <- (1-p2)
  if(cond == 2) p12 <- (1-p1)

                       }




    
    
if(intervals == TRUE){


bs1 <- rMVN(n.sim, mean = x$gam1$coefficients, sigma=x$gam1$Vp)
bs2 <- rMVN(n.sim, mean = x$gam2$coefficients, sigma=x$gam2$Vp)


if(!missing(newdata)){

p1s <- probm( predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%t(bs1[,1:x$X1.d2]), x$VC$margins[1])$pr
p2s <- probm( predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%t(bs2[,1:x$X2.d2]), x$VC$margins[2])$pr

}



if(missing(newdata)){

p1s <- probm( predict.SemiParBIVProbit(x, eq = 1, type = "lpmatrix")%*%t(bs1[,1:x$X1.d2]), x$VC$margins[1])$pr
p2s <- probm( predict.SemiParBIVProbit(x, eq = 2, type = "lpmatrix")%*%t(bs2[,1:x$X2.d2]), x$VC$margins[2])$pr


}



if(y1 == 1 && y2 == 1){ 

  p12s <- p1s*p2s     
  if(cond == 1) p12s <- p2s
  if(cond == 2) p12s <- p1s

                      }

if(y1 == 1 && y2 == 0){ 

  p12s <- p1s*(1-p2s)     
  if(cond == 1) p12s <- 1-p2s
  if(cond == 2) p12s <- p1s

                      }

if(y1 == 0 && y2 == 1){ 

  p12s <- (1-p1s)*p2s     
  if(cond == 1) p12s <- p2s
  if(cond == 2) p12s <- 1-p1s

                      }

if(y1 == 0 && y2 == 0){ 

  p12s <- (1-p1s)*(1-p2s)      
  if(cond == 1) p12s <- (1-p2s)
  if(cond == 2) p12s <- (1-p1s)

                      }



}    # inter

} # indep
    
    
    
  


} ## cont no and binary binary case

























if(intervals == TRUE){

CIp12 <- rowQuantiles(p12s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)


if(length(p12) > 1)  {res <- data.frame(p12, CIp12, p1, p2);       names(res)[2:3] <- names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2)))}
if(length(p12) == 1) {res <- data.frame(t(c(p12, CIp12, p1, p2))); names(res) <- c("p12",names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2))),"p1","p2")}


#names(res)[2:3] <- names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2)))

}else{

res <- data.frame(p12, p1, p2)

}


return(res)















}



