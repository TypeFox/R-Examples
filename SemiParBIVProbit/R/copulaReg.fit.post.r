copulaReg.fit.post <- function(SemiParFit, VC, qu.mag=NULL, 
                                      gam1, gam2, gam3, gam4, gam5, gam6, gam7){

Ve <- R <- theta <- edf <- edf1 <- NULL
theta.a <- sigma21 <- sigma22 <- sigma2.1.a <- sigma2.2.a <- nu1 <- nu2 <- nu1.a <- nu2.a <- NULL

cont2par  <- VC$m2 
cont3par  <- VC$m3 


He <- SemiParFit$fit$hessian
if(VC$margins[1] != "LN" && VC$margins[2] != "LN") logLik <- -SemiParFit$fit$l
if(VC$margins[1] == "LN" || VC$margins[2] == "LN") logLik <- -SemiParFit$fit$l.ln


epsilon <- 0.0000001 
max.p   <- 0.9999999
                                                                                                         
    He.eig <- eigen(He, symmetric=TRUE)
    if(min(He.eig$values) < sqrt(.Machine$double.eps) && sign( min( sign(He.eig$values) ) ) == -1) He.eig$values <- abs(He.eig$values)  
    if(min(He.eig$values) < sqrt(.Machine$double.eps) ) { pep <- which(He.eig$values < sqrt(.Machine$double.eps)); He.eig$values[pep] <- epsilon }
    Vb <- He.eig$vectors%*%tcrossprod(diag(1/He.eig$values),He.eig$vectors)  # this could be taken from magic as well 
    Vb <- (Vb + t(Vb) ) / 2 
    
                                     
if( (VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0 || VC$l.sp4!=0 || VC$l.sp5!=0 || VC$l.sp6!=0 || VC$l.sp7!=0) && VC$fp==FALSE){

    HeSh <- He - SemiParFit$fit$S.h
    F <- Vb%*%HeSh # diag(SemiParFit$magpp$edf)   # this could be taken from magic as well
    F1 <- diag(SemiParFit$magpp$edf1)             # needed for testing
    R <- SemiParFit$bs.mgfit$R                    # needed for testing
    Ve <- F%*%Vb                                  # diag(SemiParFit$magpp$Ve) and diag(SemiParFit$magpp$Vb) but need to be careful with dispersion parameters
                                          
}else{ 

HeSh <- He
Ve <- Vb
F <- F1 <- diag(rep(1,dim(Vb)[1]))
R <- SemiParFit$bs.mgfit$R

} 

t.edf <- sum(diag(F))

dimnames(SemiParFit$fit$hessian)[[1]] <- dimnames(SemiParFit$fit$hessian)[[2]] <- dimnames(Vb)[[1]] <- dimnames(Vb)[[2]] <- dimnames(HeSh)[[1]] <- dimnames(HeSh)[[2]] <- dimnames(F)[[1]] <- dimnames(F)[[2]] <- dimnames(He)[[1]] <- dimnames(He)[[2]] <- names(SemiParFit$fit$argument)   


############################################################################################
# SIGMAs
############################################################################################  
  
sigma21 <- esp.tr(SemiParFit$fit$etas1, VC$margins[1])$vrb  
sigma22 <- esp.tr(SemiParFit$fit$etas2, VC$margins[2])$vrb    
  
  

if( is.null(VC$X3) ) {names(sigma21) <- "sigma21"
                      names(sigma22) <- "sigma22" } 
  
  
  
  
  sigma2.1.a <- mean(sigma21); sigma2.2.a <- mean(sigma22)  
  if( length(sigma21)==1 ) sigma2.1.a <- sigma21
  if( length(sigma22)==1 ) sigma2.2.a <- sigma22  

############################################################################################
############################################################################################


############################################################################################
# NUs
############################################################################################  



if(VC$margins[1] %in% cont3par ){  


nu1 <- esp.tr(SemiParFit$fit$etan1, VC$margins[1])$vrb    
if( is.null(VC$X3) ) names(nu1) <- "nu1"

nu1.a <- mean(nu1)
if( length(nu1)==1 ) nu1.a <- nu1  

}




if(VC$margins[2] %in% cont3par ){  

nu2 <- esp.tr(SemiParFit$fit$etan2, VC$margins[2])$vrb    
if( is.null(VC$X3) ) names(nu2) <- "nu2"
                       
nu2.a <- mean(nu2)
if( length(nu2)==1 ) nu2.a <- nu2  

}


############################################################################################
# THETA
############################################################################################

dep <- SemiParFit$fit$etad
if( is.null(VC$X3) ) names(dep) <- "theta"
  

resT  <- teta.tr(VC, dep)
theta <- resT$teta        
        
theta.a  <- mean(theta) 
if( length(theta)==1 ) theta.a <- theta 

######################
# Association measures
######################

if(!(VC$BivD %in% c("AMH","FGM"))) tau <- BiCopPar2Tau(family = VC$nCa, par = theta)
if(VC$BivD == "AMH")               tau <- 1 - (2/3)/theta^2*(theta + (1-theta)^2*log(1-theta))
if(VC$BivD == "FGM")               tau <- 2/9*theta 
tau.a  <- mean(tau) 
if( length(tau)==1 ) tau.a <- tau 

######################

if(VC$gc.l == TRUE) gc()  


if( (VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0 || VC$l.sp4!=0 || VC$l.sp5!=0 || VC$l.sp6!=0 || VC$l.sp7!=0) ){

  edf <- edf1 <- list(0, 0, 0, 0, 0, 0, 0)
        
     for(i in 1:7){

       if(i==1) {mm <- VC$l.sp1; if(mm==0) next}
       if(i==2) {mm <- VC$l.sp2; if(mm==0) next} 
       if(i==3) {mm <- VC$l.sp3; if(mm==0) next} 
       if(i==4) {mm <- VC$l.sp4; if(mm==0) next} 
       if(i==5) {mm <- VC$l.sp5; if(mm==0) next}        
       if(i==6) {mm <- VC$l.sp6; if(mm==0) next} 
       if(i==7) {mm <- VC$l.sp7; if(mm==0) break}        

          for(k in 1:mm){

              if(i==1){ gam <- gam1; ind <- (gam$smooth[[k]]$first.para):(gam$smooth[[k]]$last.para) } 
              if(i==2){ gam <- gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 } 
              if(i==3){ gam <- gam3; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 + VC$X2.d2 } 
              if(i==4){ gam <- gam4; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 + VC$X2.d2 + VC$X3.d2 } 
              if(i==5){ gam <- gam5; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 + VC$X2.d2 + VC$X3.d2  + VC$X4.d2 } 
              if(i==6){ gam <- gam6; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 + VC$X2.d2 + VC$X3.d2  + VC$X4.d2 + VC$X5.d2 } 
              if(i==7){ gam <- gam7; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 + VC$X2.d2 + VC$X3.d2  + VC$X4.d2 + VC$X5.d2 + VC$X6.d2} 
              
              
	      edf[[i]][k]  <-  sum(diag(F)[ind])
	      edf1[[i]][k] <- sum(diag(F1)[ind])
                        }
                  }
         
  if(VC$l.sp1!=0) names(edf[[1]]) <- names(edf1[[1]]) <- names(gam1$sp)  
  if(VC$l.sp2!=0) names(edf[[2]]) <- names(edf1[[2]]) <- names(gam2$sp)   
  if(VC$l.sp3!=0) names(edf[[3]]) <- names(edf1[[3]]) <- names(gam3$sp)  
  if(VC$l.sp4!=0) names(edf[[4]]) <- names(edf1[[4]]) <- names(gam4$sp) 
  if(VC$l.sp5!=0) names(edf[[5]]) <- names(edf1[[5]]) <- names(gam5$sp) 
  if(VC$l.sp6!=0) names(edf[[6]]) <- names(edf1[[6]]) <- names(gam6$sp) 
  if(VC$l.sp7!=0) names(edf[[7]]) <- names(edf1[[7]]) <- names(gam7$sp)  

}
  
 
  sp <- SemiParFit$sp 
   
    
                 list(SemiParFit = SemiParFit, He = He, logLik = logLik, Vb = Vb, HeSh = HeSh, 
                      F = F, F1 = F1, t.edf = t.edf, edf = edf, 
                      edf11=edf1,
                      edf1 = edf[[1]], edf2 = edf[[2]], edf3 = edf[[3]], edf4 = edf[[4]], 
                      edf5 = edf[[5]], edf6 = edf[[6]], edf7 = edf[[7]],
                      edf1.1 = edf1[[1]], edf1.2 = edf1[[2]], edf1.3 = edf1[[3]], edf1.4 = edf1[[4]], 
                      edf1.5 = edf1[[5]], edf1.6 = edf1[[6]], edf1.7 = edf1[[7]],
                      theta = theta, theta.a = theta.a, tau = tau, tau.a= tau.a,
                      sigma21 = sigma21, sigma22 = sigma22, sigma21.a = sigma2.1.a, sigma22.a = sigma2.2.a,
                      nu1 = nu1, nu1.a = nu1.a, nu2= nu2, nu2.a = nu2.a,
                      sp = sp, R = R, Ve = Ve) 

}



