SemiParTRIVProbit.fit.post <- function(SemiParFit, VC, Model, 
                                       gam1, gam2, gam3, gam4, gam5, gam6){

Ve <- R <- X2s <- lambda1s <- X3s <- lambda2s <- eta1S <- eta2S <- theta <- edf <- edf1 <- theta.a <- sigma2 <- sigma2.a <- OR <- GM <- p1n <- p2n <- p3n <- nu <- nu.a <- NULL

He <- SemiParFit$fit$hessian
logLik <- -SemiParFit$fit$l

epsilon <- 0.0000001 
max.p   <- 0.9999999
                                                                                                         
    He.eig <- eigen(He, symmetric=TRUE)
    if(min(He.eig$values) < sqrt(.Machine$double.eps) && sign( min( sign(He.eig$values) ) ) == -1) He.eig$values <- abs(He.eig$values)  
    if(min(He.eig$values) < sqrt(.Machine$double.eps) ) { pep <- which(He.eig$values < sqrt(.Machine$double.eps)); He.eig$values[pep] <- epsilon }  
    Vb <- He.eig$vectors%*%tcrossprod(diag(1/He.eig$values),He.eig$vectors)  
    Vb <- (Vb + t(Vb) ) / 2 
    
                                     
if( (VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0 || VC$l.sp4!=0 || VC$l.sp5!=0 || VC$l.sp6!=0) && VC$fp==FALSE){

    HeSh <- He - SemiParFit$fit$S.h
    F <- Vb%*%HeSh 
    F1 <- diag(SemiParFit$magpp$edf1)             
    R <- SemiParFit$bs.mgfit$R                    
    Ve <- F%*%Vb                                  
                                          
}else{ 

HeSh <- He
Ve <- Vb
F <- F1 <- diag(rep(1,dim(Vb)[1]))
R <- SemiParFit$bs.mgfit$R

} 


t.edf <- sum(diag(F))

dimnames(SemiParFit$fit$hessian)[[1]] <- dimnames(SemiParFit$fit$hessian)[[2]] <- dimnames(Vb)[[1]] <- dimnames(Vb)[[2]] <- dimnames(HeSh)[[1]] <- dimnames(HeSh)[[2]] <- dimnames(F)[[1]] <- dimnames(F)[[2]] <- dimnames(He)[[1]] <- dimnames(He)[[2]] <- names(SemiParFit$fit$argument)   

if(VC$hess == FALSE) SemiParFit$fit$Fisher <- SemiParFit$fit$hessian

############################################################################################
# THETAs
############################################################################################

theta12 <- SemiParFit$fit$theta12     
theta13 <- SemiParFit$fit$theta13     
theta23 <- SemiParFit$fit$theta23  

names(theta12) <- "theta12"
names(theta13) <- "theta13" 
names(theta23) <- "theta23" 

theta12.a  <- mean(theta12) 
theta13.a  <- mean(theta13)
theta23.a  <- mean(theta23)   

# useful for when we will have predictors on correlations

############################################################################################




   
  if(Model=="TSS"){

  SemiParFit$fit$eta2 <- VC$X2s%*%SemiParFit$fit$argument[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  SemiParFit$fit$eta3 <- VC$X3s%*%SemiParFit$fit$argument[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
  
  
  p1n <- predict.gam(gam1, type="response")
  p2n <- probm(VC$X2s%*%coef(gam2), VC$margins[2])$pr 
  p3n <- probm(VC$X3s%*%coef(gam3), VC$margins[3])$pr 
 
  #p1 <- probm(SemiParFit$fit$eta1, VC$margins[1])$pr 
  #p2 <- probm(SemiParFit$fit$eta2, VC$margins[2])$pr 
  #p3 <- probm(SemiParFit$fit$eta3, VC$margins[2])$pr   
  
   #p11 <- BiCDF(p1, p2, VC$nC, theta)
  #
   #SemiParFit$fit$p10 <- p1 - p11
   #SemiParFit$fit$p11 <- p11
   #SemiParFit$fit$p00 <- (1 - p2) - ( p1 - p11 )
   #SemiParFit$fit$p01 <- p2 - p11
   #
   #SemiParFit$fit$p1 <- p1
   #SemiParFit$fit$p2 <- p2

}



if(VC$gc.l == TRUE) gc()  


VC$l.sp4 <- 0 # this needs to be amended if we have smooths on corrs

if( (VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0 || VC$l.sp4!=0 || VC$l.sp5!=0 || VC$l.sp6!=0) ){

  edf <- edf1 <- list(0, 0, 0, 0, 0, 0)
        
     for(i in 1:6){

       if(i==1) {mm <- VC$l.sp1; if(mm==0) next}
       if(i==2) {mm <- VC$l.sp2; if(mm==0) next} 
       if(i==3) {mm <- VC$l.sp3; if(mm==0) next} 
       if(i==4) {mm <- VC$l.sp4; if(mm==0) next} 
       if(i==5) {mm <- VC$l.sp5; if(mm==0) next}        
       if(i==6) {mm <- VC$l.sp6; if(mm==0) break} 

          for(k in 1:mm){

              if(i==1){ gam <- gam1; ind <- (gam$smooth[[k]]$first.para):(gam$smooth[[k]]$last.para) } 
              if(i==2){ gam <- gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 } 
              if(i==3){ gam <- gam3; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 + VC$X2.d2 } 
              if(i==4){ gam <- gam4; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 + VC$X2.d2 + VC$X3.d2 } 
              if(i==5){ gam <- gam5; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 + VC$X2.d2 + VC$X3.d2  + VC$X4.d2 } 
              if(i==6){ gam <- gam6; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 + VC$X2.d2 + VC$X3.d2  + VC$X4.d2 + VC$X5.d2 } 
              
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

}
  
 
  sp <- SemiParFit$sp 
  
                 list(SemiParFit = SemiParFit, He = He, logLik = logLik, Vb = Vb, HeSh = HeSh, F = F, F1 = F1, t.edf = t.edf, edf = edf, 
                      edf11=edf1,
                      edf1 = edf[[1]], edf2 = edf[[2]], edf3 = edf[[3]], edf4 = edf[[4]], edf5 = edf[[5]], edf6 = edf[[6]],
                      edf1.1 = edf1[[1]], edf1.2 = edf1[[2]], edf1.3 = edf1[[3]], edf1.4 = edf1[[4]], edf1.5 = edf1[[5]], edf1.6 = edf1[[6]],
                      theta12 = theta12, theta12.a = theta12.a, 
                      theta13 = theta13, theta13.a = theta13.a,
                      theta23 = theta23, theta23.a = theta23.a,
                      sp = sp, R = R, Ve = Ve,
                      p1n = p1n, p2n = p2n, p3n = p3n) 

}



