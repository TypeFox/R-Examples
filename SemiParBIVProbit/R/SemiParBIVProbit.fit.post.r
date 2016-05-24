SemiParBIVProbit.fit.post <- function(SemiParFit, Model, VC, qu.mag=NULL, 
                                      gam1, gam2, gam3, gam4, gam5, gam6){

Ve <- R <- X2s <- lambda1s <- X3s <- lambda2s <- eta1S <- eta2S <- theta <- edf <- edf1 <- theta.a <- sigma2 <- sigma2.a <- OR <- GM <- p1n <- p2n <- nu <- nu.a <- NULL

cont2par  <- VC$m2 
cont3par  <- VC$m3 
bin.link  <- VC$bl  


He <- SemiParFit$fit$hessian
if(VC$margins[2] != "LN") logLik <- -SemiParFit$fit$l else logLik <- -SemiParFit$fit$l.ln 

epsilon <- 0.0000001 
max.p   <- 0.9999999
                                                                                                         
    He.eig <- eigen(He, symmetric=TRUE)
   if(min(He.eig$values) < sqrt(.Machine$double.eps) && sign( min( sign(He.eig$values) ) ) == -1) He.eig$values <- abs(He.eig$values)  
    if(min(He.eig$values) < sqrt(.Machine$double.eps) ) { pep <- which(He.eig$values < sqrt(.Machine$double.eps)); He.eig$values[pep] <- epsilon }
    Vb <- He.eig$vectors%*%tcrossprod(diag(1/He.eig$values),He.eig$vectors)  # this could be taken from magic as well 
    Vb <- (Vb + t(Vb) ) / 2 
    
                                     
if( (VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0 || VC$l.sp4!=0 || VC$l.sp5!=0 || VC$l.sp6!=0) && VC$fp==FALSE){

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

  
if(VC$hess == FALSE) SemiParFit$fit$Fisher <- SemiParFit$fit$hessian



##################################################################

if(VC$margins[2] %in% cont2par ){

sigma2 <- esp.tr(SemiParFit$fit$etas, VC$margins[2])$vrb 
if( is.null(VC$X4) ) names(sigma2) <- "sigma2"

  sigma2.a <- mean(sigma2) 
  if( length(sigma2)==1 ) sigma2.a <- sigma2  

}


#################################################################

if(VC$margins[2] %in% cont3par ){

if(VC$margins[2] %in% c("DAGUM","SM")){

sigma2 <- esp.tr(SemiParFit$fit$etas, VC$margins[2])$vrb 
nu     <- esp.tr(SemiParFit$fit$etan, VC$margins[2])$vrb 

if( is.null(VC$X4) && is.null(VC$X5) ) {names(sigma2) <- "sigma2"; names(nu) <- "nu"}
           
                           }  
                           
  sigma2.a <- mean(sigma2)
  nu.a     <- mean(nu)
  if( length(sigma2)==1 ) {sigma2.a <- sigma2; nu.a <- nu}  

}



############################################################################################
# THETA
############################################################################################


if(VC$Model == "BPO0" ) dep <- theta <- 0

if(!(VC$Model %in% c("BPO0","BSS")) ){

  dep <- SemiParFit$fit$etad
  if( is.null(VC$X3) ) names(dep) <- "theta" 

theta <- teta.tr(VC, dep)$teta   
        
theta.a  <- mean(theta)   

}


############################################################################################

   
  if(Model=="BSS"){

  SemiParFit$fit$eta2 <- VC$X2s%*%SemiParFit$fit$argument[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  
  p1n <- predict.gam(gam1, type="response")
  p2n <- probm(VC$X2s%*%coef(gam2), VC$margins[2])$pr
 

if(!is.null(VC$X3)) {

SemiParFit$fit$etad <- VC$X3s%*%SemiParFit$fit$argument[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]

theta <- teta.tr(VC, SemiParFit$fit$etad)$teta 

}


if(is.null(VC$X3)){

  dep <- SemiParFit$fit$etad
  names(dep) <- "theta" 
  theta <- teta.tr(VC, dep)$teta   
        
                  }
        
  theta.a <- mean(theta)   
                   
 

}



if(Model=="BSS" || Model=="BPO" || Model=="BPO0"){

  p1 <- probm(SemiParFit$fit$eta1, VC$margins[1])$pr 
  p2 <- probm(SemiParFit$fit$eta2, VC$margins[2])$pr # eta2 modified above
  
   p11 <- BiCDF(p1, p2, VC$nC, theta)
  
   SemiParFit$fit$p10 <- p1 - p11
   SemiParFit$fit$p11 <- p11
   SemiParFit$fit$p00 <- (1 - p2) - ( p1 - p11 )
   SemiParFit$fit$p01 <- p2 - p11
   
   SemiParFit$fit$p1 <- p1
   SemiParFit$fit$p2 <- p2

}


######################
# Association measures
######################

if(!(VC$BivD %in% c("AMH","FGM"))) tau <- BiCopPar2Tau(family = VC$nCa, par = theta)
if(VC$BivD == "AMH")               tau <- 1 - (2/3)/theta^2*(theta + (1-theta)^2*log(1-theta))
if(VC$BivD == "FGM")               tau <- 2/9*theta 
tau.a <- mean(tau) 
if( length(tau)==1 ) tau.a <- tau 


####

if(VC$margins[2] %in% bin.link){

p00 <- SemiParFit$fit$p00 
p01 <- SemiParFit$fit$p01 
p11 <- SemiParFit$fit$p11 
p10 <- SemiParFit$fit$p10 

p1 <- SemiParFit$fit$p1
p2 <- SemiParFit$fit$p2

OR <- (p00*p11)/(p01*p10)

OR  <- ifelse(OR  ==  Inf,  8.218407e+307, OR ) 
OR  <- ifelse(OR  == -Inf, -8.218407e+307, OR ) 

GM <- mean((OR - 1)/(OR + 1))
OR <- mean(OR)


rm(p00,p01,p10,p11,p1,p2)

}

if(VC$gc.l == TRUE) gc()  




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
  
  if( length(theta)==1 ) theta.a <- theta  
  
  
                 list(SemiParFit = SemiParFit, He = He, logLik = logLik, Vb = Vb, HeSh = HeSh, F = F, F1 = F1, t.edf = t.edf, edf = edf, 
                      edf11=edf1,
                      edf1 = edf[[1]], edf2 = edf[[2]], edf3 = edf[[3]], edf4 = edf[[4]], edf5 = edf[[5]], edf6 = edf[[6]],
                      edf1.1 = edf1[[1]], edf1.2 = edf1[[2]], edf1.3 = edf1[[3]], edf1.4 = edf1[[4]], edf1.5 = edf1[[5]], edf1.6 = edf1[[6]],
                      theta = theta, theta.a = theta.a, sigma2 = sigma2, sigma2.a = sigma2.a,
                      nu = nu, nu.a = nu.a, tau = tau, tau.a= tau.a,
                      sp = sp, OR = OR, GM = GM, p1n=p1n, p2n=p2n, R = R, Ve = Ve) 

}



