summary.copulaReg <- function(object, n.sim = 100, prob.lev = 0.05, 
                                     cm.plot = FALSE, ylab = "Margin 2", 
                                     xlab = "Margin 1", 
                                     n.grid = 1000, n.dig = 2, ...){


  testStat <- getFromNamespace("testStat", "mgcv")
  liu2   <- getFromNamespace("liu2", "mgcv") 
  filled.contour <- getFromNamespace("filled.contour", "graphics")      

  bs <- SE <- Vb <- epds <- sigma2.st <- sigma2 <- nu.st <- nu <- est.RHOb <- XX <- Xt <- V <- 1
  
cont2par  <- object$VC$m2 
cont3par  <- object$VC$m3 
  
  n <- object$n
  
  tableN <- table <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  CIkt <- CIsig21 <- CIsig22 <- CInu1 <- CInu2 <- CIrs <- NULL  
  epsilon <- 0.0000001 
  max.p   <- 0.9999999
  
  est.RHOb <- rep(NA,n.sim) 

  lf <- length(object$coefficients)
  Vb <- object$Vb 
  SE <- sqrt(diag(Vb)) 
  bs <- rMVN(n.sim, mean = object$coefficients, sigma=Vb)  


#######
# theta
#######

  if(  is.null(object$X3) ) epds <- bs[,lf]
  
  	if( !is.null(object$X3) ){ 
  if(object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont2par) epds <- object$X5%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2)])
 if((object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont2par) || (object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont3par) ) epds <- object$X6%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2)])
  if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont3par) epds <- object$X7%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2)])
  	                          }

   
   

est.RHOb <- teta.tr(object$VC, epds)$teta
  
   
   if( is.null(object$X3) ) est.RHOb <- t(as.matrix(est.RHOb))
   
   CIrs <- rowQuantiles(est.RHOb, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ) CIrs <- t(CIrs) 


######################
# Association measures
######################

if(!(object$VC$BivD %in% c("AMH","FGM"))) tau <- BiCopPar2Tau(family = object$VC$nCa, par = est.RHOb) 
if(object$VC$BivD == "AMH")               tau <- 1 - (2/3)/est.RHOb^2*(est.RHOb + (1-est.RHOb)^2*log(1-est.RHOb))
if(object$VC$BivD == "FGM")               tau <- 2/9*est.RHOb 

CIkt <- rowQuantiles(tau, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)

if( is.null(object$X3) ) CIkt <- t(CIkt) 


#############  
# sigmas
#############  



      if( is.null(object$X3) ) {
      
if(object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont2par ){ ps1 <- lf - 2; ps2 <- lf - 1 }
if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont3par ){ ps1 <- lf - 4; ps2 <- lf - 3 }
if((object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont3par) || (object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont2par) ){ ps1 <- lf - 3; ps2 <- lf - 2 }
      
      
                                sigma2.1.star <- t(as.matrix(bs[, ps1]))
                                sigma2.2.star <- t(as.matrix(bs[, ps2]))
                                
                                }
  
  
      if( !is.null(object$X3) ) {

       sigma2.1.star <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)]) 
       sigma2.2.star <- object$X4%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)]) 

                                }  

    sigma21 <- esp.tr(sigma2.1.star, object$VC$margins[1])$vrb  
    sigma22 <- esp.tr(sigma2.2.star, object$VC$margins[2])$vrb  

   CIsig21 <- rowQuantiles(sigma21, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   CIsig22 <- rowQuantiles(sigma22, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ){ CIsig21 <- t(CIsig21); CIsig22 <- t(CIsig22) }
  
  
  
#############  
# NUs
#############    
  
  
if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont3par ){  
    
  if( is.null(object$X3) )  {    pn1 <- lf - 2 
                                 pn2 <- lf - 1
                                 nu1.st <- t(as.matrix(bs[, pn1]))
                                 nu2.st <- t(as.matrix(bs[, pn2]))  } 
  
  if( !is.null(object$X3) ) {  
 
       nu1.st <- object$X5%*%t(bs[,(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + 1):(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2)]) 
       nu2.st <- object$X6%*%t(bs[,(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + 1):(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + object$X6.d2)])
  
                }   
   
nu1 <- esp.tr(nu1.st, object$VC$margins[1])$vrb   
nu2 <- esp.tr(nu2.st, object$VC$margins[2])$vrb   


   CInu1 <- rowQuantiles(nu1, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   CInu2 <- rowQuantiles(nu2, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   
      if( is.null(object$X3) ){ CInu1 <- t(CInu1); CInu2 <- t(CInu2) }

   
} 


if(object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont3par ){  
  

  if( is.null(object$X3) )  {  pn2 <- lf - 1; nu2.st <- t(as.matrix(bs[, pn2]))  } 
  
  if( !is.null(object$X3) ) {  
 
       nu2.st <- object$X5%*%t(bs[,(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + 1):(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2)]) 
  
                }   
  
nu2 <- esp.tr(nu2.st, object$VC$margins[2])$vrb   
 
   
   CInu2 <- rowQuantiles(nu2, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ) CInu2 <- t(CInu2) 

   
} 
  
  
  
if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont2par ){  
    
  if( is.null(object$X3) )  {  pn1 <- lf - 1; nu1.st <- t(as.matrix(bs[, pn1]))  } 
  
  if( !is.null(object$X3) ) {  
 
       nu1.st <- object$X5%*%t(bs[,(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + 1):(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2)]) 
  
                }   
  
nu1 <- esp.tr(nu1.st, object$VC$margins[1])$vrb  

   CInu1 <- rowQuantiles(nu1, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ) CInu1 <- t(CInu1)

   
}  

   
#########################
         

  if(object$VC$gc.l == TRUE) gc()

  index <- 1:2
  ind1 <- 1:object$gp1
  ind2 <- object$X1.d2 + (1:object$gp2)
  ind3 <- ind4 <- ind5 <- ind6 <- ind7 <- NULL 
  
  if(!is.null(object$X3) ) {
  
       ind3 <- object$X1.d2 + object$X2.d2 + (1:object$gp3)
       index <- 1:3
       
       if(!is.null(object$X4) ) {
       ind4 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + (1:object$gp4)
       index <- 1:4
       }
                                
       if(!is.null(object$X5) ) {
       ind5 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + (1:object$gp5)
       index <- 1:5  
       }     
                                
       if(!is.null(object$X6) ) {
       ind6 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + (1:object$gp6)
       index <- 1:6  
       }  
       
       if(!is.null(object$X7) ) {
       ind7 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + object$X6.d2 + (1:object$gp7)
       index <- 1:7  
       }        
                            
  }
                            
  ind <- list( ind1 = ind1,
               ind2 = ind2,
               ind3 = ind3, 
               ind4 = ind4,
               ind5 = ind5,
               ind6 = ind6,
               ind7 = ind7)
                

  for(i in index){
  estimate <- coef(object)[ind[[i]]]
  se       <- SE[ind[[i]]]
  ratio    <- estimate/se
  pv       <- 2*pnorm(abs(ratio), lower.tail = FALSE)
  table[[i]] <- cbind(estimate,se,ratio,pv)
  dimnames(table[[i]])[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }

  
  if( object$l.sp1!=0 || object$l.sp2!=0 || object$l.sp3!=0 || object$l.sp4!=0 || object$l.sp5!=0 || object$l.sp6!=0 || object$l.sp7!=0){

  	pTerms.df <- pTerms.chi.sq <- pTerms.pv <- tableN <- list(0, 0, 0, 0, 0, 0, 0)
        XX <- object$R
        
           for(i in index){

             if(i==1) {mm <- object$l.sp1; if(mm==0) next}
             if(i==2) {mm <- object$l.sp2; if(mm==0) next} 
             if(i==3) {mm <- object$l.sp3; if(mm==0) next} 
             if(i==4) {mm <- object$l.sp4; if(mm==0) next} 
             if(i==5) {mm <- object$l.sp5; if(mm==0) next} 
             if(i==6) {mm <- object$l.sp6; if(mm==0) next} 
             if(i==7) {mm <- object$l.sp7; if(mm==0) break} 
  
		for(k in 1:mm){

                        if(i==1){ gam <- object$gam1; ind <-  gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para                                } 
                        if(i==2){ gam <- object$gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2                } 
                        if(i==3){ gam <- object$gam3; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 }
                        if(i==4){ gam <- object$gam4; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 }
                        if(i==5){ gam <- object$gam5; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 }
                        if(i==6){ gam <- object$gam6; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 }
                        if(i==7){ gam <- object$gam7; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + object$X6.d2 }
                          
                          
                        gam$sig2            <- 1
                        gam$scale.estimated <- FALSE                          
                          
                        if(gam$smooth[[k]]$null.space.dim == 0){
                        
                        LRB <- rbind(XX, t(mroot(object$fit$S.h)))
			LRB <- cbind(LRB[, -ind], LRB[, ind])
			ind1 <- (ncol(LRB) - length(ind) + 1):ncol(LRB)
			Rm <- qr.R(qr(LRB, tol = 0, LAPACK = FALSE))[ind1, ind1]
                        B <- mroot(object$Ve[ind, ind, drop = FALSE])
                          
			b.hat <- coef(object)[ind]
			d <- Rm %*% b.hat
			stat <- sum(d^2)
			ev <- eigen(crossprod(Rm %*% B), symmetric = TRUE, only.values = TRUE)$values
			ev[ev < 0] <- 0
			rank <- sum(ev > max(ev) * .Machine$double.eps^0.8)
			pval <- liu2(stat, ev)                          
                        Tp <- list(stat = stat, pval = pval, rank = rank)  
                          
                        }
                          
			if(gam$smooth[[k]]$null.space.dim != 0){
			
			b  <- coef(object)[ind]
			V  <- Vb[ind,ind, drop = FALSE]
			Xt <- XX[, ind, drop = FALSE] 
			pTerms.df[[i]][k] <- min(ncol(Xt), object$edf11[[i]][k])
			Tp <- testStat(b, Xt, V, pTerms.df[[i]][k], type = 0, res.df = -1)
			
			}
			
			
			pTerms.chi.sq[[i]][k] <- Tp$stat 
			pTerms.df[[i]][k] <- Tp$rank
                        pTerms.pv[[i]][k] <- Tp$pval
			                 
                }
                
              tableN[[i]] <- cbind(object$edf[[i]], pTerms.df[[i]], pTerms.chi.sq[[i]], pTerms.pv[[i]])
              dimnames(tableN[[i]])[[2]] <- c("edf", "Ref.df", "Chi.sq", "p-value")
              
            }

  if(object$VC$gc.l == TRUE) gc()

  }
  



 
if(cm.plot == TRUE){
  
  c1   <- c("N","GU","rGU","LO")
  c2   <- c("LN","WEI","iG","GA","GAi","DAGUM","SM","BE","FISK") 
   
  m2 <- s2 <- nu2 <- m1 <- s1 <- nu1 <- 0 
  par1 <- object$theta.a 
   
  if(object$margins[2] %in% cont2par ){ 
  m2 <- mean(object$eta2) 
  s2 <- object$sigma22.a
  nu2 <- 0
                                      } 
  
  if(object$margins[2] %in% cont3par ){ 
  m2 <- mean(object$eta2) 
  s2 <- object$sigma22.a
  nu2 <- object$nu2.a
                                      } 
                                      
  if(object$margins[1] %in% cont2par ){ 
  m1 <- mean(object$eta1) 
  s1 <- object$sigma21.a
  nu1 <- 0
                                      } 
  
  if(object$margins[1] %in% cont3par ){ 
  m1 <- mean(object$eta1) 
  s1 <- object$sigma21.a
  nu1 <- object$nu1.a
                                      }                                       
                                      
                                      
                                      
  if(object$BivD=="AMH")  {cop <- bquote(paste("AMH (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="FGM")  {cop <- bquote(paste("FGM (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="N")    {cop <- bquote(paste("Gaussian (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="F")    {cop <- bquote(paste("Frank (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="C0")   {cop <- bquote(paste("Clayton (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="C90")  {cop <- bquote(paste("90",degree," Clayton (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="C180") {cop <- bquote(paste("180",degree, " Clayton (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="C270") {cop <- bquote(paste("270",degree, " Clayton (",hat(theta)," = ",.(round(par1,2)),")",sep=""))} 
  if(object$BivD=="J0")   {cop <- bquote(paste("Joe (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="J90")  {cop <- bquote(paste("90",degree," Joe (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="J180") {cop <- bquote(paste("180",degree," Joe (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="J270") {cop <- bquote(paste("270",degree," Joe (",hat(theta)," = ",.(round(par1,2)),")",sep=""))} 
  if(object$BivD=="G0")   {cop <- bquote(paste("Gumbel (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="G90")  {cop <- bquote(paste("90",degree," Gumbel (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="G180") {cop <- bquote(paste("180",degree," Gumbel (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
  if(object$BivD=="G270") {cop <- bquote(paste("270",degree," Gumbel (",hat(theta)," = ",.(round(par1,2)),")",sep=""))} 
  
  
  
  
  
 
   Cplot <- function (BivD, par1, mar1, mar2, m1, m2, s1, s2, nu1, nu2, resp1, resp2, ...){   
   
     size <- 100
                                         
           x1 <- seq(from = min(resp1), to = max(resp1), length.out = size)                                                      
           x2 <- seq(from = min(resp2), to = max(resp2), length.out = size)                        
            
           x11 <- rep(x1, each = size)
           x22 <- rep(x2, times = size)
 
           resf1 <- distrHsAT(x11, m1, s1, nu1, mar1)
           d.x1  <- resf1$pdf2
           p1    <- resf1$p2
           
           
           resf2 <- distrHsAT(x22, m2, s2, nu2, mar2)
           d.x2  <- resf2$pdf2
           p2    <- resf2$p2
           
           
           md <- copgHsAT(p1, p2, par1, BivD, Ln = TRUE)$c.copula2.be1be2*d.x1*d.x2    
           z  <- matrix(data = md, nrow = size, byrow = TRUE)
       
           filled.contour(x1, x2, z, color = topo.colors, nlevels = 16, ...) 
                
  }
  
 
 
   
   if(object$margins[1] %in% c1) test.resp <- seq(from = min(object$y1), to = max(object$y1), length.out = n.grid)
   if(object$margins[1] %in% c2) test.resp <- seq(from = 0.0001, to = max(object$y1), length.out = n.grid)
 
   
   test.dens <- round( distrHsAT(test.resp, m1, s1, nu1, object$VC$margins[1])$pdf2, n.dig)
   resp1 <- test.resp[which(test.dens != 0)]
   
   if( length(resp1) < 2 ) stop("Increase/decrease n.dig and/or n.grid values.")
   
   
   if(object$margins[2] %in% c1) test.resp <- seq(from = min(object$y2), to = max(object$y2), length.out = n.grid)
   if(object$margins[2] %in% c2) test.resp <- seq(from = 0.0001, to = max(object$y2), length.out = n.grid)
 
   
   test.dens <- round( distrHsAT(test.resp, m2, s2, nu2, object$VC$margins[2])$pdf2, n.dig)
   resp2 <- test.resp[which(test.dens != 0)]
   
   if( length(resp2) < 2 ) stop("Increase/decrease n.dig and/or n.grid values.")
   
  
 Cplot(BivD = object$BivD, par1 = par1, mar1 = object$margins[1], mar2 = object$margins[2], 
       m1 = m1, m2 = m2, s1 = s1, s2 = s2, nu1 = nu1, nu2 = nu2, resp1 = resp1, resp2 = resp2, main = cop, ylab = ylab, xlab = xlab, ...)  
       
       
       
  }
  



 
 
rm(bs, SE, Vb, est.RHOb, XX, Xt, V) 
 
  res <- list(tableP1=table[[1]], tableP2=table[[2]], tableP3=table[[3]], 
              tableP4=table[[4]], tableP5=table[[5]], tableP6=table[[6]], tableP7=table[[7]],
              tableNP1=tableN[[1]], tableNP2=tableN[[2]], tableNP3=tableN[[3]], 
              tableNP4=tableN[[4]], tableNP5=tableN[[5]], tableNP6=tableN[[6]], tableNP7=tableN[[7]], 
              n=n, theta=object$theta, theta.a=object$theta.a, 
              sigma21=object$sigma21, sigma22=object$sigma22, 
              nu1=object$nu1, nu2=object$nu2, tau=object$tau, 
              sigma21.a=object$sigma21.a, sigma22.a=object$sigma22.a, 
              nu1.a=object$nu1.a, nu2.a=object$nu2.a, 
              tau.a=object$tau.a, formula = object$formula,
              formula1=object$gam1$formula, formula2=object$gam2$formula, formula3=object$gam3$formula,
              formula4=object$gam4$formula, formula5=object$gam5$formula, formula6=object$gam6$formula, formula7=object$gam7$formula,
              t.edf=object$t.edf, 
              CItheta=CIrs, CIsig21=CIsig21, CIsig22=CIsig22, CInu1=CInu1, CInu2=CInu2, CIkt = CIkt,
              BivD=object$BivD, margins = object$margins, 
              l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3, 
              l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6, l.sp7 = object$l.sp7,
              X3.null = is.null(object$X3)
              )
              
              
              
  class(res) <- "summary.copulaReg"
      
                                        

res

}

