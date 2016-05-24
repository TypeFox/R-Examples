summary.copulaSampleSel <- function(object, n.sim = 100, prob.lev = 0.05, 
                                     cm.plot = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), 
                                     ylab = "Margin 2", xlab = "Margin 1", 
                                     n.grid = 1000, n.dig = 2, ...){


  testStat <- getFromNamespace("testStat", "mgcv")
  liu2   <- getFromNamespace("liu2", "mgcv") 
  filled.contour <- getFromNamespace("filled.contour", "graphics")      

  bs <- SE <- Vb <- epds <- sigma2.st <- sigma2 <- nu.st <- nu <- est.RHOb <- et1s <- et2s <- p1s <- p2s <- p11s <- p10s <- p00s <- p01s <- ORs <- GMs <- XX <- Xt <- V <- 1
  
cont2par  <- object$VC$m2 
cont3par  <- object$VC$m3  
bin.link  <- object$VC$bl  

  
  n <- object$n; n.sel <- object$n.sel
  
  
  
  tableN <- table <- list(NULL, NULL, NULL, NULL, NULL, NULL)
  CIkt <- CIor <- CIgm <- CIsig2 <- CInu <- NULL  
  epsilon <- 0.0000001 
  max.p   <- 0.9999999
  
  
  est.RHOb <- rep(NA,n.sim) 

 
  lf <- length(object$coefficients)
  Vb <- object$Vb 
  SE <- sqrt(diag(Vb)) 

  
  bs <- rMVN(n.sim, mean = object$coefficients, sigma=Vb)  
  
  
  
  if(object$VC$margins[2] %in% cont2par ){
  
  if( !is.null(object$X3) ) sigma2.st <- object$X3s%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)]) 
  if(  is.null(object$X3) ) sigma2.st <- bs[, lf-1]
  

   sigma2 <- esp.tr(sigma2.st, object$VC$margins[2])$vrb  
   
   if(  is.null(object$X3) ) sigma2 <- t(as.matrix(sigma2))
   
   CIsig2 <- rowQuantiles(sigma2, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ) CIsig2 <- t(CIsig2) 

  if( !is.null(object$X4) ) epds <- object$X4s%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2 + 1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)])
  if(  is.null(object$X4) ) epds <- bs[, lf]  
   
  } 
  
  
  
    if(object$VC$margins[2] %in% cont3par ){
    
    if( !is.null(object$X3) ) sigma2.st <- object$X3s%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)]) 
    if(  is.null(object$X3) ) sigma2.st <- bs[, lf - 2]
    
    if( !is.null(object$X4) ) nu.st     <- object$X4s%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)]) 
    if(  is.null(object$X4) ) nu.st     <- bs[, lf - 1]    
    
    sigma2 <- esp.tr(sigma2.st, object$VC$margins[2])$vrb  
 
     
     if(  is.null(object$X3) ) sigma2 <- t(as.matrix(sigma2))
     
     CIsig2 <- rowQuantiles(sigma2, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
     if( is.null(object$X3) ) CIsig2 <- t(CIsig2) 

     
     if(object$VC$margins[2] %in% c("DAGUM","SM")){
     
     nu <- esp.tr(nu.st, object$VC$margins[2])$vrb  

     
     }
     

     
     if(  is.null(object$X4) ) nu <- t(as.matrix(nu))
     
     CInu <- rowQuantiles(nu, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
     if( is.null(object$X3) ) CInu <- t(CInu) 

    if( !is.null(object$X5) ) epds <- object$X5s%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2 + 1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2)])
    if(  is.null(object$X5) ) epds <- bs[, lf]  
     
  }
  
  
  
  
   est.RHOb <- teta.tr(object$VC, epds)$teta 
 
   if(  is.null(object$X3) ) est.RHOb <- t(as.matrix(est.RHOb))
   
   CIrs <- rowQuantiles(est.RHOb, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ) CIrs <- t(CIrs) 


######################
# Association measure
######################

if(!(object$VC$BivD %in% c("AMH","FGM"))) tau <- BiCopPar2Tau(family = object$VC$nCa, par = est.RHOb) 
if(object$VC$BivD == "AMH")               tau <- 1 - (2/3)/est.RHOb^2*(est.RHOb + (1-est.RHOb)^2*log(1-est.RHOb))
if(object$VC$BivD == "FGM")               tau <- 2/9*est.RHOb 

CIkt <- rowQuantiles(tau, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)

if( is.null(object$X3) ) CIkt <- t(CIkt)  

#####################
  

  
  if(object$VC$gc.l == TRUE) gc()



  index <- 1:2
  ind1 <- 1:object$gp1
  ind2 <- object$X1.d2 + (1:object$gp2)
  ind3 <- ind4 <- ind5 <- ind6 <- NULL 
  
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
                            
  }
                            
  ind <- list( ind1 = ind1,
               ind2 = ind2,
               ind3 = ind3, 
               ind4 = ind4,
               ind5 = ind5,
               ind6 = ind6)
                

  for(i in index){
  estimate <- coef(object)[ind[[i]]]
  se       <- SE[ind[[i]]]
  ratio    <- estimate/se
  pv       <- 2*pnorm(abs(ratio), lower.tail = FALSE)
  table[[i]] <- cbind(estimate,se,ratio,pv)
  dimnames(table[[i]])[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }

  
  if( object$l.sp1!=0 || object$l.sp2!=0 || object$l.sp3!=0 || object$l.sp4!=0 || object$l.sp5!=0 || object$l.sp6!=0){

  	pTerms.df <- pTerms.chi.sq <- pTerms.pv <- tableN <- list(0, 0, 0, 0, 0, 0)
        XX <- object$R
        
           for(i in index){

             if(i==1) {mm <- object$l.sp1; if(mm==0) next}
             if(i==2) {mm <- object$l.sp2; if(mm==0) next} 
             if(i==3) {mm <- object$l.sp3; if(mm==0) next} 
             if(i==4) {mm <- object$l.sp4; if(mm==0) next} 
             if(i==5) {mm <- object$l.sp5; if(mm==0) next} 
             if(i==6) {mm <- object$l.sp6; if(mm==0) break} 
  
		for(k in 1:mm){

                        if(i==1){ gam <- object$gam1; ind <-  gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para                                } 
                        if(i==2){ gam <- object$gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2                } 
                        if(i==3){ gam <- object$gam3; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 }
                        if(i==4){ gam <- object$gam4; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 }
                        if(i==5){ gam <- object$gam5; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 }
                        if(i==6){ gam <- object$gam6; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 }
                          
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
  
 m2 <- s2 <- nu2 <- 0 
 par1 <- object$theta.a 
  
 if(object$margins[2] %in% cont2par ){ 
 m2 <- mean(object$eta2) 
 s2 <- object$sigma2.a
 nu2 <- 0
                                     } 
 
 if(object$margins[2] %in% cont3par ){ 
 m2 <- mean(object$eta2) 
 s2 <- object$sigma2.a
 nu2 <- object$nu.a
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
 
 
 
 
 

  Cplot <- function (BivD, par1, mar2, m2, s2, nu2, resp, margins, ...){   
  
    size <- 100
          
                             x1 <- seq(from = xlim[1], to = xlim[2], length.out = size)         
                             
          if(mar2 %in% bin.link) x2 <- seq(from = ylim[1], to = ylim[2], length.out = size)                            
        if(!(mar2 %in% bin.link)) x2 <- seq(from = min(resp), to = max(resp), length.out = size)                        
          
          x11 <- rep(x1, each = size)
          x22 <- rep(x2, times = size)
          
          px11 <- probm(x11, margins[1], only.pr = FALSE)
          
          d.x1 <- px11$d.n 
          p1   <- px11$pr
          
          resf <- distrHsAT(x22, m2, s2, nu2, mar2)
          d.x2 <- resf$pdf2
          p2   <- resf$p2
          
          
          md <- copgHsAT(p1, p2, par1, BivD, Ln = TRUE)$c.copula2.be1be2*d.x1*d.x2    
          z  <- matrix(data = md, nrow = size, byrow = TRUE)
      
          filled.contour(x1, x2, z, color = topo.colors, nlevels = 16, ...) 
               
 }
 


  if(object$margins[2] %in% c(c1,c2)){
  
  if(object$margins[2] %in% c1) test.resp <- seq(from = min(object$y2), to = max(object$y2), length.out = n.grid)
  if(object$margins[2] %in% c2) test.resp <- seq(from = 0.0001, to = max(object$y2), length.out = n.grid)

  
  test.dens <- round( distrHsAT(test.resp, m2, s2, nu2, object$VC$margins[2])$pdf2, n.dig)
  resp <- test.resp[which(test.dens != 0)]
  
  if( length(resp) < 2 ) stop("Increase/decrease n.dig and/or n.grid values.")
  
 
  } else resp <- object$y2
  
  
  
 
Cplot(BivD = object$BivD, par1 = par1, mar2 = object$margins[2], 
      m2 = m2, s2 = s2, nu2 = nu2, resp = resp, main = cop, ylab = ylab, xlab = xlab, object$margins, ...)  
      
      
      
 }
 
 
 
rm(bs, SE, Vb, epds, sigma2.st, sigma2, est.RHOb, et1s, et2s, p1s, p2s, p11s, p10s, p00s, p01s, ORs, GMs, XX, Xt, V) 
 
  res <- list(tableP1=table[[1]], tableP2=table[[2]], tableP3=table[[3]], formula = object$formula,
              tableP4=table[[4]], tableP5=table[[5]], tableP6=table[[6]],
              tableNP1=tableN[[1]], tableNP2=tableN[[2]], tableNP3=tableN[[3]], 
              tableNP4=tableN[[4]], tableNP5=tableN[[5]], tableNP6=tableN[[6]], 
              n=n, theta.a=object$theta.a, sigma2.a=object$sigma2.a, nu.a=object$nu.a, 
              theta=object$theta, sigma2=object$sigma2, nu=object$nu, 
              formula1=object$gam1$formula, formula2=object$gam2$formula, formula3=object$gam3$formula,
              formula4=object$gam4$formula, formula5=object$gam5$formula, formula6=object$gam6$formula,
              t.edf=object$t.edf, CItheta=CIrs, CIsig2=CIsig2, CInu=CInu,  
              n.sel=n.sel, tau=object$tau, tau.a=object$tau.a, CIkt = CIkt, 
              BivD=object$BivD, margins = object$margins, bin.link = bin.link, 
              l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3, 
              l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6
              )
  class(res) <- "summary.copulaSampleSel"
      
                                        

res

}
