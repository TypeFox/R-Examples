summary.SemiParSampleSel <- function(object, n.sim=1000, s.meth="svd", prob.lev=0.05, cm.plot = FALSE, xlim = c(-3, 3), ylab = "Outcome margin", xlab = "Selection margin", ...){

  testStat <- getFromNamespace("testStat", "mgcv")
  reTest   <- getFromNamespace("reTest", "mgcv") 
  filled.contour <- getFromNamespace("filled.contour", "graphics")

  bd <- object$bd
  
  Vb <- object$Vb
          
  SE <- sqrt(diag(Vb))
  n  <- object$n 

  epsilon <- sqrt(.Machine$double.eps)

  bs <- rmvnorm(n.sim, mean = coef(object), sigma=Vb, method=s.meth)

  d.rho <- dim(Vb)[1]
  d.sig <- d.nu <- CIphi <- CInu <- NULL
  
  est.NUb <- est.SIGb <- est.THETAb <- NULL              ## est.KeTb  rep(NA,n.sim)  CIkt  <- as.numeric(quantile(est.KeTb,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE))


  if(object$margins[2] %in% c("G", "N", "NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2")) d.sig <- d.rho - 1 
  if(object$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG"))                            {            d.sig <- d.rho - 2; d.nu  <- d.rho - 1}
  
  if( is.null(object$X3) ) { 
     if(object$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE") )                                                                                                      etatt <- bs[,d.rho]    
     if(object$margins[2] %in% c("G", "N", "NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") ) etasg <- bs[,d.sig];                     etatt <- bs[,d.rho] 
     if(object$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG") )                                         etasg <- bs[,d.sig]; etanu <- bs[,d.nu]; etatt <- bs[,d.rho]          
  }  
  
  
  
  if( !is.null(object$X3) && object$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE")  ) {
     etatt <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)])
  }
  
  
  if( !is.null(object$X3) && object$margins[2] %in% c("G", "N", "NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2")  ) {
     etasg <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)])
     etatt <- object$X4%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)])
  } 
  
  
  if( !is.null(object$X3) && object$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG")  ) {
     etasg <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)])
     etanu <- object$X4%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)])       
     etatt <- object$X5%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2)])
  }   
  
  if(object$margins[2] %in% c("G", "N", "NB", "PIG", "D", "S", "BB", "NBII", "WARING", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG")) est.SIGb <- exp(etasg)
  if(object$margins[2] %in% c("ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") ) est.SIGb <- plogis(etasg)
  if(object$margins[2] %in% c("D", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG") ) est.NUb <- plogis(etanu)
  if(object$margins[2] == "S") est.NUb <- etanu 
  
        if(object$BivD=="N")                                  est.THETAb <- tanh(etatt)  
        if(object$BivD %in% c("C0", "C180"))                  est.THETAb <-   exp(etatt) + epsilon
        if(object$BivD %in% c("C90","C270"))                  est.THETAb <- -(exp(etatt) + epsilon)
        if(object$BivD %in% c("J0", "J180"))                  est.THETAb <- 1 + exp(etatt) + epsilon
        if(object$BivD %in% c("J90", "J270"))                 est.THETAb <- -(1 + exp(etatt) + epsilon)
  if(object$BivD=="FGM")                                est.THETAb <- tanh(etatt)
  if(object$BivD=="F")                                  est.THETAb <- etatt + epsilon
	if(object$BivD=="AMH")                                est.THETAb <- tanh(etatt)
        if(object$BivD %in% c("G0", "G180"))                  est.THETAb <- 1 + exp(etatt) 
        if(object$BivD %in% c("G90","G270"))                  est.THETAb <- -(1 + exp(etatt)) 


  if(object$margins[2] %in% c("G", "N", "NB", "PIG", "D", "S", "BB", "NBII", "WARING", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2")) 
                                                                                      { CIphi <- as.numeric(quantile(est.SIGb,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE)) }
  if(object$margins[2] %in% c("D", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG", "S") )    CInu <- as.numeric(quantile(est.NUb,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE)) 

                                                                                        CIth  <- as.numeric(quantile(est.THETAb,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE))
  
  tableN <- list(NULL, NULL, NULL, NULL, NULL)
  table  <- list(NULL, NULL, NULL, NULL, NULL)


  index <- 1:2
  ind1 <- 1:object$gp1
  ind2 <- object$X1.d2 + (1:object$gp2)
  ind3 <- ind4 <- ind5 <- ind6 <- NULL 
  
  ind1gam <- 1:object$X1.d2
  ind2gam <- object$X1.d2 + (1:object$X2.d2) 
  ind3gam <- ind4gam <- ind5gam <- ind6gam <- NULL   
  
  indsp1 <- 1
  
  
  if(!is.null(object$X3) ) {
  
       ind3 <- object$X1.d2 + object$X2.d2 + (1:object$gp3)
       index <- 1:3
       ind3gam <- object$X1.d2 + object$X2.d2 + (1:object$X3.d2)
       
       if(!is.null(object$X4) ) {
       
       ind4 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + (1:object$gp4)
       index <- 1:4
       ind4gam <- object$X1.d2 + object$X2.d2 + object$X3.d2 + (1:object$X4.d2)
       
                                }
                                
       if(!is.null(object$X5) ) {
       
       ind5 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + (1:object$gp5)
       index <- 1:5  
       ind5gam <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + (1:object$X5.d2)
                                }     
                                
       #if(!is.null(object$X6) ) {
       #
       #ind6 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + (1:object$gp6)
       #index <- 1:6  
       #ind6gam <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + (1:object$X6.d2)
       #                         }                                 
                            
  }
                            
  ind <- list( ind1 = ind1,
               ind2 = ind2,
               ind3 = ind3, 
               ind4 = ind4,
               ind5 = ind5,
               ind6 = ind6)
               
  indgam <- list( ind1gam = ind1gam,
                  ind2gam = ind2gam,
                  ind3gam = ind3gam, 
                  ind4gam = ind4gam,
                  ind5gam = ind5gam,
                  ind6gam = ind6gam)  
   

sp1 <- sp2 <- sp3 <- sp4 <- sp5 <- sp6 <- 0
indsp1 <- indsp2 <- indsp3 <- indsp4 <- indsp5 <- indsp6 <- 0

if(object$l.sp1!=0) {indsp1 <- 1:object$l.sp1;                              sp1 <- object$l.sp1  } 
if(object$l.sp2!=0) {indsp2 <- 1:object$l.sp2 + sp1;                        sp2 <- object$l.sp2  }
if(object$l.sp3!=0) {indsp3 <- 1:object$l.sp3 + sp1 + sp2;                  sp3 <- object$l.sp3  }
if(object$l.sp4!=0) {indsp4 <- 1:object$l.sp4 + sp1 + sp2 + sp3;            sp4 <- object$l.sp4  }
if(object$l.sp5!=0) indsp5 <- 1:object$l.sp5 + sp1 + sp2 + sp3 + sp4 #;      sp5 <- object$l.sp5  }
#if(object$l.sp6!=0) {indsp6 <- 1:object$l.sp6 + sp1 + sp2 + sp3 + sp4 + sp5                      }

indsp <- list(indsp1 = indsp1,
              indsp2 = indsp2,
              indsp3 = indsp3,
              indsp4 = indsp4,
              indsp5 = indsp5) #,
              # indsp6 = indsp6)


  for(i in index){
  estimate <- coef(object)[ind[[i]]]
  se       <- SE[ind[[i]]]
  ratio    <- estimate/se
  pv       <- 2*pnorm(abs(ratio), lower.tail = FALSE)
  table[[i]] <- cbind(estimate,se,ratio,pv)
  dimnames(table[[i]])[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }

  
  if( object$l.sp1!=0 || object$l.sp2!=0 || object$l.sp3!=0 || object$l.sp4!=0 || object$l.sp5!=0){

  	pTerms.df <- pTerms.chi.sq <- pTerms.pv <- tableN <- list(0,0,0,0,0)
        XX <- object$R
        
        
           for(i in index){

             if(i==1) {mm <- object$l.sp1; if(mm==0) next}
             if(i==2) {mm <- object$l.sp2; if(mm==0) next} 
             if(i==3) {mm <- object$l.sp3; if(mm==0) next} 
             if(i==4) {mm <- object$l.sp4; if(mm==0) next} 
             if(i==5) {mm <- object$l.sp5; if(mm==0) break} 
             #if(i==6) {mm <- object$l.sp6; if(mm==0) break} 
  
		for(k in 1:mm){

                        if(i==1){ gam <- object$gam1; ind <-  gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para                                } 
                        if(i==2){ gam <- object$gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2                } 
                        if(i==3){ gam <- object$gam3; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 }
                        if(i==4){ gam <- object$gam4; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 }
                        if(i==5){ gam <- object$gam5; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 }
                        #if(i==6){ gam <- object$gam6; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 }
                          
                        gam$sig2         <- 1
                        gam$scale.estimated <- FALSE                          
                          
                        if(gam$smooth[[k]]$null.space.dim == 0){
                        
                        gam$Ve           <- object$Ve[indgam[[i]], indgam[[i]], drop = FALSE] 
                        gam$sp           <- object$sp[ indsp[[i]] ]                       
                        gam$R            <- XX[indgam[[i]], indgam[[i]], drop = FALSE] # here we have a problem, need to enlarge gam object in R and S1 to have consistent output
                        gam$coefficients <- coef(object)[ indgam[[i]] ]
                        Tp <- reTest(gam, k)
                        
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


  }
  

    if (cm.plot == TRUE) {
    
       if (object$margins[2] %in% c("P", "NB", "D", "PIG", "S", "BB", "BI", "GEOM", "LG", 
            "NBII", "WARING", "YULE", "ZIBB", "ZABB", "ZABI", "ZIBI", 
            "ZALG", "ZANBI", "ZINBI", "ZAP", "ZIP", "ZIP2", "ZIPIG"))  {

	Mixed.pdf <- function(x1, x2, par1, BivD, outcome.margin, mu_o=mu_o, sigma=sigma, nu=nu, bd=bd) {
   		          F1 <- pnorm(x1)
                precision <- 10^(-7)
   			     if (outcome.margin=="P") {
        			  F2 <-  pPO(x2, mu=mu_o)
        			  F22 <-  pPO(x2, mu=mu_o) - dPO(x2, mu=mu_o)  
              	F2 <- ifelse(F2<(1-precision), F2, 1-precision)
              	F2 <- ifelse(F2>precision, F2, precision)    
              	F22 <- ifelse(F22<(1-precision), F22, 1-precision)
              	F22 <- ifelse(F22>precision, F22, precision)
            } else if (outcome.margin=="NB") {
        	  		F2  <-  pNBI(x2, mu=mu_o, sigma=sigma)
        		  	F22 <-  pNBI(x2, mu=mu_o, sigma=sigma) - dNBI(x2, mu=mu_o, sigma=sigma)  
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
            		F22 <- ifelse(F22>precision, F22, precision)    
     			  } else if (outcome.margin=="D") {
        			  F2   <- pDEL(x2, mu=mu_o, sigma=sigma, nu=nu)
        			  F22  <- pDEL(x2, mu=mu_o, sigma=sigma, nu=nu) - dDEL(x2, mu=mu_o, sigma=sigma, nu=nu)  
              	F2 <- ifelse(F2<(1-precision), F2, 1-precision)
              	F2 <- ifelse(F2>precision, F2, precision)    
         		    F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      			    F22 <- ifelse(F22>precision, F22, precision)   
     			  } else if (outcome.margin=="PIG") {
        			  F2  <-  pPIG(x2, mu=mu_o, sigma=sigma)
        			  F22 <-  pPIG(x2, mu=mu_o, sigma=sigma) - dNBI(x2, mu=mu_o, sigma=sigma)   
            		F2 <- ifelse(F2<(1-precision), F2, 1-precision)
            		F2 <- ifelse(F2>precision, F2, precision)   
            	 	F22 <- ifelse(F22<(1-precision), F22, 1-precision)
            		F22 <- ifelse(F22>precision, F22, precision)
            } else if (outcome.margin=="S") {
        			  F2  <-  pSICHEL(x2, mu=mu_o, sigma=sigma, nu=nu)
        			  F22  <- pSICHEL(x2, mu=mu_o, sigma=sigma, nu=nu) - dSICHEL(x2, mu=mu_o, sigma=sigma, nu=nu)   
              	F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
            		F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
     	      } else if (outcome.margin=="BB") {
            		F2  <-  pBB(x2, bd=bd, mu=mu_o, sigma=sigma)
          			F22  <- pBB(x2, bd=bd, mu=mu_o, sigma=sigma) - dBB(x2, bd=bd, mu=mu_o, sigma=sigma) 
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
            		F2 <- ifelse(F2>precision, F2, precision)    
              	F22 <- ifelse(F22<(1-precision), F22, 1-precision)
              	F22 <- ifelse(F22>precision, F22, precision)
            } else if (outcome.margin=="BI") {
            	  F2  <-  pBI(x2, bd=bd, mu=mu_o)
        			  F22  <- pBI(x2, bd=bd, mu=mu_o) - dBI(x2, bd=bd, mu=mu_o)  
           			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
            		F22 <- ifelse(F22>precision, F22, precision)
  	         } else if (outcome.margin=="GEOM") {
                F2  <-  pGEOM(x2, mu=mu_o)
        			  F22  <- pGEOM(x2, mu=mu_o) - dGEOM(x2, mu=mu_o)  
            		F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
     		     } else if (outcome.margin=="LG") {
                F2  <-  pLG(x2, mu=mu_o)
            		F22  <- pLG(x2, mu=mu_o) - dLG(x2, mu=mu_o)   
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
     		     } else if (outcome.margin=="NBII") {
            		F2  <-  pNBII(x2, mu=mu_o, sigma=sigma)
          			F22  <- pNBII(x2, mu=mu_o, sigma=sigma) - dNBII(x2, mu=mu_o, sigma=sigma)  
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
     		     } else if (outcome.margin=="WARING") {
              	F2  <-  pWARING(x2, mu=mu_o, sigma=sigma)
          			F22  <- pWARING(x2, mu=mu_o, sigma=sigma) - dWARING(x2, mu=mu_o, sigma=sigma)   
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
             } else if (outcome.margin=="YULE") {
                F2  <-  pYULE(x2, mu=mu_o)
        	  		F22  <- pYULE(x2, mu=mu_o) - dYULE(x2, mu=mu_o)  
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
    		     } else if (outcome.margin=="ZABB") {
            		F2  <-  pZABB(x2, bd=bd, mu=mu_o, sigma=sigma, nu=nu)
          			F22  <- pZABB(x2, bd=bd, mu=mu_o, sigma=sigma, nu=nu) - dZABB(x2, bd=bd, mu=mu_o, sigma=sigma, nu=nu)  
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
              	F22 <- ifelse(F22<(1-precision), F22, 1-precision)
              	F22 <- ifelse(F22>precision, F22, precision)
    		     } else if (outcome.margin=="ZIBB") {
              	F2  <-  pZIBB(x2, bd=bd, mu=mu_o, sigma=sigma, nu=nu)
          			F22  <- pZIBB(x2, bd=bd, mu=mu_o, sigma=sigma, nu=nu) - dZIBB(x2, bd=bd, mu=mu_o, sigma=sigma, nu=nu)  
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
             } else if (outcome.margin=="ZABI") {
                F2  <-  pZABI(x2, bd=bd, mu=mu_o, sigma=sigma)
        	  		F22  <- pZABI(x2, bd=bd, mu=mu_o, sigma=sigma) - dZABI(x2, bd=bd, mu=mu_o, sigma=sigma)  
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
             } else if (outcome.margin=="ZIBI") {
                F2  <-  pZIBI(x2, bd=bd, mu=mu_o, sigma=sigma)
            		F22  <- pZIBI(x2, bd=bd, mu=mu_o, sigma=sigma) - dZIBI(x2, bd=bd, mu=mu_o, sigma=sigma)  
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
     		     } else if (outcome.margin=="ZALG") {
                F2  <-  pZALG(x2, mu=mu_o, sigma=sigma)
              	F22  <- pZALG(x2, mu=mu_o, sigma=sigma) - dZALG(x2, mu=mu_o, sigma=sigma)   
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
            		F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
    		     } else if (outcome.margin=="ZANBI") {
            		F2  <-  pZANBI(x2, mu=mu_o, sigma=sigma, nu=nu)
          			F22  <- pZANBI(x2, mu=mu_o, sigma=sigma, nu=nu) - dZANBI(x2, mu=mu_o, sigma=sigma, nu=nu)   
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
              	F22 <- ifelse(F22<(1-precision), F22, 1-precision)
            		F22 <- ifelse(F22>precision, F22, precision)
     		     } else if (outcome.margin=="ZINBI") {
              	F2  <-  pZINBI(x2, mu=mu_o, sigma=sigma, nu=nu)
          			F22  <- pZINBI(x2, mu=mu_o, sigma=sigma, nu=nu) - dZINBI(x2, mu=mu_o, sigma=sigma, nu=nu)  
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
     		     } else if (outcome.margin=="ZAP") {
                F2  <-  pZAP(x2, mu=mu_o, sigma=sigma)
          			F22  <- pZAP(x2, mu=mu_o, sigma=sigma) - dZAP(x2, mu=mu_o, sigma=sigma)  
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
            		F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
              } else if (outcome.margin=="ZIP") {
                F2  <-  pZIP(x2, mu=mu_o, sigma=sigma)
          		  F22  <- pZIP(x2, mu=mu_o, sigma=sigma) - dZIP(x2, mu=mu_o, sigma=sigma)   
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
       		     } else if (outcome.margin=="ZIP2") {
                F2  <-  pZIP2(x2, mu=mu_o, sigma=sigma)
              	F22  <- pZIP2(x2, mu=mu_o, sigma=sigma) - dZIP2(x2, mu=mu_o, sigma=sigma)  
          			F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
       		     } else if (outcome.margin=="ZIPIG") {
                F2  <-  pZIPIG(x2, mu=mu_o, sigma=sigma, nu=nu)
                F22  <- pZIPIG(x2, mu=mu_o, sigma=sigma, nu=nu) - dZIPIG(x2, mu=mu_o, sigma=sigma, nu=nu)  
            		F2 <- ifelse(F2<(1-precision), F2, 1-precision)
          			F2 <- ifelse(F2>precision, F2, precision)    
          			F22 <- ifelse(F22<(1-precision), F22, 1-precision)
          			F22 <- ifelse(F22>precision, F22, precision)
                    }
  			if(BivD=="FGM") {
    				d.Cop.1<- F2*((1+par1*(1-F1)*(1-F2))-F1*par1*(1-F2))
    				d.Cop.2 <- F22*((1+par1*(1-F1)*(1-F22))-F1*par1*(1-F22))
   			} else  if (BivD=="N") {
    				d.Cop.1<- pnorm((qnorm(F2)-par1*qnorm(F1))/sqrt(1-par1^2))
    				d.Cop.2 <- pnorm((qnorm(F22)-par1*qnorm(F1))/sqrt(1-par1^2))
    			} else  if (BivD=="AMH") {
				    d.Cop.1<- (F2*(1-par1*(1-F1)*(1-F2))-F1*F2*(par1*(1-F2)))/(1-par1*(1-F1)*(1-F2))^2
    				d.Cop.2 <- (F22*(1-par1*(1-F1)*(1-F22))-F1*F22*(par1*(1-F22)))/(1-par1*(1-F1)*(1-F22))^2
  			} else  if (BivD=="C0") {
				    d.Cop.1<- F1^(-par1-1)*(F1^(-par1)+F2^(-par1)-1)^(-(1/par1)-1)
    				d.Cop.2 <- F1^(-par1-1)*(F1^(-par1)+F22^(-par1)-1)^(-(1/par1)-1)
  			} else  if (BivD=="C90") {
    				par1 <- -par1
    				d.Cop.1<- (1-F1)^(-par1-1)*((1-F1)^(-par1)+F2^(-par1)-1)^(-(1/par1)-1)
    				d.Cop.2 <- (1-F1)^(-par1-1)*((1-F1)^(-par1)+F22^(-par1)-1)^(-(1/par1)-1)
  			} else  if (BivD=="C180") {
    				d.Cop.1<- 1-(1-F1)^(-par1-1)*((1-F1)^(-par1)+(1-F2)^(-par1)-1)^(-(1/par1)-1)
    				d.Cop.2 <-1-(1-F1)^(-par1-1)*((1-F1)^(-par1)+(1-F22)^(-par1)-1)^(-(1/par1)-1)
  			} else  if (BivD=="C270") {
    				par1 <- -par1
    				d.Cop.1<- 1-F1^(-par1-1)*(F1^(-par1)+(1-F2)^(-par1)-1)^(-(1/par1)-1)
    				d.Cop.2 <- 1-F1^(-par1-1)*(F1^(-par1)+(1-F22)^(-par1)-1)^(-(1/par1)-1)
  			} else  if (BivD=="F") {
 				    d.Cop.1<- -(1/par1)*(1/(1+(exp(-par1*F1)-1)*(exp(-par1*F2)-1)/(exp(-par1)-1)))*((exp(-par1*F2)-1)/(exp(-par1)-1))*exp(-par1*F1)*(-par1)
    				d.Cop.2 <- -(1/par1)*(1/(1+(exp(-par1*F1)-1)*(exp(-par1*F22)-1)/(exp(-par1)-1)))*((exp(-par1*F22)-1)/(exp(-par1)-1))*exp(-par1*F1)*(-par1) 
  			} else  if (BivD=="G0") {
				    d.Cop.1<- exp(-((-log(F1))^(par1)+(-log(F2))^(par1))^(1/par1))*((-1/par1)*((-log(F1))^(par1)+(-log(F2))^(par1))^(1/par1-1))*(par1*(-log(F1))^(par1-1))*(-1/F1) 
    				d.Cop.2 <- exp(-((-log(F1))^(par1)+(-log(F22))^(par1))^(1/par1))*((-1/par1)*((-log(F1))^(par1)+(-log(F22))^(par1))^(1/par1-1))*(par1*(-log(F1))^(par1-1))*(-1/F1)  
  			} else  if (BivD=="G90") {
				    par1 <- -par1
    				d.Cop.1<- exp(-((-log((1-F1)))^(par1)+(-log(F2))^(par1))^(1/par1))*((-1/par1)*((-log((1-F1)))^(par1)+(-log(F2))^(par1))^(1/par1-1))*(par1*(-log((1-F1)))^(par1-1))*(-1/(1-F1)) 
    				d.Cop.2 <- exp(-((-log((1-F1)))^(par1)+(-log(F22))^(par1))^(1/par1))*((-1/par1)*((-log((1-F1)))^(par1)+(-log(F22))^(par1))^(1/par1-1))*(par1*(-log((1-F1)))^(par1-1))*(-1/(1-F1)) 
  			} else  if (BivD=="G180") {
				    d.Cop.1<- 1-exp(-((-log((1-F1)))^(par1)+(-log((1-F2)))^(par1))^(1/par1))*((-1/par1)*((-log((1-F1)))^(par1)+(-log((1-F2)))^(par1))^(1/par1-1))*(par1*(-log((1-F1)))^(par1-1))*(-1/(1-F1)) 
    				d.Cop.2 <- 1-exp(-((-log((1-F1)))^(par1)+(-log((1-F22)))^(par1))^(1/par1))*((-1/par1)*((-log((1-F1)))^(par1)+(-log((1-F22)))^(par1))^(1/par1-1))*(par1*(-log((1-F1)))^(par1-1))*(-1/(1-F1))  
  			} else  if (BivD=="G270") {
				    par1 <- -par1
    				d.Cop.1<- 1-exp(-((-log(F1))^(par1)+(-log((1-F2)))^(par1))^(1/par1))*((-1/par1)*((-log(F1))^(par1)+(-log((1-F2)))^(par1))^(1/par1-1))*(par1*(-log(F1))^(par1-1))*(-1/F1) 
    				d.Cop.2 <- 1-exp(-((-log(F1))^(par1)+(-log((1-F22)))^(par1))^(1/par1))*((-1/par1)*((-log(F1))^(par1)+(-log((1-F22)))^(par1))^(1/par1-1))*(par1*(-log(F1))^(par1-1))*(-1/F1)  
  			} else  if (BivD=="J0") {
   				d.Cop.1<- (-1/par1)*((1-F1)^par1+(1-F2)^par1-(1-F1)^par1*(1-F2)^par1)^(1/par1-1)*(-par1*(1-F1)^(par1-1)+(1-F2)^par1*par1*(1-F1)^(par1-1))
    				d.Cop.2<- (-1/par1)*((1-F1)^par1+(1-F22)^par1-(1-F1)^par1*(1-F22)^par1)^(1/par1-1)*(-par1*(1-F1)^(par1-1)+(1-F22)^par1*par1*(1-F1)^(par1-1))   
    			} else  if (BivD=="J90") {
				    par1 <- -par1
    				d.Cop.1<- (-1/par1)*((1-(1-F1))^par1+(1-F2)^par1-(1-(1-F1))^par1*(1-F2)^par1)^(1/par1-1)*(-par1*(1-(1-F1))^(par1-1)+(1-F2)^par1*par1*(1-(1-F1))^(par1-1))
    				d.Cop.2<- (-1/par1)*((1-(1-F1))^par1+(1-F22)^par1-(1-(1-F1))^par1*(1-F22)^par1)^(1/par1-1)*(-par1*(1-(1-F1))^(par1-1)+(1-F22)^par1*par1*(1-(1-F1))^(par1-1))   
  			} else  if (BivD=="J180") {
  				d.Cop.1<- 1-(-1/par1)*((1-(1-F1))^par1+(1-(1-F2))^par1-(1-(1-F1))^par1*(1-(1-F2))^par1)^(1/par1-1)*(-par1*(1-(1-F1))^(par1-1)+(1-(1-F2))^par1*par1*(1-(1-F1))^(par1-1))
    				d.Cop.2 <- 1-(-1/par1)*((1-(1-F1))^par1+(1-(1-F22))^par1-(1-(1-F1))^par1*(1-(1-F22))^par1)^(1/par1-1)*(-par1*(1-(1-F1))^(par1-1)+(1-(1-F22))^par1*par1*(1-(1-F1))^(par1-1))   
  			} else if (BivD=="J270") {
				    par1 <- -par1
    				d.Cop.1<- 1-(-1/par1)*((1-F1)^par1+(1-(1-F2))^par1-(1-F1)^par1*(1-(1-F2))^par1)^(1/par1-1)*(-par1*(1-F1)^(par1-1)+(1-(1-F2))^par1*par1*(1-F1)^(par1-1))
    				d.Cop.2 <- 1-(-1/par1)*((1-F1)^par1+(1-(1-F22))^par1-(1-F1)^par1*(1-(1-F22))^par1)^(1/par1-1)*(-par1*(1-F1)^(par1-1)+(1-(1-F22))^par1*par1*(1-F1)^(par1-1))   
  			} 
 
 			res <- d.Cop.1*dnorm(x1) - d.Cop.2*dnorm(x1)
 			res
 
		}
         theta <- mean(object$theta)
        if (object$margins[2]  %in% c("P", "GEOM", "YULE")) {
              mu_o <- mean(exp(object$X2s%*%object$coefficients[(object$X1.d2 + 1):(object$X1.d2 + object$X2.d2)]))
        } else if (object$margins[2] %in% c("BI", "LG")) {
              mu_o <- mean(plogis(object$X2s%*%object$coefficients[(object$X1.d2 + 1):(object$X1.d2 + object$X2.d2)]))
        } else if (object$margins[2] %in% c("NB", "PIG", "NBII", "BB", "WARING", "ZAP", "ZIP", "ZIP2")) {
              mu_o <- mean(exp(object$X2s%*%object$coefficients[(object$X1.d2 + 1):(object$X1.d2 + object$X2.d2)]))
              sigma <- mean(object$sigma)
        } else if (object$margins[2] %in% c("ZABI", "ZIBI", "ZALG")) {
              mu_o <- mean(plogis(object$X2s%*%object$coefficients[(object$X1.d2 + 1):(object$X1.d2 + object$X2.d2)]))
              sigma <- mean(object$sigma)
        } else if (object$margins[2] %in% c("D", "S", "ZANBI", "ZINBI", "ZIPIG")) {
              mu_o <- mean(exp(object$X2s%*%object$coefficients[(object$X1.d2 + 1):(object$X1.d2 + object$X2.d2)]))
              sigma <- mean(object$sigma)
              nu <- mean(object$nu)
        } else if (object$margins[2] %in% c("ZABB", "ZIBB")) {
              mu_o <- mean(plogis(object$X2s%*%object$coefficients[(object$X1.d2 + 1):(object$X1.d2 + object$X2.d2)]))
              sigma <- mean(object$sigma)
              nu <- mean(object$nu)
        } 

  		  Cplot1 <- function (BivD, par1, xlim, resp, outcome.margin="P", mu_o=1, sigma=1, nu=0.5, bd=1, ...){    
          		  size <- 100
            		x  <- seq(from = xlim[1], to = xlim[2], length.out = size)
                y  <- min(resp):max(resp)
            		x1 <- rep(x = x, each = (max(resp)-min(resp)+1))
            		y1  <- y - mean(y)
          			x2 <- rep(x = y, times = size)
          			md <- Mixed.pdf(x1=x1, x2=x2, par1=par1, BivD=BivD, outcome.margin=outcome.margin, mu_o=mu_o, sigma=sigma, nu=nu, bd=bd)
          			z  <- matrix(data = md, nrow = size, byrow = T)
          			filled.contour <- getFromNamespace("filled.contour", "graphics")
          			filled.contour(x=x, y=y1, z, color = heat.colors, nlevels = 16, ...)
          			
  }

	if (object$BivD == "N") {
            cop <- bquote(paste("Gaussian (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "F") {
            cop <- bquote(paste("Frank (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "FGM") {
            cop <- bquote(paste("FGM (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "AMH") {
            cop <- bquote(paste("AMH (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "C0") {
            cop <- bquote(paste("Clayton (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "C90") {
            cop <- bquote(paste("Rotated Clayton - 90 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "C180") {
            cop <- bquote(paste("Survival Clayton (", hat(theta), 
                " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "C270") {
            cop <- bquote(paste("Rotated Clayton - 270 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "J0") {
            cop <- bquote(paste("Joe (", hat(theta), " = ", .(round(object$theta, 
                2)), ")", sep = ""))
        }
        if (object$BivD == "J90") {
            cop <- bquote(paste("Rotated Joe - 90 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "J180") {
            cop <- bquote(paste("Survival Joe (", hat(theta), 
                " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "J270") {
            cop <- bquote(paste("Rotated Joe - 270 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "G0") {
            cop <- bquote(paste("Gumbel (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "G90") {
            cop <- bquote(paste("Rotated Gumbel - 90 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "G180") {
            cop <- bquote(paste("Survival Gumbel (", hat(theta), 
                " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "G270") {
            cop <- bquote(paste("Rotated Gumbel - 270 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
               
         
  		test.resp  <- 0:max(object$y2)
      if (object$margins[2]=="LG") { test.resp  <- 1:max(object$y2) }
      if (object$margins[2]=="P") {
        f2 <- dPO(test.resp, mu=mu_o)
     } else if (object$margins[2]=="NB") {
        f2 <- dNBI(test.resp, mu=mu_o, sigma=sigma)
     } else if (object$margins[2]=="D") {
        f2 <- dDEL(test.resp, mu=mu_o, sigma=sigma, nu=nu)
     } else if (object$margins[2]=="PIG") {
        f2 <- dPIG(test.resp, mu=mu_o, sigma=sigma)
     } else if (object$margins[2]=="S") {
        f2 <- dSICHEL(test.resp, mu=mu_o, sigma=sigma, nu=nu)
     }  else if (object$margins[2]=="BB") {
        f2 <- dBB(test.resp, bd=bd, mu=mu_o, sigma=sigma)
     } else if (object$margins[2]=="BI") {
        f2 <- dBI(test.resp, bd=bd, mu=mu_o)
     } else if (object$margins[2]=="GEOM") {
        f2 <- dGEOM(test.resp, mu=mu_o)
     } else if (object$margins[2]=="LG") {
        f2 <- dLG(test.resp, mu=mu_o)
     } else if (object$margins[2]=="NBII") {
        f2 <- dNBII(test.resp, mu=mu_o, sigma=sigma)
     } else if (object$margins[2]=="WARING") {
        f2 <- dWARING(test.resp, mu=mu_o, sigma=sigma)
     } else if (object$margins[2]=="YULE") {
        f2 <- dYULE(test.resp, mu=mu_o)
     } else if (object$margins[2]=="ZABB") {
        f2 <- dZABB(test.resp, bd=bd, mu=mu_o, sigma=sigma, nu=nu)
     } else if (object$margins[2]=="ZIBB") {
        f2 <- dZIBB(test.resp, bd=bd, mu=mu_o, sigma=sigma, nu=nu)
     } else if (object$margins[2]=="ZABI") {
        f2 <- dZABI(test.resp, bd=bd, mu=mu_o, sigma=sigma)
     } else if (object$margins[2]=="ZIBI") {
        f2 <- dZIBI(test.resp, bd=bd, mu=mu_o, sigma=sigma)
     } else if (object$margins[2]=="ZALG") {
        f2 <- dZALG(test.resp, mu=mu_o, sigma=sigma)
     } else if (object$margins[2]=="ZANBI") {
        f2 <- dZANBI(test.resp, mu=mu_o, sigma=sigma, nu=nu)
     } else if (object$margins[2]=="ZINBI") {
        f2 <- dZINBI(test.resp, mu=mu_o, sigma=sigma, nu=nu)
     } else if (object$margins[2]=="ZAP") {
        f2 <- dZAP(test.resp, mu=mu_o, sigma=sigma)
     } else if (object$margins[2]=="ZIP") {
        f2 <- dZIP(test.resp, mu=mu_o, sigma=sigma)
     } else if (object$margins[2]=="ZIP2") {
        f2 <- dZIP2(test.resp, mu=mu_o, sigma=sigma)
     } else if (object$margins[2]=="ZIPIG") {
        f2 <- dZIPIG(test.resp, mu=mu_o, sigma=sigma, nu=nu)
     }   
      precision <- 10^-8
      f2 <- ifelse(f2>precision, f2, precision)  
      test.dens <- round(f2, 2)
      resp <- test.resp[which(test.dens != 0)]
      
  		  
      Cplot1(BivD=object$BivD, par1=theta, xlim=xlim, resp=resp, main=cop, xlab=xlab, ylab=ylab, outcome.margin=object$margins[2], mu_o=mu_o, sigma=sigma, nu=nu, bd=bd, ...)

       
       } else if (object$margins[2] %in% c("N", "G")) {
       

		    Cont.pdf <- function(x1, x2, par1, BivD, outcome.margin, mu.o, sigma.o, shape, rate = 1) {
        
        precision <- 10^(-7)
 				F1 <- pnorm(x1)
        f1 <- dnorm(x1)
      
  				if (outcome.margin=="N") {
        				F2 <-  pnorm(x2, mean=mu.o, sd=sigma.o)
        				f2 <-  dnorm(x2, mean=mu.o, sd=sigma.o)
                F2 <- ifelse(F2<(1-precision), F2, 1-precision)
              	F2 <- ifelse(F2>precision, F2, precision)    
              	f2 <- ifelse(f2>precision, f2, precision)
     				} else if (outcome.margin=="G") {
        				F2  <-  pgamma(x2, shape=shape, rate=rate)
        				f2  <-  dgamma(x2, shape=shape, rate=rate) 
                F2 <- ifelse(F2<(1-precision), F2, 1-precision)
                F2 <- ifelse(F2>precision, F2, precision)    
              	f2 <- ifelse(f2>precision, f2, precision)
    				} 

      
  				if(BivD=="FGM") {
    				d2Cdcop1dcop2 <- (1+par1*(1-F1)*(1-F2)-F1*par1*(1-F2))+F2*(-par1*(1-F1)+F1*par1)	
  				} else  if (BivD=="N") {
    				d2Cdcop1dcop2 <- dnorm((qnorm(F2)-par1*qnorm(F1))/sqrt(1-par1^2))*(dnorm(qnorm(F2)))^(-1)*(1-par1^2)^(-1/2)	
  				} else  if (BivD=="AMH") {
   					d2Cdcop1dcop2 <- (((1-par1*(1-F1)*(1-F2))+F2*par1*(1-F1)-(F1*par1-2*par1*F1*F2))*(1-par1*(1-F1)*(1-F2))^2-(F2*(1-par1*(1-F1)*(1-F2))-F1*F2*(par1*(1-F2)))*(1-par1*(1-F1)*(1-F2))*2*par1*(1-F1))/(1-par1*(1-F1)*(1-F2))^4
 				  } else  if (BivD=="C0") {
    				d2Cdcop1dcop2 <- F1^(-par1-1)*(F1^(-par1)+F2^(-par1)-1)^(-(1/par1)-2)*(-1/par1-1)*(-par1)*F2^(-par1-1)	
 				  } else  if (BivD=="C90") {
    				par1 <- -par1
    				d2Cdcop1dcop2 <- (1-F1)^(-par1-1)*((1-F1)^(-par1)+F2^(-par1)-1)^(-(1/par1)-2)*(-1/par1-1)*(-par1)*F2^(-par1-1)
  				} else  if (BivD=="C180") {
    				d2Cdcop1dcop2 <-  ( (1-F1)^(-par1-1)*((1-F1)^(-par1)+(1-F2)^(-par1)-1)^(-(1/par1)-2)*(-1/par1-1)*(-par1)*(1-F2)^(-par1-1)  )	
  				} else  if (BivD=="C270") {
    				par1 <- -par1
    				d2Cdcop1dcop2 <- F1^(-par1-1)*(F1^(-par1)+(1-F2)^(-par1)-1)^(-(1/par1)-2)*(-1/par1-1)*(-par1)*(1-F2)^(-par1-1)	
  				} else  if (BivD=="F") {
    				d2Cdcop1dcop2 <- (exp(-par1*F1)/(exp(-par1)-1)*exp(-par1*F2)*(-par1)*(1+(exp(-par1*F1)-1)*(exp(-par1*F2)-1)/(exp(-par1)-1))-(exp(-par1*F2)-1)/(exp(-par1)-1)*exp(-par1*F1)*(exp(-par1*F1)-1)/(exp(-par1)-1)*exp(-par1*F2)*(-par1))/(1+(exp(-par1*F1)-1)*(exp(-par1*F2)-1)/(exp(-par1)-1))^2	
  				} else  if (BivD=="G0") {
            dcop2 <- exp(-((-log(F1))^(par1)+(-log(F2))^(par1))^(1/par1))*((-1/par1)*((-log(F1))^(par1)+(-log(F2))^(par1))^(1/par1-1))*(par1*(-log(F2))^(par1-1))*(-1/F2)
    				d2Cdcop1dcop2 <- dcop2*((-1/par1)*((-log(F1))^(par1)+(-log(F2))^(par1))^(1/par1-1))*(par1*(-log(F1))^(par1-1))*(-1/F1)+exp(-((-log(F1))^(par1)+(-log(F2))^(par1))^(1/par1))*( (-1/par1)*(1/par1-1)*((-log(F1))^(par1)+(-log(F2))^(par1))^(1/par1-2)*par1*(-log(F2))^(par1-1)*(-1/F2))*par1*(-log(F1))^(par1-1)*(-1/F1)
  				} else  if (BivD=="G90") {
    				par1 <- -par1
            dcop2 <- 1-exp(-((-log((1-F1)))^(par1)+(-log(F2))^(par1))^(1/par1))*((-1/par1)*((-log((1-F1)))^(par1)+(-log(F2))^(par1))^(1/par1-1))*(par1*(-log(F2))^(par1-1))*(-1/F2)
    				d2Cdcop1dcop2 <- exp(-((-log((1-F1)))^(par1)+(-log(F2))^(par1))^(1/par1))*((-1/par1)*((-log((1-F1)))^(par1)+(-log(F2))^(par1))^(1/par1-1))*(par1*(-log(F2))^(par1-1))*(-1/F2)*((-1/par1)*((-log((1-F1)))^(par1)+(-log(F2))^(par1))^(1/par1-1))*(par1*(-log((1-F1)))^(par1-1))*(-1/(1-F1))+exp(-((-log((1-F1)))^(par1)+(-log(F2))^(par1))^(1/par1))*( (-1/par1)*(1/par1-1)*((-log((1-F1)))^(par1)+(-log(F2))^(par1))^(1/par1-2)*par1*(-log(F2))^(par1-1)*(-1/F2))*par1*(-log((1-F1)))^(par1-1)*(-1/(1-F1))	 	
  				} else  if (BivD=="G180") {
      			par1 <- -par1
            dcop2 <- 1-exp(-((-log((1-F1)))^(par1)+(-log((1-F2)))^(par1))^(1/par1))*((-1/par1)*((-log((1-F1)))^(par1)+(-log((1-F2)))^(par1))^(1/par1-1))*(par1*(-log((1-F2)))^(par1-1))*(-1/(1-F2))
    				d2Cdcop1dcop2 <- ( exp(-((-log((1-F1)))^(par1)+(-log((1-F2)))^(par1))^(1/par1))*((-1/par1)*((-log((1-F1)))^(par1)+(-log((1-F2)))^(par1))^(1/par1-1))*(par1*(-log((1-F2)))^(par1-1))*(-1/(1-F2))*((-1/par1)*((-log((1-F1)))^(par1)+(-log((1-F2)))^(par1))^(1/par1-1))*(par1*(-log((1-F1)))^(par1-1))*(-1/(1-F1))+exp(-((-log((1-F1)))^(par1)+(-log((1-F2)))^(par1))^(1/par1))*( (-1/par1)*(1/par1-1)*((-log((1-F1)))^(par1)+(-log((1-F2)))^(par1))^(1/par1-2)*par1*(-log((1-F2)))^(par1-1)*(-1/(1-F2)))*par1*(-log((1-F1)))^(par1-1)*(-1/(1-F1))     )	
          } else  if (BivD=="G270") {
    				par1 <- -par1
            dcop2 <- exp(-((-log(F1))^(par1)+(-log((1-F2)))^(par1))^(1/par1))*((-1/par1)*((-log(F1))^(par1)+(-log((1-F2)))^(par1))^(1/par1-1))*(par1*(-log((1-F2)))^(par1-1))*(-1/(1-F2))
    				d2Cdcop1dcop2 <- dcop2*((-1/par1)*((-log(F1))^(par1)+(-log((1-F2)))^(par1))^(1/par1-1))*(par1*(-log(F1))^(par1-1))*(-1/F1)+exp(-((-log(F1))^(par1)+(-log((1-F2)))^(par1))^(1/par1))*( (-1/par1)*(1/par1-1)*((-log(F1))^(par1)+(-log((1-F2)))^(par1))^(1/par1-2)*par1*(-log((1-F2)))^(par1-1)*(-1/(1-F2)))*par1*(-log(F1))^(par1-1)*(-1/F1)	
  				} else  if (BivD=="J0") {
    				d2Cdcop1dcop2 <- -1/par1*(1/par1-1)*((1-F1)^(par1)+(1-F2)^(par1)-(1-F1)^(par1)*(1-F2)^(par1))^(1/par1-2)*(-par1*(1-F2)^(par1-1)+(1-F1)^(par1)*par1*(1-F2)^(par1-1))*(-par1*(1-F1)^(par1-1)+(1-F2)^(par1)*par1*(1-F1)^(par1-1))-1/par1*((1-F1)^(par1)+(1-F2)^(par1)-(1-F1)^(par1)*(1-F2)^(par1))^(1/par1-1)*(-par1^2*(1-F2)^(par1-1)*(1-F1)^(par1-1))	
  				} else  if (BivD=="J90") {
    				par1 <- -par1
    				d2Cdcop1dcop2 <- -1/par1*(1/par1-1)*((1-(1-F1))^(par1)+(1-F2)^(par1)-(1-(1-F1))^(par1)*(1-F2)^(par1))^(1/par1-2)*(-par1*(1-F2)^(par1-1)+(1-(1-F1))^(par1)*par1*(1-F2)^(par1-1))*(-par1*(1-(1-F1))^(par1-1)+(1-F2)^(par1)*par1*(1-(1-F1))^(par1-1))-1/par1*((1-(1-F1))^(par1)+(1-F2)^(par1)-(1-(1-F1))^(par1)*(1-F2)^(par1))^(1/par1-1)*(-par1^2*(1-F2)^(par1-1)*(1-(1-F1))^(par1-1))	
  				} else  if (BivD=="J180") {
    				d2Cdcop1dcop2 <- ( -1/par1*(1/par1-1)*((1-(1-F1))^(par1)+(1-(1-F2))^(par1)-(1-(1-F1))^(par1)*(1-(1-F2))^(par1))^(1/par1-2)*(-par1*(1-(1-F2))^(par1-1)+(1-(1-F1))^(par1)*par1*(1-(1-F2))^(par1-1))*(-par1*(1-(1-F1))^(par1-1)+(1-(1-F2))^(par1)*par1*(1-(1-F1))^(par1-1))-1/par1*((1-(1-F1))^(par1)+(1-(1-F2))^(par1)-(1-(1-F1))^(par1)*(1-(1-F2))^(par1))^(1/par1-1)*(-par1^2*(1-(1-F2))^(par1-1)*(1-(1-F1))^(par1-1))   )	
  				} else if (BivD=="J270") {
    				par1 <- -par1
            d2Cdcop1dcop2 <- -1/par1*(1/par1-1)*((1-F1)^(par1)+(1-(1-F2))^(par1)-(1-F1)^(par1)*(1-(1-F2))^(par1))^(1/par1-2)*(-par1*(1-(1-F2))^(par1-1)+(1-F1)^(par1)*par1*(1-(1-F2))^(par1-1))*(-par1*(1-F1)^(par1-1)+(1-(1-F2))^(par1)*par1*(1-F1)^(par1-1))-1/par1*((1-F1)^(par1)+(1-(1-F2))^(par1)-(1-F1)^(par1)*(1-(1-F2))^(par1))^(1/par1-1)*(-par1^2*(1-(1-F2))^(par1-1)*(1-F1)^(par1-1))
    					
  } 
 
       
      
 			res <-   d2Cdcop1dcop2*f1*f2
   
 			res
 
		}
 
 
          theta <- mean(object$theta)
          	   if (object$margins[2]=="N") {
                      mu.o <- mean(object$X2s%*%object$coefficients[(object$X1.d2 + 1):(object$X1.d2 + object$X2.d2)])
                      sigma <- mean(object$sigma)
                       shape <- NULL
                      rate <- NULL
                   } else if (object$margins[2]=="G") {
                      mu.o <- NULL
                      # sigma <- NULL
                      sigma <- shape <- mean(object$sigma)
                      rate <- mean(1/object$phi)*mean(exp(-object$eta2))
                   } 


  	Cplot2 <- function (resp, BivD, par1, xlim, outcome.margin, mu.o, sigma.o, shape, rate, ...){
              size <- 100    
          		x  <- seq(from = xlim[1], to = xlim[2], length.out = size)
              y <- seq(from = min(resp), to = max(resp), length.out = size)
              y1 <- y - mean(y)
          	  x1 <- rep(x = x, each = size)
          		x2 <- rep(x = y, times = size)
          		md <- Cont.pdf(x1=x1, x2=x2, par1=par1, BivD=BivD, outcome.margin=outcome.margin, mu.o=mu.o, sigma.o=sigma.o, shape=shape, rate = rate)
          		z  <- matrix(data = md, nrow = size, byrow = TRUE)
          		filled.contour(x=x, y=y1, z, color = heat.colors, nlevels = 16, ...)

  		}
	if (object$BivD == "N") {
            cop <- bquote(paste("Gaussian (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "F") {
            cop <- bquote(paste("Frank (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "FGM") {
            cop <- bquote(paste("FGM (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "AMH") {
            cop <- bquote(paste("AMH (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "C0") {
            cop <- bquote(paste("Clayton (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "C90") {
            cop <- bquote(paste("Rotated Clayton - 90 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "C180") {
            cop <- bquote(paste("Survival Clayton (", hat(theta), 
                " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "C270") {
            cop <- bquote(paste("Rotated Clayton - 270 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "J0") {
            cop <- bquote(paste("Joe (", hat(theta), " = ", .(round(object$theta, 
                2)), ")", sep = ""))
        }
        if (object$BivD == "J90") {
            cop <- bquote(paste("Rotated Joe - 90 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "J180") {
            cop <- bquote(paste("Survival Joe (", hat(theta), 
                " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "J270") {
            cop <- bquote(paste("Rotated Joe - 270 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "G0") {
            cop <- bquote(paste("Gumbel (", hat(theta), " = ", 
                .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "G90") {
            cop <- bquote(paste("Rotated Gumbel - 90 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "G180") {
            cop <- bquote(paste("Survival Gumbel (", hat(theta), 
                " = ", .(round(object$theta, 2)), ")", sep = ""))
        }
        if (object$BivD == "G270") {
            cop <- bquote(paste("Rotated Gumbel - 270 degrees (", 
                hat(theta), " = ", .(round(object$theta, 2)), ")", sep = ""))
        }

  	   if (object$margins[2] %in% c("G")) {
          test.resp <- seq(from = 1e-04, to = max(object$y2), length.out = 1000)
          f2 <- dgamma(test.resp, shape = shape, rate = rate)
	      }
       if (object$margins[2] %in% c("N")) {
          test.resp <- seq(from = min(object$y2), to = max(object$y2), length.out = 1000)
	        f2 <- dnorm(test.resp, mean = mu.o, sd = sigma)
        }
  	      test.dens <- round(f2, 2)
          resp <- test.resp[which(test.dens != 0)]
  	
       Cplot2(resp=resp, BivD=object$BivD, par1=theta, xlim=xlim, xlab=xlab, ylab=ylab, main=cop, outcome.margin=object$margins[2], mu.o=mu.o, sigma=sigma, shape=shape, rate=rate, ...)


      }
    }


  
  

     res <- list(tableP1=table[[1]], tableP2=table[[2]], tableP3=table[[3]], tableP4=table[[4]], tableP5=table[[5]], 
                 tableNP1=tableN[[1]], tableNP2=tableN[[2]], tableNP3=tableN[[3]], tableNP4=tableN[[4]], tableNP5=tableN[[5]], 
                 n=n, phi=object$phi.a, sigma=object$sigma.a, theta=object$theta.a,  
                 nu = object$nu.a, tau= NULL,  # object$tau, 
                 formula1=object$gam1$formula, formula2=object$gam2$formula,
                 formula3=object$gam3$formula, formula4=object$gam4$formula, 
                 formula5=object$gam5$formula, formula = object$formula,             
                 l.sp1=object$l.sp1, l.sp2=object$l.sp2, l.sp3=object$l.sp3, l.sp4=object$l.sp4, l.sp5=object$l.sp5, 
                 t.edf=object$t.edf, CIsig=CIphi, CIth=CIth, CInu=CInu, CIkt= NULL, #CIkt, 
                 BivD=object$BivD, margins=object$margins, n.sel=object$n.sel)
  

  class(res) <- "summary.SemiParSampleSel"

  res


}