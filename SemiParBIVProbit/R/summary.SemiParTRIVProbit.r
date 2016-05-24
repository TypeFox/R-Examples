summary.SemiParTRIVProbit <- function(object, n.sim = 100, prob.lev = 0.05, ...){


  testStat <- getFromNamespace("testStat", "mgcv")
  liu2   <- getFromNamespace("liu2", "mgcv") 
  filled.contour <- getFromNamespace("filled.contour", "graphics")      

  bs <- SE <- Vb <- epds <- sigma2.st <- sigma2 <- nu.st <- nu <- est.RHOb <- et1s <- et2s <- p1s <- p2s <- p11s <- p10s <- p00s <- p01s <- ORs <- GMs <- XX <- Xt <- V <- 1
  
  n <- object$n
  
  tableN <- table <- list(NULL, NULL, NULL, NULL, NULL, NULL) 
  epsilon <- 0.0000001 
  max.p   <- 0.9999999
  
  lf <- length(object$coefficients)
  Vb <- object$Vb 
  SE <- sqrt(diag(Vb)) 

  bs <- rMVN(n.sim, mean = object$coefficients, sigma=Vb)  

  #epd12s <- t(as.matrix(   teta.tr(object$VC, bs[, lf-2])$teta  ) # t matrix needed when dealing with vectors
  #epd13s <- t(as.matrix(   teta.tr(object$VC, bs[, lf-1])$teta  )
  #epd23s <- t(as.matrix(   teta.tr(object$VC, bs[, lf])$teta    )
  
  epd12s <- teta.tr(object$VC, bs[, lf-2])$teta  
  epd13s <- teta.tr(object$VC, bs[, lf-1])$teta  
  epd23s <- teta.tr(object$VC, bs[, lf]  )$teta  
  
  CI12s <- rowQuantiles(t(epd12s), probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
  CI13s <- rowQuantiles(t(epd13s), probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
  CI23s <- rowQuantiles(t(epd23s), probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
   

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
  

 
rm(bs, SE, Vb, epds, sigma2.st, sigma2, est.RHOb, et1s, et2s, p1s, p2s, p11s, p10s, p00s, p01s, ORs, GMs, XX, Xt, V) 
 
  res <- list(tableP1=table[[1]], tableP2=table[[2]], tableP3=table[[3]], formula = object$formula,
              tableP4=table[[4]], tableP5=table[[5]], tableP6=table[[6]],
              tableNP1=tableN[[1]], tableNP2=tableN[[2]], tableNP3=tableN[[3]], 
              tableNP4=tableN[[4]], tableNP5=tableN[[5]], tableNP6=tableN[[6]], 
              n=n, 
              theta12.a=object$theta12.a, theta12=object$theta12, 
              theta13.a=object$theta13.a, theta13=object$theta13,
              theta23.a=object$theta23.a, theta23=object$theta23,
              formula1=object$gam1$formula, formula2=object$gam2$formula, formula3=object$gam3$formula,
              formula4=object$gam4$formula, formula5=object$gam5$formula, formula6=object$gam6$formula,
              t.edf=object$t.edf, CI12s=CI12s, CI13s=CI13s, CI23s=CI23s,
              l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3, 
              l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6
              )
  class(res) <- "summary.SemiParTRIVProbit"
      
                                        

res

}

