
##############################################
# Q3 and model fit statistics for objects of class tam
tam.modelfit <- function( tamobj , progress=TRUE ){
  # extract data from tam object
  # zz0 <- Sys.time()
  mod1 <- tamobj
  resp0 <- resp <- mod1$resp
  jkunits <- 1
  maxKi <- apply( resp , 2 , max , na.rm=TRUE)
  resp.ind <- mod1$resp.ind
  rprobs <- mod1$rprobs
  theta <- mod1$theta
  hwt <-  mod1$hwt
  maxK <- dim(rprobs)[2]
  TP <- dim(rprobs)[3]
  residM <- matrix( 0 , nrow=nrow(resp) , ncol=ncol(resp) )
  N <- nrow(resp)
  I <- ncol(resp)
  # calculate Q3 and residuals
  #*** Rcpp call
  if (progress){
    cat("**** Calculate Residuals \n") ; utils::flush.console()
    #		cat(paste0("     |",paste0(rep("*", 10),collapse="") , "|\n     |"))
  }
  RR <- I*(I-1) / 2 
  res0 <- .Call("tam_q3_calc_residM" , as.vector( rprobs ) , as.matrix(resp) , 
                I , TP , maxK ,  maxKi , hwt , PACKAGE="TAM" )								
  residM <- res0$residM
  resp[ resp.ind == 0 ] <- NA
  residM <- resp - residM
  # cat("calc residM Rcpp") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1		
  res0 <- .Call("tam_q3_calc_V2q3jack", residM , as.matrix(resp.ind)  ,
                PACKAGE="TAM")
  # cat("calc tam_q3_calc_V2q3jack") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1		
  dfr <- as.data.frame( res0$dfr )
  colnames(dfr) <- c("index1" , "index2" , "Q3" , "aQ3" )
  dfr$aQ3 <- dfr$Q3 - mean(dfr$Q3, na.rm=TRUE)
  # cat("calc jack Rcpp") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1		
  dfr2 <- data.frame( "index1" = dfr$index1  , "index2"=dfr$index2 )
  dfr2$Q3 <- dfr$Q3
  dfr2$aQ3 <- dfr$aQ3
  # compute p value
  N <- nrow(resp)
  se1 <- - abs( dfr2$aQ3 * sqrt( N -3 ) )
  dfr2$p <- 2 * stats::pnorm( se1  )
  dfr <- dfr2
  dfr <- dfr[ order( dfr$aQ3 , decreasing=TRUE) , ]
  dfr$p.holm <- stats::p.adjust( dfr$p , method="holm")    
   # include sample size of each item pair
   resp_ind <- 1 - is.na(resp)
   cp1 <- crossprod( resp_ind )
   dfr$N_itempair <- cp1[ as.matrix(dfr[ , c("index1" , "index2" ) ]) ]         
  cn <- colnames(resp)
  # cat("p values") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1					
  #**** fit statistic
  stat.MADaQ3 <- data.frame( "MADaQ3" = mean( abs(dfr$aQ3), na.rm=TRUE ) )	
  # cat("fit statistic") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1						
  #**** order item pair table
  dfr <- data.frame( "item1" = cn[ dfr$index1 ] , "item2" = cn[dfr$index2] , dfr )
  dfr <- dfr[ order(dfr$p) , ]		
  stat.MADaQ3$maxaQ3 <- abs(dfr$aQ3[1])
  stat.MADaQ3$p <- dfr$p.holm[1]
  # cat("before calc counts") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1							
  #***** calculate observed and expected counts
  if (progress){
    cat("**** Calculate Counts \n") ; utils::flush.console()
  }
  res1 <- .Call("tam_q3_calc_V2counts" , as.matrix(resp0) , as.matrix(resp.ind) , 
                as.vector(rprobs) , as.matrix(hwt) , maxKi , maxK ,
                PACKAGE="TAM")
  obs_counts <- res1$obs_counts
  exp_counts <- res1$exp_counts
  pair_exists <- rowSums(obs_counts) > 0
    
  # cat("calc counts V2") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1  					
  chi2.stat <- data.frame(res1$maxKiM)
  colnames(chi2.stat) <- c("index1","index2", "maxK1" , "maxK2" , "df")
  eps <- 1e-10
  chi2.stat$chi2[pair_exists] <- rowSums( ( obs_counts - exp_counts )^2 / ( exp_counts + eps ) )[pair_exists]
  chi2.stat$p[pair_exists] <- 1- stats::pchisq( chi2.stat$chi2[pair_exists] , df=chi2.stat$df[pair_exists]  )
  chi2.stat$p.holm[pair_exists] <- stats::p.adjust( chi2.stat$p[pair_exists] , method="holm")
  #*****
  # calculate covariance and correlation
  scorematrix <- cbind( rep( 0:(maxK-1) , maxK ) ,
                        rep(0:(maxK-1) , each=maxK) )
  # observation covariances and correlations
  if (progress){
    cat("**** Calculate Covariances \n") ; utils::flush.console()
  }
  res2 <- .Call("tam_calccov" , obs_counts , scorematrix , adjust_= 1 , PACKAGE="TAM")
  res2e <- .Call("tam_calccov" , exp_counts , scorematrix , adjust_= 0 , PACKAGE="TAM" )
  
  # compute fit statistics
  residcov <- 100*mean( abs( res2$cov_ij - res2e$cov_ij ), na.rm=TRUE )
  #    residcov_jack <- 100*colMeans( abs( res3$cov_ij_jack - res3e$cov_ij_jack ) )
  
  srmr <- mean( abs( res2$cor_ij - res2e$cor_ij ), na.rm=TRUE  )
  #    srmr_jack <- colMeans( abs( res3$cor_ij_jack - res3e$cor_ij_jack ) )
  srmsr <- sqrt( mean( ( res2$cor_ij - res2e$cor_ij )^2, na.rm=TRUE  ) )
  
  fitstat <- c( residcov , srmr , srmsr )
  #	fitstat_jack <- rbind( residcov_jack , srmr_jack )
  #	resjj <- .tam.q3.jackknife2( fitstat , fitstat_jack )	
  names(fitstat) <- c("100*MADCOV" , "SRMR","SRMSR")
  # cat("calccov") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1							
  
  # create matrix of Q3 and adjusted Q3 statistics
  dfr1 <- dfr
  dfr1 <- dfr1[ order( 10000*dfr1$index1 + dfr1$index2 ) , ]
  Q3.matr <- matrix(NA , I , I )	
  rownames(Q3.matr) <- colnames(Q3.matr) <- colnames(resp)
  aQ3.matr <-  Q3.matr
  RR <- nrow(dfr1)
  for (rr in 1:RR){
    ii1 <- dfr1$index1[rr]
    ii2 <- dfr1$index2[rr]		
    Q3.matr[ ii1 , ii2 ] <- Q3.matr[ii2,ii1] <- dfr1$Q3[rr]
    aQ3.matr[ ii1 , ii2 ] <- aQ3.matr[ii2,ii1] <- dfr1$aQ3[rr]		
  }	
  # inclusion of item fit statistics
  chisquare.itemfit <- data.frame("item" = colnames(resp) ,
                                  "index" = 1:I )
  for (ii in 1:I){				
    # ii <- 1
    h1 <- dfr[ dfr$index1 == ii | dfr$index2 == ii , ]
    h1 <- h1[!is.nan(h1$p),]
    chisquare.itemfit$p[ii] <- min( stats::p.adjust( h1$p , method="holm") )
  }		
  chisquare.itemfit$p.holm <- stats::p.adjust( chisquare.itemfit$p , method="holm") 
  
  # maximum chi square
  modelfit.test <- data.frame(
        "maxX2" = max( chi2.stat$chi2) , 
		"Npairs" = nrow(chi2.stat) , 
		"p.holm" = min( chi2.stat$p.holm[pair_exists] ) 
				)
  
  #******
  # modelfit.stat
  modelfit.stat <- fitstat
  modelfit.stat["MADaQ3"] <- stat.MADaQ3["MADaQ3"]
  modelfit.stat["pmaxX2"] <- modelfit.test["p.holm"]
  modelfit.stat <- as.data.frame(modelfit.stat)
  
  
  #****************
  # calculate summary of Q3 statistics
  Q3_summary <- data.frame( "type" = c("Q3" , "aQ3" ) )
  diag(cp1) <- NA
  Q3_summary[1,"M"] <- sum( Q3.matr * cp1 , na.rm=TRUE ) / sum( cp1 , na.rm=TRUE ) 
  Q3_summary[2,"M"] <- sum( aQ3.matr * cp1 , na.rm=TRUE ) / sum( cp1 , na.rm=TRUE ) 
  Q3_summary[1,"SD"] <- sum( Q3.matr^2 * cp1 , na.rm=TRUE ) / sum( cp1 , na.rm=TRUE ) 
  Q3_summary[2,"SD"] <- sum( aQ3.matr^2 * cp1 , na.rm=TRUE ) / sum( cp1 , na.rm=TRUE )   
  Q3_summary$SD <- sqrt( Q3_summary$SD - Q3_summary$M^2 )
  Q3_summary[1,"min"] <- min( Q3.matr , na.rm=TRUE ) 
  Q3_summary[1,"max"] <- max( Q3.matr , na.rm=TRUE ) 
  Q3_summary[2,"min"] <- min( aQ3.matr , na.rm=TRUE ) 
  Q3_summary[2,"max"] <- max( aQ3.matr , na.rm=TRUE )
  cp10 <- 1 - is.na(cp1)    
  Q3_summary[1,"SGDDM"] <- sum( abs(Q3.matr) * cp10 , na.rm=TRUE ) / 
				sum( cp10 , na.rm=TRUE )   
  Q3_summary[2,"SGDDM"] <- sum( abs(aQ3.matr) * cp10 , na.rm=TRUE ) / 
						sum( cp10 , na.rm=TRUE )   
  Q3_summary[1,"wSGDDM"] <- sum( abs(Q3.matr) * cp1 , na.rm=TRUE ) / 
				sum( cp1 , na.rm=TRUE )   
  Q3_summary[2,"wSGDDM"] <- sum( abs(aQ3.matr) * cp1 , na.rm=TRUE ) / 
						sum( cp1 , na.rm=TRUE ) 						

  
  #******* OUTPUT *******
  res <- list( "stat.MADaQ3"=stat.MADaQ3 , "chi2.stat"=chi2.stat ,
               "fitstat" = fitstat , "modelfit.test"= modelfit.test , 
               "stat.itempair"=dfr , 
               "chisquare.itemfit"=chisquare.itemfit , 
               "residuals" = residM ,
               "Q3.matr"=Q3.matr  ,  "aQ3.matr"=aQ3.matr ,
			   Q3_summary = Q3_summary , 
			   N_itempair = cp1 , 
               "statlist" = modelfit.stat )
  class(res) <- "tam.modelfit"
  return(res)
}
################################################
