
#############################################
# Q3 statistic based on weighted likelihood estimates
tam.Q3 <- function( tamobj , ... ){
    object <- tamobj
	res1 <- IRT.residuals( object , ...)	
	resp_ind <-  1 - is.na( res1$residuals )
	cp20 <- cp1 <- crossprod( resp_ind )
	# calculate Q3
	Q3.matr <- stats::cor( res1$residuals , use="pairwise.complete.obs" )
	diag(Q3.matr) <- NA
	aQ3.matr <- Q3.matr - mean( Q3.matr , na.rm=TRUE )	
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
	
	
	
	res <- list( Q3_summary = Q3_summary , "Q3.matr" = Q3.matr , 
				aQ3.matr = aQ3.matr , 
			"N_itempair" = cp20 , "residuals" = res1 )
	class(res) <- "tam.Q3"		
	return(res)
		}
#############################################	