
############################################################
# Item fit according to the S-X^2 statistic 
# Orlando & Thissen (2000, 2003)
itemfit.sx2 <- function( object , Eik_min=1 , progress=TRUE ){
	mod <- object 
	#*****************************
	# object of class gdm
	if ( class(mod) == "gdm" ){
#		n.ik <- mod$n.ik[,,,1]
		if ( dim( mod$n.ik)[4] > 1){ 
			stop("Only applicable in One-group case!")
							}
		pi.k <- mod$pi.k[,1]
		pjk <- mod$pjk
		data <- mod$data
		I <- ncol(data)
		if ( mod$irtmodel == "1PL" ){ npars <- rep( 1, I ) }
		if ( mod$irtmodel == "2PL" ){ npars <- rep( 2, I ) }		
				}
	#*****************************
	# object of class smirt (sirt)
	if ( class(mod) == "smirt" ){
#		n.ik <- mod$n.ik[,,,1]
		pi.k <- mod$pi.k[,1]
		pjk <- aperm( mod$probs , c(3,1,2) )
		D <- ncol( object$Qmatrix )
		dimQ <- dim(object$Qmatrix)
		Qmatrix <- matrix( object$Qmatrix , nrow= dimQ[1] , ncol=dimQ[2] )
		irtmodel <- object$irtmodel
		if ( is.null( mod$se.a ) ){ 
		   if (irtmodel == "noncomp"){ 	
				npars <- rowSums(Qmatrix)
								}		   
		   if (irtmodel != "noncomp"){ 	
				npars <- 1+0*rowSums(Qmatrix)
								}		 								
					}
		if ( ! is.null( mod$se.a ) ){ 
		   if (irtmodel == "noncomp"){ 	
				npars <- 2*rowSums(Qmatrix)
								}		   
		   if (irtmodel != "noncomp"){ 	
				npars <- 1+rowSums(Qmatrix)
								}		
					}
		data <- mod$dat
				}
	#*****************************************
	# rasch.mml (in sirt)
	if ( class(object) == "rasch.mml" ){
		pi.k <- object$trait.distr[,2]
		pjk0 <- object$pjk
		pjk <- array( 0 , dim=c( dim(pjk0) , 2 ) )
		pjk[,,2] <- pjk0
		pjk[,,1] <- 1 - pjk0
		data <- object$dat
		npars <- 1*( object$item$est.a > 0 ) + 1*( object$item$est.b > 0 ) +
			1*( object$item$est.c > 0 ) + 1*( object$item$est.d > 0 )				
					}	
    #*********************************************
	# din (in CDM)	
	if (class(object)=="din"){					
		data <- object$data
		pi.k <-  object$attribute.patt$class.prob
		# $ pjk                    : num [1:9, 1:2, 1:8] 0.914 0.891 0.871 0.774 0.94 ...		
		pjk <- aperm( object$pjk , c(3,1,2) )
		npars <- rep(2,ncol(data))					
						}		
    #*********************************************
	# gdina (in CDM)	
	if (class(object)=="gdina"){					
		data <- object$data
		pi.k <-  object$attribute.patt$class.prob
		pjk <- aperm( object$pjk , c(3,1,2) )
		npars <- unlist( lapply( object$delta , FUN = function(ll){ length(ll) } ) )				
						}	
    #***********************************************
	# tam.mml (in TAM)
	if (class(object)=="tam.mml"){					
		data <- object$resp
		I <- ncol(data)
		pi.k <-  object$pi.k
		if ( is.matrix(pi.k) ){			
			if (ncol(pi.k) > 1){
				cat("Used first group for assessment of fit.\n")
					}
			pi.k <- pi.k[,1]
				}
		pjk <- aperm( object$rprobs , c(3,1,2) )
		npars <- rep(1,I)
		if ( object$irtmodel == "2PL"){ npars <- rep(2,I) }
						}	

						
    #************************************						
    #************************************
	# data preparation	
	I <- ncol(data)
	sumscore <- rowSums( data )
	N <- nrow(data)
	P1 <- pjk[,,2]
	Q1 <- pjk[,,1]	
	pi.kI <- matrix( pi.k , nrow= dim(pjk)[1] , ncol=dim(pjk)[2]+1 )
	# check input data
	if ( sum( is.na(data) ) > 0 ){ stop("No missing responses are allowed!") }
	if ( max(data)  > 1 ){ stop("Only dichotomous responses are allowed!") }	
	# distribution sum scores
	sumscore.distribution <- sapply( 0:I , FUN = function(ss){ sum( sumscore == ss)  } )
	# score distribution
	scoredistribution <- .calc.scoredistribution.cdm( pjk )
	itemtable <- NULL
	itemfit.stat <- data.frame( "item" = colnames(data) , "itemindex" = 1:I  )

	if (progress){ 
		i3 <- c(1,diff( floor( 10 * ( 1:I )/ (I+1) )+1 ))
		cat( paste0( "|" , paste0( rep("*" , 10) , collapse="") , "|\n|") )
				}
	# calculate fit for item ii
	eps <- 10^(-10)
	for (ii in 1:I){
		  #ii <- 3
#		res <- .calc.itemfit.oneitem( data , ii , pjk , pi.k , P1 , I , 
		res <- .calc.itemfit.oneitem( ii , pjk , pi.k , P1 , I , 
				Eik_min , sumscore.distribution,  scoredistribution , data , sumscore )
		itemtable <- rbind( itemtable , res$table2.ii )

		r1 <- res$table2.ii
		itemfit.stat[ ii , "S-X2" ] <- sum( r1$Nik * ( r1$oik - r1$eik )^2 / ( r1$eik * ( 1 - r1$eik)  + eps ) )
		itemfit.stat[ ii , "df" ] <- nrow(r1) - npars[ii]
		itemfit.stat[ii,"p"] <- 1 - stats::pchisq( itemfit.stat[ ii , "S-X2" ] , df= itemfit.stat[ ii , "df" ] )		
		itemfit.stat[ ii , "S-X2_df" ] <- itemfit.stat[ ii , "S-X2" ] /     itemfit.stat[ ii , "df" ]    
		xg <-  itemfit.stat[ ii , "S-X2" ] - itemfit.stat[ ii , "df" ]    
		itemfit.stat[ ii , "RMSEA" ] <- sqrt( (  ifelse( xg > 0 , xg , 0 )    ) / ( N - 1) / itemfit.stat[ ii , "df" ]  )
		itemfit.stat[ii,"Nscgr"] <- nrow(r1)
		itemfit.stat[ii,"Npars"] <- npars[ii]
			
		if (progress){ if (i3[ii] == 1 ){ cat("-") ; utils::flush.console() } }
				}
	if (progress){ cat("|\n") }
    itemfit.stat[,"p.holm"] <- stats::p.adjust( itemfit.stat[,"p"] , method="holm")			
	res <- list( "itemfit.stat" = itemfit.stat , "itemtable" = itemtable , "I" = I )
	class(res) <- "itemfit.sx2"
	return(res)
	}
#############################################################################
# summary of item fit
summary.itemfit.sx2 <- function(object,...){
    itemfit.stat <- object$itemfit.stat
    i1 <- itemfit.stat
    i1[,-1] <- round( itemfit.stat[,-1] , 3 )
	cat("Please check degrees of freedom (number of estimated paramters) carefully!\n")
	cat("They are maybe not correctly calculated.\n")	
	cat("******    df = Nscgr - Npars   ******* \n\n")
    print(i1)	
	cat("\n-- Average Item Fit Statistics --\n")
	cat( paste0( "S-X2 = " , round( mean( itemfit.stat[,"S-X2"] ) , 3 ) ) )
	cat( paste0( " | S-X2_df = " , round( mean( itemfit.stat[,"S-X2_df"] ) , 3 ) )  , "\n")
	
        }
#############################################################################


###################################################################################
plot.itemfit.sx2 <- function(x,ask=TRUE,...){
    object <- x
    itemtable <- object$itemtable
    itemfit.stat <- object$itemfit.stat
    I <- object$I
    # loop over all items
    for (ii in 1:I){        
        # ii <- 3
        descii <- itemfit.stat[ ii , ]
        title.ii <- paste0( "Item " , descii$item , " | S-X2(df=" , descii$df , ") = " , round( descii[ , "S-X2"] , 3) , 
                        ", p = " , round( descii$p,3) ,  "\n S-X2/df = " , round( descii[,"S-X2_df"] , 3 )  ,
                        " | RMSEA = " , round( descii[,"RMSEA"] , 3 ) )
        itemtable.ii <- itemtable[ itemtable$itemindex == ii , ]
        graphics::plot( itemtable.ii$score , itemtable.ii$oik , xlim = c(1,I-1) , ylim=c(0,1) , 
                type="o" , pch=16 , xlab = "Score group" , ylab="Observed and expected frequency" ,
                main = title.ii)
        graphics::lines( itemtable.ii$score , itemtable.ii$eik , lty=2 , lwd=2)
        graphics::legend( 1 , 1 , c("Observed frequency" , "Expected frequency" ) , pch=c(16,NA) ,
                lty=1:2 , lwd=c(1,2) )
        graphics::par( ask=ask)
        }
    }
###################################################################################
