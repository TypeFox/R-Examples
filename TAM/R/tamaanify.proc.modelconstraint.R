


#########################################################
# process model constraint
tamaanify.proc.modelconstraint <- function( res ){ 	
	tam1 <- res$tammodel.dfr
	ind1 <- which( paste(tam1$syn) == "MODELCONSTRAINT:" )
	index1 <- tam1$part_begin[ ind1 ]
	syncon <- paste( tam1[ which( tam1$part_begin == index1 )[-1] , "syn" ] )
	# extract newly defined parameters
	res$MODELCONSTRAINT <- syncon
	lavpartable <- res$lavpartable

	dfr.syncon <- NULL
	
	if ( length(syncon)>0){
		dfr.syncon <- data.frame( "index" = seq(1,length(syncon) ) ,
						"syn" = syncon )						
		dfr.syncon$derived <- 0
		dfr.syncon$trafopar <- ""		
		ind2 <- grep("==",dfr.syncon$syn )
		if ( length(ind2) > 0 ){
			dfr.syncon[ ind2 , "derived"] <- 1
			v1 <- strsplit( syncon[ind2] , split="==" , fixed=TRUE )
			v1 <- unlist( lapply( v1 , FUN = function(vv){ vv[1] } ) )
			dfr.syncon[ ind2 , "trafopar" ] <- v1				
						}
						
		dfr.syncon <- dfr.syncon[ order( dfr.syncon$trafopar ) , ]
		dfr.syncon <- dfr.syncon[ order(dfr.syncon$index) , ]
        dfr.syncon$syn <- paste(dfr.syncon$syn)
        dfr.syncon$trafopar <- paste(dfr.syncon$trafopar)

		#******************
        # add "==" operator if they are not included in model constraints
		# grep( paste(dfr.syncon$syn) )
		N <- nrow(dfr.syncon)
		vv0 <- paste0( dfr.syncon$trafopar)[1]
		ii0 <- 1
		del0 <- NULL
		for (ii in 2:N){
			ind <- grep(  "==" ,  paste0( dfr.syncon$syn )[ii] ,fixed=TRUE )
			if ( length(ind) == 0 ){
			   dfr.syncon[ii0,"syn"] <- paste0( dfr.syncon$syn[ii0] , 
										    paste0( dfr.syncon$syn )[ii] )
			   del0 <- c( del0 , ii )			  
						} else {
				ii0 <- ii
						}
				      }
					  
         if ( length(del0) > 0 ){
		      dfr.syncon <- dfr.syncon[ - del0 , ]
								}
						
								
		#*************
		# if there is some "__", then transform dfr.syncon
		# Revalpr("lavpartable")        
		ind2 <- dfr.syncon[ grep( "__" , paste( dfr.syncon$trafopar )  , fixed=TRUE ) , "index" ]
		labs <- unique( paste(lavpartable$label) ) 
		dfr.syncon$trafopar <- paste(dfr.syncon$trafopar)
		dfr.syncon$expanded <- 0	
		
        if ( length(ind2)>0){   
		dfr.syncon$syn_orig <- dfr.syncon$syn		
        for (ii in ind2){
	#		ii <- ind2[1]
			parms.ii <- unlist( strsplit( paste( dfr.syncon$trafopar )[ dfr.syncon$index == ii  ] , split="__" ) )
			
			ind3 <- c( which( labs == parms.ii[1] ) , which( labs == parms.ii[2] ) )
			ind3 <- seq( ind3[1] , ind3[2] )
			labs_sub <- labs[ ind3 ]
			LS <- length(labs_sub)
			dfr.syn1 <- dfr.syncon[ dfr.syncon$index == ii , , drop=FALSE ]
			dfr.syn1 <- dfr.syn1[ rep(1,LS) , ]
			dfr.syn1$trafopar <- labs_sub
			dfr.syn1$expanded <- 1
			m1 <- unlist( strsplit( paste(dfr.syn1[1,"syn"]) , split="==" , fixed=TRUE) )[2]			
			dfr.syn1$syn <- paste0( dfr.syn1$trafopar , "==" , m1 )
   		    rownames(dfr.syn1) <- NULL
			dfr.syncon <- rbind( dfr.syncon[ dfr.syncon$index != ii , , drop=FALSE ] , dfr.syn1 )
							}
						}
		lavpartable <- res$lavpartable
		lab1 <- paste(dfr.syncon$trafopar)
		lab1 <- lab1[ lab1 != "" ]
		lav1 <- lavpartable[ paste(lavpartable$label) %in% lab1 ,
						c("lhs","op","rhs","label") ]			
		lav1 <- lav1[ order( lav1$label ) , ]
		dfr.syncon$rhs <- dfr.syncon$op <- dfr.syncon$lhs <- ""
		dfr.syncon[ dfr.syncon$derived == 1, c("lhs","op","rhs") ] <- lav1[,1:3]
		dfr.syncon$fullsyn <- paste0( dfr.syncon$lhs , dfr.syncon$op , dfr.syncon$rhs)
		
					}		

	res$MODELCONSTRAINT.dfr <- dfr.syncon	
	return(res)
	   }
######################################################################	   
