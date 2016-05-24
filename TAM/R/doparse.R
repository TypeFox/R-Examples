
####################################################
# doparse function
doparse <- function(model){
   syn <- model	
   syn <- gsub( ";" , "\n" , paste(syn) )
   syn <- unlist(strsplit( syn , split="\n" ))  
   # delete empty entries
   dels <- c("")
   for (vv in 1:10){ 
		dels <- c( dels , paste0( rep( " " , vv ) , collapse="")  )
		}
   syn <- syn[ ! ( syn %in% dels )]
   
    
   #***************
   # delete lines beginning with "#" and so on
#   lin <- c("#" , " #" , "  #" , "    #")
#   L1 <- 4
   
   ind1 <- which( substring( syn , 1 , 1 ) == "#" )
   ind2 <- which( substring( syn , 1 , 2 ) == " #" )
   ind3 <- which( substring( syn , 1 , 3 ) == "  #" )
   ind4 <- which( substring( syn , 1 , 4 ) == "   #" )   
   ind <- unique( c( ind1 , ind2 , ind3 , ind4 ) ) 
   if ( length(ind) > 0 ){
     syn <- syn[ - ind ]
					}

   #***************
   # process syntax	
   S1 <- length(syn)
   dfr <- data.frame( "index" = 1:S1 , "syn" = syn )
   dfr$DO <- 0
   dfr$DO2 <- 0 
   dfr$DOTYPE <- 0   
   dfr$DOEND <- 0
   
   ind <- grep( "DO(" , syn , fixed=TRUE )
   if ( length(ind) > 0){   dfr$DO[ind] <- 1	}
   ind <- grep( "DO2(" , syn , fixed=TRUE )
   if ( length(ind) > 0){   dfr$DO2[ind] <- 1	}
   ind <- grep( "DOEND" , syn , fixed=TRUE )
   if ( length(ind) > 0){   dfr$DOEND[ind] <- 1	}   
   dfr$DOTYPE <- dfr$DO + dfr$DO2
   
   dfr$DOINDEX <- 0
   N1 <- nrow(dfr)
   vv <- 0
   hh <- 0
   # search for different DO parts in the syntax
   for (ii in 1:N1){
      if ( dfr$DOTYPE[ii] == 1 ){
			vv <- vv + 1
			dfr$DOINDEX[ii] <- vv
			hh <- 1
				}
      if ( dfr$DOEND[ii] == 1 ){
			hh <- 0
				}								
	   if (hh==1){
	       dfr$DOINDEX[ii] <- vv
				}
			}
	# extract conditions of do statements
    NDO <- max( dfr$DOINDEX)
	if (NDO == 0 ){ syn2 <- model }
	if (NDO > 0 ){	
		dfr1 <- dfr[ dfr$DOTYPE > 0 , ]
		dolist <- as.list(1:(NDO+1))
		dotypelist <- as.list(1:(NDO+1))
		for (uu in 1:NDO){
			# uu <- 1
			syn.uu <- gsub( "DO(" , "" , paste( dfr1$syn[uu] ) , fixed=TRUE )
			syn.uu <- gsub( "DO2(" , "" , paste( syn.uu ) , fixed=TRUE )
			syn.uu <- gsub( ")" , "" , paste( syn.uu ) , fixed=TRUE )    
			syn.uu <- unlist(strsplit( syn.uu , split="," , fixed=TRUE ))
			syn.uu <- gsub(" " , "" , syn.uu)
			if ( length(syn.uu) > 3 ){
			    if ( syn.uu[4] == "%1" ){ syn.uu[4] <- "Inf" }
									}
			syn.uu <- as.numeric(syn.uu[ syn.uu != "" ])
			dolist[[uu]] <- syn.uu
			dotypelist[[uu]] <- "DO"
			if ( dfr1$DO2[uu] == 1 ){ dotypelist[[uu]] <- "DO2" }
						}
	   syn2 <- NULL 
	   hh <- 0
	   for (ii in 1:N1){
		  if ( dfr$DOEND[ii] + dfr$DOINDEX[ii] == 0 ){
			  syn2 <- c( syn2 , paste(dfr$syn[ii]) )
					}
		  if ( ( dfr$DOTYPE[ii] == 1 ) & ( dfr$DOINDEX[ii] > 0 ) ){			
				hh <- hh+1	
					}
		  if ( ( dfr$DOTYPE[ii]  == 0 ) & ( dfr$DOINDEX[ii] > 0 ) ){			
				syn2.ii <- paste(dfr$syn[ii])
				if ( dotypelist[[hh]] == "DO" ){
				   dom <- dolist[[hh]]
				   for ( zz in seq( dom[1] , dom[2] , dom[3] ) ){
					   syn2 <- c( syn2 , gsub( "%" , zz , syn2.ii) )			   
									}
								}

				if ( dotypelist[[hh]] == "DO2" ){
				   dom <- dolist[[hh]]
				   for ( zz in seq( dom[1] , dom[2] , dom[3] ) ){
				   if ( dom[4] == Inf ){ bb <- zz+1 } else { bb <- dom[4] }
				   for ( zz2 in seq( bb , dom[5] , dom[6] ) ){
					   syn2 <- c( syn2 , gsub( "%2" , zz2 , gsub( "%1" , zz , syn2.ii) ) )			   
									}
								}
							}			
								
								
							}			
					}	
		syn2 <- paste0( syn2 , collapse="\n" )				
			}
    return(syn2)
		}
####################################################		