


###################################################################
lavaanify.sirt.v1 <- function( lavmodel ){
 z0 <- Sys.time()
	syn <- lavmodel
	syn <- strsplit( syn , " " )[[1]]
	syn <- syn[ syn != "" ]
	syn <- gsub( ";" , "\n" , syn )
	
	#*****	
	syn <- split_syn_string( syn , "\\n" )		
	syn[ syn == "\\n" ] <- "\n"
	
	#***
	dfr1 <- data.frame( "index" = 1:length(syn) , "syntax"=syn )
	# look for specific strings and breaks
	dfr1$eqind <- 0
	N1 <- nrow(dfr1)
	vv <- 1
	for (ii in 1:N1){
	#    ii <- 1
		dfr1[ ii , "eqind" ] <- vv
		if ( length( grep( "\n" , dfr1$syntax[ii] ) ) > 0 ){ vv <- vv + 1 }
					}										
	syn0 <- lavmodel
	
	#***************************************************************
	# handling of guessing and slipping parameters						
	dfr1$guess_slip <- 0
	ind <- grep( "\\?=" , dfr1$syntax , perl=FALSE)
	if ( length(ind) > 0 ){
		dfr1$guess_slip[ ind ] <- 1
							}
	eqgroups <- dfr1$eqind[ which( dfr1$guess_slip == 1 ) ]
	dfr1$guess_slip[ dfr1$eqind %in% eqgroups ] <- 1

	
	# create "normal" lavaan syntax
	dfr2 <- dfr1[ dfr1$guess_slip == 0 , ]
	lavmodel1 <- paste0( dfr2$syntax , collapse="")
	
	lavpartable1 <- lavaan::lavaanify( as.character(lavmodel1 ) , warn = FALSE , debug=FALSE ,
						fixed.x=FALSE)
						
    # lavpartable1 <- lavaanify_in_sirt( as.character(lavmodel1 ) , warn = FALSE , debug=FALSE )	
	res1 <- change.grep.lavpartable( lavpartable1 )
		
# cat("**** change grep") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1		
	# create a new model syntax here!!!
    if ( res1$changed ){
	   syn0 <- lavpartable2lavsyntax( res1$lavpartable )	   
	   
# cat("**** lavpartable2lavsyntax") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1			   
	   lavpartable1 <- lavaan::lavaanify( as.character( syn0 ) , warn = FALSE , debug=FALSE ,
						    fixed.x = FALSE)
 
#  cat("**** lavaanify changed sirt.v1") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1		   
						}
	

	# create lavaan parameter table for guessing/slipping parameters
	dfr2 <- dfr1[ dfr1$guess_slip == 1 , ]
	ug <- unique( dfr2$eqind)
	vecstr <- c("\\+" , "\\\n" , "\\?=" , "\\*" )

	
	for (uu in ug){
#		uu <- ug[1]
# cat("\n\n---------------" , uu , "--------------\n")
		syn.temp <- paste0( dfr2$syntax[ dfr2$eqind == uu ] , collapse="")
		syn.temp <- split_syn_string_vec( syn=syn.temp , vecstr = vecstr )
        syn.temp[ syn.temp == "g1" ] <- "t1"
        syn.temp[ syn.temp == "s1" ] <- "t2"
		syn.temp[ syn.temp == "\\?=" ] <- "|"
        syn.temp[ syn.temp == "\\+" ] <- "+"
		syn.temp[ syn.temp == "\\*" ] <- "*"
		syn.temp[ syn.temp == "\\\n" ] <- "\n"
        syn.temp <- paste0( syn.temp , collapse="")
		h1 <- lavaan::lavaanify( syn.temp)	
		h1 <- h1[ h1$op == "|" , ]
        h1$op <- "?="
	
        h1[ h1$rhs == "t1" ,"rhs"] <- "g1"		
        h1[ h1$rhs == "t2" ,"rhs"] <- "s1"	
		h0 <- h1
		h1$free <- h1$free + max(lavpartable1$free)
		h1$free[ h0$free == 0 ] <- 0
#		h1$eq.id <- h1$eq.id + max(lavpartable1$eq.id)
#		h1$eq.id[ h0$eq.id == 0 ] <- 0
#		h1$unco <- h1$unco + max(lavpartable1$unco)
#		h1$unco[ h0$unco == 0 ] <- 0
		h1$id <- h1$id + max(lavpartable1$id)	
		lavpartable1 <- rbind( lavpartable1 , h1 )
				}	
					
#cat("**** guess / slip") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1	
	res <- list("lavpartable" = lavpartable1 , "lavaan.syntax"=syn0 )			
	return(res)	
			}
##################################################################			
