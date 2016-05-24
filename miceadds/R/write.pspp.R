write.pspp <- function (data , datafile, pspp.path , decmax=6 ,
  as.factors=TRUE ,  use.bat=FALSE) {
    data <- as.data.frame(data) 
    df <- data 
	codefile <- paste0( datafile , ".sps" )
    adQuote <- function(x){paste("\"", x, "\"", sep = "")}
	varnames <- colnames(df)
	if ( as.factors ){ 
		 dfn <- lapply(df, function(x) if (is.factor(x))
			 as.numeric(x)
			 else x) } else {
		 dfn <- lapply(df, function(x) 
			if (is.factor(x)) paste(x)  else x
					) }
#     write.table(dfn, file = paste0(datafile,".dat") , row = FALSE, col = FALSE , na =".")
     if(is.null(attributes(df)$variable.labels)) varlabels <- names(df) else varlabels <- attributes(df)$variable.labels
     if (is.null(varnames)) {
         varnames <- abbreviate(names(df), 8)
         if (any(sapply(varnames, nchar) > 8))
             stop("I cannot abbreviate the variable names to eight or fewer letters")
         if (any(varnames != names(df)))
             warning("some variable names were abbreviated")
     }
#     cat("DATA LIST FILE=", dQuote(datafile), " list\n", file = codefile)
	# log dfn
	eps2 <- .001
	ldfn <- lapply( dfn ,  FUN = function(vv){
			if ( is.numeric(vv) ){	
				floor( max( log(abs(vv)+1,10) , na.rm=TRUE ) )+2  } else { max( nchar(vv) ) }
						}  
						)
	# number of decimals after digits
	V <- length(dfn)
	lafter_dfn <-  rep(0,V)	
	stringentry <- rep(0,V)
	for (vv in 1:V){
	  if ( is.numeric( as.vector(dfn[[vv]] ) ) ){
		dvv <- abs( as.numeric( as.vector(dfn[[vv]] )) )
		dd <- 0
		hh <- 1
		while( ( hh == 1 ) & ( dd < decmax ) ){	
			yvv <- 10^dd * dvv  - floor( 10^dd* dvv )
			if ( max(yvv,na.rm=TRUE) == 0 ){ hh <- 0 ; break } else { dd <- dd+1 }
					}
					} else { 
				dd <- max( nchar( paste(dfn[[vv]] ))) 
				stringentry[vv] <- 1
					}
		lafter_dfn[vv] <- dd
						
				}
#	pformat <- paste0( "F" , max( 1, unlist(ldfn) + 1 + lafter_dfn ) , "." , lafter_dfn , "" )	
	xf <- unlist(ldfn) + 1 + lafter_dfn
	xf <- ifelse( xf == "0" , "1" , xf )
    pformat <- paste0( "F" , xf , "." , lafter_dfn , "" )	
	pformat <- ifelse( stringentry==1 , paste0( "A" , lafter_dfn  ) , pformat )
	
	vars2 <- paste( paste( varnames , pformat ) , collapse="\n " )

	dfn1 <- as.data.frame( dfn )
	
	utils::write.csv2( dfn1 , paste0( datafile , ".csv" ) , row.names=FALSE , 
			quote= FALSE, na ="")
	
	
		
#     cat(paste0( "DATA LIST FILE='", gsub( "\\" , "//" , getwd(), fixed=TRUE )  , "/" , datafile , 
#				".dat' list\n" ) , file = codefile)
     cat(paste0( "GET DATA \n /TYPE=TXT \n /FILE='", gsub( "\\" , "//" , getwd(), fixed=TRUE )  , "/" , datafile , 
				".csv' \n" ,
				"/IMPORTCASES=ALL\n" ,
				"/ARRANGEMENT=DELIMITED\n" ,
				"/DELCASE=LINE\n" ,				
				"/FIRSTCASE=2\n" ,								
				"/DELIMITERS=';'\n" ,												
				"/QUALIFIER=''\n" ,					
#				"/QUALIFIER=\"\"\n" ,					
				"/ESCAPE \n /VARIABLES= \n" 
				) , file = codefile)
     cat( paste0( vars2 , " .\n\n" ), file = codefile, append = TRUE)
     cat("VARIABLE LABELS\n", file = codefile, append = TRUE)
     cat(paste(varnames, adQuote(varlabels), "\n"), ".\n", file = codefile,
         append = TRUE)
#     factors <- sapply(df, is.factor)
	if ( as.factors){
		 factors <- sapply(dfn1, is.factor)
		 if (any(factors)) {
			 for (v in which(factors)) {
			 cat("\nVALUE LABELS", file = codefile, append = TRUE)
				 cat("\n", file = codefile, append = TRUE)
				 cat(varnames[v], " \n", file = codefile, append = TRUE)
				 levs <- levels(df[[v]])
				 cat(paste(1:length(levs), adQuote(levs), "\n", sep = " "),
					 file = codefile, append = TRUE)
			 }
			 cat(".\n", file = codefile, append = TRUE)
		 }
		 }
	 cat("\n",file=codefile,append=TRUE)
	 ###########################################################
	 # write value labels
	 varnames <- colnames(df)
	 for (vv in varnames){ 
		# vv <- varnames[1]
		avv <- attr( df[,vv] , "value.labels" )
		if ( length(avv) > 0 ){				
			cat("VALUE LABELS\n" , file=codefile, append=TRUE )
			pvv <- paste0( avv , " '" ,  names(avv)  , "'" )	
			pvv <- paste( paste( vv , paste( pvv , collapse=" ") )  )
			cat( pvv , ".\n", file = codefile, append = TRUE)
					}
					
					
				}
	 ############################################################
     cat("\nEXECUTE.\n", file = codefile, append = TRUE) 
	 	 
	 cat( paste0( "\n\n save outfile='"  , getwd() , "/" , datafile , ".sav'.\n execute.") ,
			file=codefile , append=TRUE )
	#****	
	# run PSPP		
	p1 <- paste0( "\"" , pspp.path , "pspp.exe\" "  , codefile  )
	if ( use.bat ){
		base::writeLines( p1 , "_batch_pspp.bat" )
		system( "_batch_pspp.bat" )			
				} else {
		system( p1 )
				}
        }
