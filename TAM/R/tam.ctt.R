tam.ctt <-
function( resp , wlescore=NULL, pvscores=NULL , group=NULL , progress=TRUE){
    I <- ncol(resp)
	resp0 <- resp
	wlescore0 <- wlescore
	pvscores0 <- pvscores
    # define progress bar
	if ( is.null(wlescore) & is.null(pvscores) ){			
		stop("Provide WLEs or Plausible Values!") }
	if ( is.null(group) ){ group = rep(1 , nrow(resp) ) }
    groups <- sort( unique( group )	)
    G <- length(groups)
    dfr <- NULL
	for (gg in 1:G){
	    if (progress){
			cat("|")
			cat( paste( rep("*" , 10 ) , collapse="") )
			cat("| Group" , groups[gg] , "\n|")
			prg <- round( seq( 1 , I , len=10 ) )
			prg[ prg == I ] <- I-1
            }	
		ind.gg <- which( group == groups[gg] )
		resp <- resp0[ ind.gg , ]
		wlescore <- wlescore0[ ind.gg ]
		if ( !is.null(pvscores0)){ pvscores <- pvscores0[ ind.gg, ]	}
		for (ii in 1:I){
			# ii <- 1
			categ.ii <- table( resp0[,ii] )
			categ.ii <- names(categ.ii)		
			CCI <- length(categ.ii)
			dfr.ii <- data.frame( "group"=groups[gg] , "itemno" = ii , "item" = colnames(resp)[ii] ,
				"N" = sum(1-is.na( resp[,ii] ) )  , "Categ"=categ.ii)
			for (cc in 1:CCI){	
				dfr.ii$AbsFreq[cc]	<- sum( 1 * ( paste(resp[ , ii ] ) == categ.ii[cc]  )	, na.rm=TRUE )
								}
				# ,  categ.ii )
#			colnames(dfr.ii)[4:5] <- c("Categ" , "AbsFreq")
			dfr.ii$RelFreq <- dfr.ii$AbsFreq / dfr.ii$N
			dfr.ii$rpb.WLE <- NA
			dfr.ii$M.WLE <- NA
			dfr.ii$SD.WLE <- NA
			dfr.ii$rpb.PV <- NA
			dfr.ii$M.PV <- NA
			dfr.ii$SD.PV <- NA        			
			for (cc in 1:CCI){
				# cc <- 1
				categ.cc <- categ.ii[cc]
				dat.cc <- 1 * ( resp[ , ii ] == categ.cc  )			
				ind.cc <- which( paste(dat.cc) == 1 )
# print( paste( ii , cc ) ) ; flush.console()						
# print(ind.cc)				
				if ( ! is.null( wlescore ) ){
					dfr.ii[cc , "rpb.WLE"] <- stats::cor( dat.cc , wlescore , use="pairwise.complete.obs")
					dfr.ii[cc , "M.WLE" ] <- mean( wlescore[ ind.cc] , na.rm=TRUE )
					dfr.ii[cc , "SD.WLE" ] <- stats::sd( wlescore[ ind.cc] , na.rm=TRUE )
								}
				if ( ! is.null( pvscores ) ){
					dfr.ii[cc , "rpb.PV"] <- rowMeans( stats::cor( dat.cc , pvscores , use="pairwise.complete.obs") )
					dfr.ii[cc , "M.PV" ] <- mean( apply( pvscores[ ind.cc ,] , 2,mean , na.rm=TRUE ) )
					dfr.ii[cc , "SD.PV" ] <- mean( apply( pvscores[ ind.cc ,] , 2, stats::sd , na.rm=TRUE ) )
								}
					}
			dfr <- rbind( dfr , dfr.ii )
			if (progress){ if ( ii %in% prg ){ cat("-" ) } ; flush.console() }
				} # end items within group
    if (progress){ cat("|\n") }
			} # end group
    dfr <- dfr[ order( paste0( 10000+ dfr$itemno , dfr$group , dfr$Categ ) ) , ]
    dfr <- data.frame( "index" = seq(1,nrow(dfr) ) , dfr )
    return(dfr)
        }
