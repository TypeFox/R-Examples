

#######################################################################
# Histogram
BIFIE.hist <- function( BIFIEobj , vars , breaks=NULL ,
		group=NULL , group_values=NULL  ){
	#****
	s1 <- Sys.time()
	cl <- match.call()		
	bifieobj <- BIFIEobj
	if (bifieobj$cdata){
		varnames <- unique( c( vars , group , "one") )
		bifieobj <- BIFIE.BIFIEcdata2BIFIEdata( bifieobj , varnames=varnames )	
						}							
	FF <- Nimp <- bifieobj$Nimp
	N <- bifieobj$N
	dat1 <- bifieobj$dat1
	wgt <- bifieobj$wgt
	wgtrep <- bifieobj$wgtrep
	varnames <- bifieobj$varnames
	RR <- bifieobj$RR
	datalistM <- bifieobj$datalistM
    fayfac <- bifieobj$fayfac	
	vars <- vars[1]	
	vars_index <- unlist( sapply( vars , FUN = function(vv){ 
						which( varnames == vv ) } ) )
	
	if ( is.null(breaks) ){
	    x <- dat1[ , vars_index ]
		breaks <- base::pretty(x , n = grDevices::nclass.Sturges(x))
				}
	
		RR <- 0
#	if (RR==1){ RR <- 0 }
#	if ( ! se ){ 
#		wgtrep <- matrix( wgt , ncol=1 )
#		RR <- 0
#				}	
	
    # vars values
	VV <- length(vars)
					
	wgt_ <- matrix( wgt , ncol=1 )
	if ( is.null( group) ){ nogroup <- TRUE } else { nogroup <- FALSE }
	cat(paste0( "|" , paste0( rep("*" , FF) , collapse="") , "|\n" ))
	if (nogroup){
	    group <- "one"
	    group_values <- c(1)
			}
			

	#@@@@***
    group_index <- match( group , varnames )
	#@@@@***

    if ( is.null(group_values ) ){ 
		t1 <- fasttable( datalistM[ , group_index ] )				  
	    group_values <- sort( as.numeric( paste( names(t1) ) ))
				}
	
	#@@@@***
	res00 <- BIFIE_create_pseudogroup( datalistM , group , group_index , group_values )				
	res00$datalistM -> datalistM 
	res00$group_index -> group_index
	res00$GR -> GR 
	res00$group_values -> group_values
	res00$group -> group
	#@@@@***			

				
				
	#**************************************************************************#
	# Rcpp call
	res <- .Call("bifie_hist" ,  datalistM , wgt_ , wgtrep , 
				vars_index - 1,    fayfac ,
				Nimp , group_index - 1 , group_values , breaks , PACKAGE="BIFIEsurvey")				
	
	# create histogram objects
	GG <- length(group_values)
	histobj <- list(1:GG)
	BB <- res$BB
	
	for (gg in 1:GG){
		h1 <- list( "breaks" = res$breaks ,
				"counts" = res$sumwgt[ ( gg-1)*BB + 1:BB ] ,
				"density" = res$density_vec[ ( gg-1)*BB + 1:BB ] ,
				"mids" = res$mids	
						)
		h1$xname <- paste0( vars , "_" , group , group_values[gg] )
		if ( stats::sd ( diff(res$mids) ) < .000001 ){ h1$equidist <- TRUE } else { h1$equidist <- FALSE }
		class(h1) <- "histogram"	
		histobj[[gg]] <- h1
		# histobj[[gg]] <- h1$xname
					}
	names(histobj) <- paste0( vars , "_" , group , group_values )	
	
	
	
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "histobj" = histobj , 
			"output" = res ,
			"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
			"NMI" = BIFIEobj$NMI , "Nimp_NMI" = BIFIEobj$Nimp_NMI , 
			"GG"=GG , "CALL"= cl)
	class(res1) <- "BIFIE.hist"
	return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.hist function
summary.BIFIE.hist <- function( object ,  ... ){
    BIFIE.summary(object)
			}
##########################################################################			
# plot function
plot.BIFIE.hist <- function( x , ask=TRUE , ... ){
    res <- x
    GG <- res$GG
    for (gg in 1:GG){
      graphics::plot( res$histobj[[gg]] , ... )
      graphics::par(ask=ask)
                    }
      }
			