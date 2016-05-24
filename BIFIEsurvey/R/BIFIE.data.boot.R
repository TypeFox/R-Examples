###########################################################
# BIFIE.data objects for bootstrap designs
BIFIE.data.boot <- function( data , wgt=NULL ,  pv_vars = NULL ,
	     Nboot = 500 , seed = .Random.seed , cdata=FALSE
		){	
	cl <- match.call()		
	#*** list of multiply imputed datasets
	if ( ( is.list(data) ) & ( ! is.data.frame(data) ) ){
		dataL <- data
		data <- dataL[[1]]
				}  else {
		dataL <- data
				}					
    data <- as.data.frame( data )		
	if ( ! is.null(seed) ){
		set.seed( seed )
				}
	
	
	# normalize weights
	if ( is.null(wgt) ){
		wgt <- "_wgt2"
		data[ , wgt ] <- rep( 1 , nrow(data) )
			}
	wgtname <- wgt
	wgt <- data[ , wgt ]
	N <- length(wgt)
	wgt <- N * wgt / sum(wgt)
	# cumulated weights
	cumwgt <- cumsum(wgt)
	# random numbers
	rand_wgt <- N*matrix( stats::runif(N*Nboot) , nrow=N , ncol=Nboot ) 
	
	#****
	# apply bootstrap subfunction
	cat("+++ Generate bootstrap samples\n"); utils::flush.console()
	# datarep <- bifie_boot( cumwgt , rand_wgt )$wgtM
	datarep <- .Call("bifie_boot", cumwgt , rand_wgt , PACKAGE="BIFIEsurvey")$wgtM
	RR <- Nboot
	addname <- 10^( floor( log( RR+.5 , 10 ) )  + 1 )
	colnames(datarep) <- paste0("w_fstr" , substring( paste0(addname +1:RR),2) )	
	datarep <- sum( wgt )/N * datarep
	
#	cat("|\n")
	#******** generate replicated datasets for datasets
	if ( is.null( pv_vars) ){ datalist <- dataL  }
	if ( ! is.null( pv_vars )){
		dfr <- NULL
		VV <- length(pv_vars)
		for (vv in 1:VV){
			vv1 <- pv_vars[vv]
			ind.vv1 <- which( substring( colnames(data) , 1 , nchar( vv1 ) ) == pv_vars[vv] )
			Nimp <- length(ind.vv1)
			dfr2 <- data.frame( "variable" = vv1 , "var_index" = vv , "data_index" = ind.vv1 ,
							 "impdata_index"=1:Nimp ) 			
			dfr <- rbind( dfr , dfr2 )
			}
		sel_ind <- setdiff( 1:( ncol(data) ) , dfr$data_index )
		data0 <- data[ , sel_ind ]
		V0 <- ncol(data0)
		newvars <- seq( V0+1 , V0+VV )
		datalist <- as.list( 1:Nimp )
		for (ii in 1:Nimp ){
			dat1 <- data.frame( data0 , data[ , dfr[ dfr$impdata_index == ii  , "data_index" ] ] )
			colnames(dat1)[ newvars ] <- pv_vars
			datalist[[ii]] <- dat1 
						}
					}			
	fayfac <- 1/Nboot
	#*** create BIFIE.data object
	bifiedat <- BIFIE.data( datalist , wgt = data[, wgtname ] , 
			      wgtrep = datarep , fayfac = fayfac , cdata=cdata)
	bifiedat$CALL <- cl				  
	return(bifiedat)
}
###############################################################################