###########################################################
# BIFIE.data objects for designs with jackknife zones
BIFIE.data.jack <- function( data , wgt=NULL , jktype="JK_TIMSS" , pv_vars = NULL ,
	jkzone=NULL , jkrep=NULL  ,	jkfac=NULL , fayfac = NULL , 
	wgtrep = "W_FSTR" , pvpre = paste0("PV",1:5) , ngr = 100 , 
	seed = .Random.seed ,
	cdata=FALSE 
		){	
	cl <- match.call()	
	
	# subroutine for preparation of nested multiple imputations
	# res0 <- BIFIE_data_nested_MI( data.list=data.list , NMI=NMI )
	# data.list <- res0$data.list
	# Nimp_NMI <- res0$Nimp_NMI	
	
	fayfac0 <- fayfac
	
	if ( ( ! is.null(wgtrep) ) & ( is.null(fayfac) ) ){
		fayfac <- 1 
				}
				
	
	#*** list of multiply imputed datasets
	if ( ( is.list(data) ) & ( ! is.data.frame(data) ) ){
		dataL <- data
		data <- dataL[[1]]
				}  else {
		dataL <- data
				}					
    data <- as.data.frame( data )		
	#*********************************************************
	# using fixed jackknife zones
	if (jktype == "JK_GROUP"){
		N <- nrow(data)
		if ( is.null(wgt) ){
			data$wgt <- rep(1,N)
			wgt <- "wgt"
				}
		data$jkrep <- rep(0,N)
		jkrep <- "jkrep"
		fayfac <- ngr / ( ngr - 1 )
        jkfac <- 0	 						
							}

	
	#**********************************************************	
	#*** defaults for jackknife creation: random groups
	if (jktype == "JK_RANDOM"){	
		N <- nrow(data)
		if ( is.null(wgt) ){
			data$wgt <- rep(1,N)
			wgt <- "wgt"
				}
		if ( ! is.null(seed) ){
			set.seed( seed )
			indzone <- base::sample(1:N) 
					} else {
			indzone <- 1:N
				}
		jkzone <- 1:N
		N1 <- N / ngr
		jkzone <- floor( jkzone / ( N1 + .00001  ) ) + 1
		jkzone <- jkzone[indzone]
		jkrep <- rep(0,N)
		data$jkzone <- jkzone
		jkzone <- "jkzone"
		data$jkrep <- jkrep
		jkrep <- "jkrep"
		fayfac <- ngr / ( ngr - 1 )
        jkfac <- 0	     
				}
							
	#**********************************************************
	#**** defaults for TIMSS	
	if (jktype == "JK_TIMSS"){
	      if ( is.null( jkrep) ){
			jkrep <- "JKREP"
						}
		  if ( is.null( jkzone ) ){ 
				jkzone <- "JKZONE"
							}
		  if ( is.null( wgt ) ){ 							
			wgt <- "TOTWGT"
					}
		  jkfac <- 2
				}
	#***********************************************************
	#**** defaults for PISA
	if (jktype == "RW_PISA"){
	      jkrep <- NULL
	      jkzone <- NULL
		  if ( is.null(wgt)){
			wgt <- "W_FSTUWT"
					}				
		  jkfac <- NULL
		  repvars <- grep( wgtrep , colnames( data ) )
		  RR <- length(repvars)
		  # select variables with plausible values
		  nc1 <- nchar( pvpre[1] )
		  pv_vars <- which( substring( colnames(data) , 1 , nc1 ) == pvpre[1] )
		  pv_vars <- gsub( pvpre[1] , "" , colnames(data)[ pv_vars ] )
#		  cdata <- FALSE
		datarep <- data[ , repvars ]
        RR <- ncol(datarep)		
		fayfac <- 1 /  RR / ( 1 - .5)^2
        data <- data[ , - repvars ]	

				}								
	#******** generate replicate weights
	if ( jktype != "RW_PISA") {
	    #**** bug fix ARb 2014-12-11
		# redefine jackknife zones
		jkzones1 <- unique( data[,jkzone] )
		data[,jkzone] <- match( data[,jkzone] , jkzones1)
		#***********
		RR <- max( data[,jkzone] )		
		prblen <- 10
		prbar <- BIFIE.progressbar( ops = RR , prblen = prblen )
		cat("+++ Generate replicate weights\n")
		cat(paste0("|" , paste0(rep("*",prblen), collapse="") , "|\n|")) ; utils::flush.console()	
		addname <- 10^( floor( log( RR+.5 , 10 ) )  + 1 )
		data[ , jkzone ] <- match( data[ , jkzone ] , unique( data[ , jkzone] ) )		
		datarep <- .Call( "bifie_jack_timss" , wgt_= data[,wgt] , data[,jkzone]-1 , data[,jkrep] , RR , jkfac ,
						prbar , PACKAGE="BIFIEsurvey" )
		colnames(datarep) <- paste0("w_fstr" , substring( paste0(addname +1:RR),2) )		
		cat("|\n")
					}

	#******** generate replicated datasets for datasets
	if ( is.null( pv_vars) ){ 
				datalist <- dataL  
							}						
	if ( ! is.null( pv_vars )){
		dfr <- NULL
		VV <- length(pv_vars)
		for (vv in 1:VV){
			vv1 <- pv_vars[vv]
			if (jktype != "RW_PISA"){ 
				ind.vv1 <- which( substring( colnames(data) , 1 , nchar( vv1 ) ) == pv_vars[vv] )
									}  else {
				varsel <- paste0( pvpre , vv1	)
				ind.vv1 <- which( colnames(data) %in% varsel )									
									}
									
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
						}  # end imputations						
					}  # end pv_vars
		
	if ( ! is.null(fayfac0) ){
		fayfac <- fayfac0	
					}
					
	#*** create BIFIE.data object
	
	bifiedat <- BIFIE.data( datalist , wgt = data[, wgt ] , wgtrep = datarep , fayfac = fayfac ,
							cdata=cdata , NMI = FALSE )
	bifiedat$CALL <- cl
	return(bifiedat)
}
###############################################################################