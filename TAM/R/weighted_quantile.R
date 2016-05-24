
#####################################################
# weighted_quantile
weighted_quantile <- function( x , w=rep(1,length(x)) , probs=seq(0,1,.25),
			type=NULL){
	dfr <- data.frame( x , w )
	dfr <- dfr[ ! is.na(x) , ]
	dfr <- dfr[ order(dfr$x) , ]
	N <- nrow(dfr)
	weights_NULL <- if( stats::sd(w) == 0 ){ TRUE } else { FALSE }
	
	#*** reweighting
	if ( weights_NULL){
		dfr$w <- dfr$w * N / sum(dfr$w)
  					   }
	if ( ! weights_NULL){
		if ( is.null(type) ){ type <- "i/n" }		
				}
	
	#*** init vector of quantiles
	PP <- length(probs)
	res <- rep(NA,PP)
	names(res) <- paste0( 100*probs , "%")		
	dfr$w_cum <- cumsum( dfr$w )
	dfr$w_cum <- dfr$w_cum / max( dfr$w_cum)
	
	for (kk in 1:PP){
		pp <- probs[kk]
		#*** specifications according to type
		
		res1 <- quantiles_type_selection( type=type , pp=pp , N=N,
					dfr = dfr , weights_NULL = weights_NULL )
		jj <- res1$jj
		GAMMA <- res1$GAMMA	
		jj1 <- res1$jj1
		quant_pp <- (1-GAMMA)*dfr[jj,1] + GAMMA * dfr[jj1,1]	
		res[kk] <- quant_pp
					}
	return(res)
			}
#####################################################			
# selection of parameters for weighted data
quantiles_type_selection <- function( type , pp , N, dfr, weights_NULL){

    # type 4 
#	if (type==4){
#		mm <- 0
#		jj <- floor(N*pp + mm)
#		gg <- N*pp + mm - jj
#		GAMMA <- gg
#			}
    eps <- 1E-10
	mm <- NULL  
	set1 <- FALSE
	if ( ! weights_NULL ){
	    type <- -9
		a1 <- dfr$w_cum <= pp
		if ( sum(a1) > 0 ){
			ind <- which( a1 )
						} else {
			ind <- 0			
						}	
#		if ( length(ind) == 0){
#			ind <- 1
#			}
		jj <- max(ind)
		jj1 <- jj + 1
		if (jj1 > N){ jj1 <- N}
		
		if (jj %in% c(0,-Inf)){ 
			jj <- 1
			set1 <- TRUE
			jj1 <- 1
				}		
		if ( jj != jj1){		
			GAMMA0 <- ( pp - dfr[jj,"w_cum"] )/ ( eps + dfr[jj1,"w_cum"] - dfr[jj,"w_cum"] )
						} else {
			GAMMA0 <- 0			
						}
        GAMMA <- GAMMA0
			}
    # type 6 
	if (type==6){
		mm <- pp
		jj <- floor(N*pp + mm)
		gg <- N*pp + mm - jj
		GAMMA <- gg
			}
	# type 7
	if (type==7){
		mm <- 1 - pp
		jj <- floor(N*pp + mm)
		gg <- N*pp + mm - jj
		GAMMA <- gg
			}

	if ( ! set1){
		jj1 <- jj+1
				}
	if (jj1 > N){ jj1 <- N}
	if (jj==0){ jj <- 1}			
	
	res1 <- list(mm=mm, jj=jj, GAMMA=GAMMA, jj1 = jj1)
	return(res1)
			}