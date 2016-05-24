################################################
# d effect size for missingness indicators
mi_dstat <- function(dat){
	resp <- is.na(dat)
	# means of missing data
	miss_vars <- colnames(resp)[ colMeans( resp ) > 0 ]
	MV <- length(miss_vars)
	V <- ncol(dat)

	dstat <- matrix( 0 , nrow=MV , ncol=V )
	rownames(dstat) <- miss_vars
	colnames(dstat) <- colnames(dat)

	for (vv in 1:MV){
		# vv <- 5
		dat_vv0 <- dat[  resp[ , miss_vars[vv]  ] , , drop=FALSE ]
		dat_vv1 <- dat[ ! resp[ , miss_vars[vv]  ] , , drop=FALSE ]    
		m0 <- colMeans( dat_vv0 , na.rm=TRUE )
		m1 <- colMeans( dat_vv1 , na.rm=TRUE )    
		sd0 <- apply( dat_vv0 , 2 , sd , na.rm=TRUE)
		sd1 <- apply( dat_vv1 , 2 , sd , na.rm=TRUE)    
		d <- (m0-m1) / sqrt( ( sd0^2 + sd1^2 ) / 2 )
		dstat[vv,] <- d
			}
	return(dstat)
			}
#####################################################			