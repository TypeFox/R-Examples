
########################################################
# lpcm initial parameters
# This function is mainly copied from the
# pcmodel function from the psychotools package
lpcm_inits <- function( dat , weights , maxK , b_const , W ,
		irtmodel , normalization ){
	I <- ncol(dat)
	m <- I
	oj_vec <- sapply( 1:I , FUN = function(ii){
			seq( 0 , maxK[ii] ) } , simplify=FALSE )
	y <- dat
	ctot <- vector("list", length = m)
    for (j in seq_len(m)){
		ctot[[j]] <- as.vector(tapply(weights, factor(y[, j], levels = oj_vec[[j]]), sum))
						 }
    start <- lapply(ctot, function(x){
			- cumsum(diff.default(log(x))) } ) 
			# delta_jk = log(x_j(k-1) - log(x_jk), beta_jk = cumsum_k=1^k(delta_jk)
    start <- unlist(start)
    start[ is.na(start) ] <- 0
	par_init <- ( start - b_const ) %*% W
    par_init <- par_init[1,]
	if ( ( irtmodel=="PCM2") & ( normalization == "sum") ){ 
			par_init <- 0 * par_init
								}
    return(par_init)
			}
########################################################