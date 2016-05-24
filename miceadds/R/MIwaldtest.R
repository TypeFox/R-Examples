
########################################
# MI Wald test
MIwaldtest <- function( qhat , u , Cdes=NULL , rdes=NULL ,
				testnull=NULL){
		
		
		# conversion of inputs
		if ( class(qhat) %in% c("array","data.frame","matrix") ){
			qhat <- qhat2list_MI(qhat)
						}				
		if ( class(u) %in% c("array") ){
			u <- u2list_MI(u)
						}
						
		if ( ! is.null(testnull) ){
			k <- length(testnull)
			pars <- names( qhat[[1]][[1]] )
			des <- create.designMatrices.waldtest( pars = pars , k=k)
			Cdes <- des$Cdes
			rdes <- des$rdes
			for (ii in 1:k){
				Cdes[ ii , testnull[ii] ] <- 1
							}
						}


		#**************************************
		# compute distribution of linear form
		NB <- length( qhat )
		NV <- length( qhat[[1]] )
		NW <- 1
		
		# qhat and u for linear forms
		qhat0 <- qhat
		u0 <- u
		for (bb in 1:NB){
				u00 <- u0[[bb]]
				q0 <- as.vector(qhat0[[bb]])	
				qhat[[bb]] <- ( Cdes %*% q0 - rdes )[,1]
				u[[bb]] <- Cdes %*% u00 %*% t(Cdes)
						}
		# res0 <- mitools::MIcombine( results=qhat , variances=u )
		#*******************
		# compute F test (D1 statistic)
		u1 <- qhat1 <- as.list( 1:NB )
		for (bb in 1:NB){
				qv1 <- list(1)
				qv1[[1]] <- qhat[[bb]] 
				qhat1[[bb]] <- qv1
				qv1[[1]] <- u[[bb]] 
				u1[[bb]] <- qv1
						}

		res0 <- NMIcombine( qhat1 , u1 )			
		# linear_hyp <- res1
		eps <- 1E-20
		ubar <- res0$ubar
		qbar <- res0$qbar
		Bm <- res0$Bm
		Wm <- res0$Wm
		k <- nrow(Cdes)
		# quadratic form
		uinv <- solve(ubar)		
		rmb <- (1+1/NB)*sum(diag( Bm %*% uinv )) / k
		rmw <- (1-1/NW)*sum(diag( Wm %*% uinv )) / k
		stat <- t(qbar) %*% uinv %*% qbar
		stat <- stat / ( k * ( 1 + rmb + rmw ) )
		stat <- stat[1,1]
		df1 <- k	
		df2 <- rmb^2 / ( NB - 1 + eps ) / ( 1 + rmw + rmb )^2 +
				   rmw^2 / ( NB*( NW - 1 + eps) ) / ( 1 + rmw + rmb )^2
		df2 <- k / df2 	

		stat <- data.frame( "F" = stat , "df1" = df1 ,
					          "df2" = df2 , 
							  "pval" = 1 - stats::pf( stat , df1=df1 , df2 = df2 ) )			
		res <- list( stat=stat , linear_hyp = res0 ,
					qhat = qhat , u=u , Cdes = Cdes , rdes=rdes )
		class(res) <- "MIwaldtest"
		return(res)										
						}

##############################################################
summary.MIwaldtest <- function(object, digits=4 ,...){
	obji <- object$stat
	V <- ncol(obji)
	for (vv in 1:V){
		obji[,vv] <- round( obji[,vv] , digits )
					}
	cat("Wald Test\n")
	print(obji)
	cat("\nLinear Hypotheses\n")
	summaryMIwaldtest_linear_hyp(object$linear_hyp , digits=digits)
			}
###############################################################


						
################################################################
# convert qhat into a list
qhat2list_MI <- function( qhat ){
			qhat0 <- qhat
			dq <- dim(qhat)
			NB <- dq[1]
			NV <- dq[2]
			qhat <- as.list(1:NB)	
			for (bb in 1:NB){
  				qhat[[bb]] <- qhat0[ bb , ]
							}
			return(qhat)
				}
################################################################
# convert u into a list
u2list_MI <- function( u ){
			u0 <- u
			dq <- dim(u)
			NB <- dq[1]
			NV <- dq[2]
			u <- as.list(1:NB)
			for (bb in 1:NB){
					u[[bb]] <- u0[ bb ,  ,]
							}
			return(u)
				}
################################################################


#################################################################
summaryMIwaldtest_linear_hyp <- function(object, digits) {
    x <- object
    table <- array( x$qbar, dim = c(length(x$qbar), 10) )
    dimnames(table) <- list(labels(x$qbar), 
			c("est", "se", "t", "df", "Pr(>|t|)", "lo 95", "hi 95", 
					"fmi" , "fmi_Betw" , "fmi_Within"))
    table[, 2] <- sqrt( diag(x$Tm) )
    table[, 3] <- table[, 1]/table[, 2]
    table[, 4] <- x$df
    table[, 5] <- if (all(x$df > 0)) 
        2 * (1 - stats::pt(abs(table[, 3]), x$df)) else NA
    table[, 6] <- table[, 1] - qt(0.975, x$df) * table[, 2]
    table[, 7] <- table[, 1] + qt(0.975, x$df) * table[, 2]
#    if (is.null(x$nmis) | is.null(names(x$qbar)))
#        table[, 8] <- NA else table[, 8] <- x$nmis[names(x$qbar)]
    table[, "fmi"] <- x$lambda
    table[, "fmi_Betw"] <- x$lambda_Between
	table[, "fmi_Within"] <- x$lambda_Within
	table <- as.data.frame(table)
	if ( is.na(table$se)[1] ){
			table$df <- NA
						}
						
	table$fmi_Betw <- NULL
	table$fmi_Within <- NULL						
	table0 <- table
	for (vv in seq(1 , ncol(table) ) ){
		table[,vv] <- round( table[,vv] , digits=digits )
							}
    print(table)
	invisible(table0)
}
#################################################################						