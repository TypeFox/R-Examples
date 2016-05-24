

##############################################################
# Wald test for nested multiply imputed datasets
NMIwaldtest <- function( qhat , u , Cdes = NULL , rdes = NULL ,
					testnull = NULL ){
		# convert qhat into a list if necessary
		if ( class(qhat) == "array" ){
			qhat <- qhat2list(qhat)
					}
		# convert u into a list			
		if ( class(u) == "array" ){
			u <- u2list(u)
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
		NW <- length( qhat[[1]] )
		NV <- length( qhat[[1]][[1]] )
		
		# qhat and u for linear forms
		qhat0 <- qhat
		u0 <- u
		for (bb in 1:NB){
		   for (ww in 1:NW){
				u00 <- u0[[bb]][[ww]]
				qhat[[bb]][[ww]] <- ( Cdes %*% qhat0[[bb]][[ww]] - rdes )[,1]
				u[[bb]][[ww]] <- Cdes %*% u00 %*% t(Cdes)
							}
						}
		
		#***********
		# statistical inference
		eps <- 1E-20
		res0 <- NMIcombine( qhat , u )
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
		class(res) <- "NMIwaldtest"
		return(res)					
					
					}
##############################################################
summary.NMIwaldtest <- function(object, digits=4 ,...){
	obji <- object$stat
	V <- ncol(obji)
	for (vv in 1:V){
		obji[,vv] <- round( obji[,vv] , digits )
					}
	cat("Wald Test\n")
	print(obji)
	cat("\nLinear Hypotheses\n")
	summary(object$linear_hyp , digits=digits)
			}
###############################################################


################################################################
# convert qhat into a list
qhat2list <- function( qhat ){
			qhat0 <- qhat
			dq <- dim(qhat)
			NB <- dq[1]
			NW <- dq[2]
			NV <- dq[3]
			qhat <- as.list(1:NB)
			for (bb in 1:NB){
				qhat1 <- as.list(1:NW)
			   for (ww in 1:NW){
					qhat1[[ww]] <- qhat0[ bb , ww , ]
							}
				qhat[[bb]] <- qhat1
						}
			return(qhat)
				}
################################################################
# convert u into a list
u2list <- function( u ){
			u0 <- u
			dq <- dim(u)
			NB <- dq[1]
			NW <- dq[2]
			NV <- dq[3]
			u <- as.list(1:NB)
			for (bb in 1:NB){
				u1 <- as.list(1:NW)
			   for (ww in 1:NW){
					u1[[ww]] <- u0[ bb , ww , ,]
							}
				u[[bb]] <- u1
						}
			return(u)
				}
################################################################