

###############################################################
## estimating F entries
.noharm.estF <- function( Fval , Pval , Fpatt , Ppatt , 
		I , D ,  b0.jk , b1.jk , b2.jk , b3.jk , wgtm , pm , 
		Psival , Psipatt ){

	# compute dj
	dj <- sqrt( diag( Fval %*% Pval %*% t(Fval) ) )
	
	# compute ej
	 ej <- sqrt( 1 + dj^2 )
	# ej <- sqrt( 1+ .5^2 )
	 ej.ek <- 1 / outer( ej , ej )
	 #ej.ek <- 1+0*ej.ek
	 diag( ej.ek ) <- 0
	v0.jk <- b0.jk
	v1.jk <- b1.jk   * ej.ek
	v2.jk <- b2.jk   * ej.ek^2
	v3.jk <- b3.jk   * ej.ek^3

	
	# compute gamma.jk = f_j' P f_k
	gamma.jk <- Fval %*% Pval %*% t(Fval ) + Psival
	gamma.jk2 <- gamma.jk^2
	gamma.jk3 <- gamma.jk^3
	# gamma.jk <- gamma.jk * ej.ek
	# compute p_d ' f_k
	pd.fk <- Fval %*% Pval
	Fval_old <- Fval
#	Pval_old <- Pval


	###########################################
	# estimation F


	# derivative with respect to f_{jd}
	for (jj in 1:I){
	for (dd in 1:D){
	  if ( Fpatt[jj,dd] > 0 ){
		# dd <- 1     # dd = 1
		# jj <- 1     # Item jj
		# zeroth derivative
			# if (update){ gamma.jk <- Fval %*% Pval %*% t(Fval )	 }
#		eps0.jj <- ( wgtm[jj,] * ( pm[jj,] - v0.jk[jj,] - v1.jk[jj,]*gamma.jk[jj,] - v2.jk[jj,]*gamma.jk[jj,]^2 - 
#					v3.jk[jj,]*gamma.jk[jj,]^3 ) )					
		eps0.jj <- ( wgtm[jj,] * ( pm[jj,] - v0.jk[jj,] - v1.jk[jj,]*gamma.jk[jj,] - v2.jk[jj,]*gamma.jk2[jj,] - 
					v3.jk[jj,]*gamma.jk3[jj,] ) )					

		# first derivative eps1(jj,dd)
		eps1.jj <- - ( v1.jk[jj,]  * pd.fk[,dd] ) - 
					2 * ( v2.jk[jj,]  * gamma.jk[jj,]  * pd.fk[ , dd] ) -
					3 * ( v3.jk[jj,]  * gamma.jk2[jj,]  * pd.fk[ , dd] )
		eps2.jj <- - 2 * ( v2.jk[jj,]  * pd.fk[ , dd]^2 ) -
						6 * ( v3.jk[jj,]  * gamma.jk[jj,]  * pd.fk[ , dd]^2 )
		# total derivative
		f1.jj <- 2* eps0.jj * eps1.jj
		f2.jj <- 2* eps1.jj^2 + 2*eps0.jj * eps2.jj
		increment <-  - sum(f1.jj) / sum(f2.jj)

		increment <- ifelse( abs(increment) > .2 , .2*sign(increment) , increment )
		Fval[jj,dd] <- Fval[ jj,dd] + increment 
							}
					}
				}
	#************* end F

			
 res <- list("Fval" = Fval , "change" = max( abs( Fval - Fval_old ) )  ) 
 return(res)
    }
######################################################################





###############################################################
# estimating P entries
.noharm.estP <- function( Fval , Pval , Fpatt , Ppatt , 
		I , D ,  b0.jk , b1.jk , b2.jk , b3.jk , wgtm , pm ,
		Psival , Psipatt ){

	# compute dj
	dj <- sqrt( diag( Fval %*% Pval %*% t(Fval) ) )

	# compute ej
	 ej <- sqrt( 1 + dj^2 )
	# ej <- sqrt( 1+ .5^2 )
	 ej.ek <- 1 / outer( ej , ej )
	 #ej.ek <- 1+0*ej.ek
	 diag( ej.ek ) <- 0
	v0.jk <- b0.jk
	v1.jk <- b1.jk   * ej.ek
	v2.jk <- b2.jk   * ej.ek^2
	v3.jk <- b3.jk   * ej.ek^3
	# compute gamma.jk = f_j' P f_k
	gamma.jk <- Fval %*% Pval %*% t(Fval ) + Psival
	# gamma.jk <- gamma.jk * ej.ek
	# compute p_d ' f_k
	pd.fk <- Fval %*% Pval
	Fval_old <- Fval
	Pval_old <- Pval


	# derivative with respect to P
	for (dd in 1:D){
	for (ee in 1:D){
	  if ( ( Ppatt[dd,ee] > 0 ) & ( dd>=ee) ){
		# dd <- 1     # dd = 1
		# jj <- 1     # Item jj
		gammajk1 <- outer( Fval[ ,dd] , Fval[ ,ee] )
		if (dd==ee){ gammajk1 <- 2*gammajk1 }
		# zeroth derivative
			# if (update){ gamma.jk <- Fval %*% Pval %*% t(Fval )	 }
		eps0.jj <- ( wgtm * ( pm - v0.jk - v1.jk*gamma.jk - v2.jk*gamma.jk^2 - 
					v3.jk*gamma.jk^3 ) )
		# first derivative eps1(jj,dd)
		eps1.jj <- - ( v1.jk  * gammajk1 ) - 
					2 * ( v2.jk  * gamma.jk  * gammajk1 ) -
					3 * ( v3.jk  * gamma.jk^2  * gammajk1 )
										
		eps2.jj <- - 2 * ( v2.jk  * gammajk1^2 ) -
						6 * ( v3.jk  * gamma.jk  * gammajk1^2 )				
		# total derivative
		f1.jj <- 2* eps0.jj * eps1.jj
		f2.jj <- 2* eps1.jj^2 + 2*eps0.jj * eps2.jj
		increment <-  - sum(f1.jj) / sum(f2.jj)
		increment <- ifelse( abs(increment) > .2 , .2*sign(increment) , increment )
		Pval[dd,ee] <- Pval[ dd,ee] + increment 
		if ( dd > ee ){ Pval[ee,dd] <- Pval[dd,ee] }
							}
					}
				}
	#************* end P
			
 res <- list(  "Pval" = Pval , "change" = max( abs( Pval - Pval_old ) ) ,
		"residuals" = eps0.jj )
 return(res)
    }
######################################################################


###############################################################
## estimating Psi entries
.noharm.estPsi <- function( Fval , Pval ,  Fpatt , Ppatt , 
		I , D ,  b0.jk , b1.jk , b2.jk , b3.jk ,  wgtm , pm ,
		Psival , Psipatt ){

	# compute dj
	dj <- sqrt( diag( Fval %*% Pval %*% t(Fval) ) )

	# compute ej
	 ej <- sqrt( 1 + dj^2 )
	# ej <- sqrt( 1+ .5^2 )
	 ej.ek <- 1 / outer( ej , ej )
	 #ej.ek <- 1+0*ej.ek
	 diag( ej.ek ) <- 0
	v0.jk <- b0.jk
	v1.jk <- b1.jk   * ej.ek
	v2.jk <- b2.jk   * ej.ek^2
	v3.jk <- b3.jk   * ej.ek^3
	# compute gamma.jk = f_j' P f_k
	gamma.jk <- Fval %*% Pval %*% t(Fval ) + Psival
	# gamma.jk <- gamma.jk * ej.ek
	# compute p_d ' f_k
#	pd.fk <- Fval %*% Pval
	Psival_old <- Psival
	###########################################
	# estimation F


	# derivative with respect to f_{jd}
	for (jj in 1:(I-1)){
	for (kk in (jj+1):I){
	  if ( Psipatt[jj,kk] > 0 ){
		# dd <- 1     # dd = 1
		# jj <- 1     # Item jj
		# zeroth derivative 
			# if (update){ gamma.jk <- Fval %*% Pval %*% t(Fval )	 }
		eps0.jj <- ( wgtm[jj,kk] * ( pm[jj,kk] - v0.jk[jj,kk] - v1.jk[jj,kk]*gamma.jk[jj,kk] - 
					v2.jk[jj,kk]*gamma.jk[jj,kk]^2 - 
					v3.jk[jj,kk]*gamma.jk[jj,kk]^3 ) )					
		# first derivative eps1(jj,dd)
		eps1.jj <- - ( v1.jk[jj,kk]  * 1 ) - 
					2 * ( v2.jk[jj,kk]  * gamma.jk[jj,kk]  * 1 ) -
					3 * ( v3.jk[jj,kk]  * gamma.jk[jj,]^2  * 1 )
#		eps2.jj <- - 2 * ( v2.jk[jj,]  * pd.fk[ , dd]^2 ) -
#						6 * ( v3.jk[jj,]  * gamma.jk[jj,]  * pd.fk[ , dd]^2 )
        eps2.jj <- 0 
		# total derivative
		f1.jj <- 2* eps0.jj * eps1.jj
		f2.jj <- 2* eps1.jj^2 + 2*eps0.jj * eps2.jj
		increment <-  - sum(f1.jj) / sum(f2.jj)

		increment <- ifelse( abs(increment) > .2 , .2*sign(increment) , increment )
		Psival[jj,kk] <- Psival[ jj,kk] + increment 
		Psival[kk,jj] <- Psival[jj,kk]
							}
					}
				}
	#************* end Psi
			
 res <- list("Psival" = Psival , "change" = max( abs( Psival - Psival_old ) ))
 return(res)
    }
######################################################################





###############################################################
# estimating P entries
.noharm.est.residuals <- function( Fval , Pval , Fpatt , Ppatt , 
		I , D ,  b0.jk , b1.jk , b2.jk , b3.jk , wgtm , pm ,
		Psival , Psipatt ){

	# compute dj
	dj <- sqrt( diag( Fval %*% Pval %*% t(Fval) ) )

	# compute ej
	 ej <- sqrt( 1 + dj^2 )
	# ej <- sqrt( 1+ .5^2 )
	 ej.ek <- 1 / outer( ej , ej )
	 #ej.ek <- 1+0*ej.ek
	 diag( ej.ek ) <- 0
	v0.jk <- b0.jk
	v1.jk <- b1.jk   * ej.ek
	v2.jk <- b2.jk   * ej.ek^2
	v3.jk <- b3.jk   * ej.ek^3
	# compute gamma.jk = f_j' P f_k
	gamma.jk <- Fval %*% Pval %*% t(Fval ) + Psival
	# gamma.jk <- gamma.jk * ej.ek
	# compute p_d ' f_k
	pd.fk <- Fval %*% Pval
	Fval_old <- Fval
	Pval_old <- Pval
	eps0.jj <- ( wgtm * ( pm - v0.jk - v1.jk*gamma.jk - v2.jk*gamma.jk^2 - 
					v3.jk*gamma.jk^3 ) )
	return(eps0.jj)			
    }
######################################################################
