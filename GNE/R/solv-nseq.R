#############################################################################
#   Copyright (c) 2012 Christophe Dutang                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################
### non-smooth equations in GNE
###
###         R functions
### 

nseq <- function(xinit, Phi, jacPhi, argfun, argjac, 
	method=c("Newton", "Broyden", "Levenberg-Marquardt"), 
	global=c("gline", "qline", "dbldog", "pwldog", "none"), 
	control=list(), silent=TRUE, ...)	
{
#	print(method)
	
	method <- match.arg(method, c("Newton", "Broyden", "Levenberg-Marquardt"))
	global <- match.arg(global, c("gline", "qline", "dbldog", "pwldog", "none"))
	
	#default control parameters
	con <- list(ftol=1e-6, maxit=100, trace=0)
	namc <- names(con)
	con[namc <- names(control)] <- control

	test.try <- try( Phi(xinit, argfun, argjac), silent=silent )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evalate Phi(init).", fvec=NA) )
	if(any(is.na(test.try) || is.nan(test.try) || is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Phi(init) has infinite, NA or NaN values.", fvec=NA) )
	
	test.try <- try( jacPhi(xinit, argfun, argjac), silent=silent )
	if(class(test.try) == "try-error")
	return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evalate jacPhi(init).", fvec=NA) )
	if(any(is.na(test.try) || is.nan(test.try) || is.infinite(test.try)) )
	return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="jacPhi(init) has infinite, NA or NaN values.", fvec=NA) )
	
	if(method != "Levenberg-Marquardt")
	test.try <- try( nleqslv(xinit, Phi, jacPhi, argfun, argjac,
		method = method, global = global, control=con, ...), silent=silent)
	
	if(method == "Levenberg-Marquardt")
	{
		LM.param <- match.arg(con$LM.param, c("merit", "jacmerit", "min", "adaptive"))
		
		if(LM.param == "adaptive")
			test.try <- try( nseq.LM.adapt(xinit, Phi, jacPhi, argfun=argfun, argjac=argjac, 
							global=global, control=con), silent=silent)	

		if(LM.param != "adaptive")
			test.try <- try( nseq.LM(xinit, Phi, jacPhi, argfun=argfun, argjac=argjac, 
							global=global, control=con), silent=silent)
	}	
	
	
	if(class(test.try) == "try-error")
		res <- list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					message=paste("Error in the non smooth problem:", test.try, "."), 
					fvec=NA)
	if(class(test.try) != "try-error")
		res <- list(par = test.try$x, value = sqrt(sum( test.try$fvec^2 )), 
					counts = c(phicnt = test.try$nfcnt, jaccnt = test.try$njcnt), 
					iter = test.try$njcnt, code = test.try$termcd, 
					message = test.try$message, fvec= test.try$fvec)
	res
}	


nseq.LM <- function(xinit, Phi, jacPhi, argfun, argjac, control, global, silent=TRUE)	
{	
	global <- match.arg(global, c("gline", "qline", "none"))
	control$LM.param <- match.arg(control$LM.param, c("merit", "jacmerit", "min"))
	
	#default control parameters
	con <- list(ftol=1e-6, xtol=1e-5, btol=1e-2, maxit=100, trace=0, sigma=1/2, 
				echofile=NULL, delta=1)
	namc <- names(con)
	con[namc <- names(control)] <- control
	

	
	xk_1 <- xinit
	xk <- xinit
	fk <- Phi(xinit, argfun=argfun)
	Jacfk <- jacPhi(xinit, argjac=argjac)
	iter <- 0
	inner.iter <- 0
	nfcnt <- 0
	njcnt <- 0
	termcd <- 0
	
#traces in R console
	if(con$trace >= 1 && is.null(con$echofile))
	cat("**** k", iter, "\n")
	if(con$trace >= 1 && is.null(con$echofile))
	cat("x_k", xk, "\n")	

	while(termcd == 0 && iter < con$maxit) 
	{
		b <- -crossprod( Jacfk, fk )
		
		if(con$LM.param == "merit")
			lambdak <- sqrt( sum( fk^2 ) )^con$delta
		if(con$LM.param == "jacmerit")
			lambdak <- sqrt( sum( b^2 ) )^con$delta
		if(con$LM.param == "min")
			lambdak <- min( sqrt( sum( fk^2 ) ), sqrt( sum( b^2 ) ) )^con$delta
		
		A <- crossprod( Jacfk, Jacfk ) + lambdak * diag( length(xk) )

		mycatch <- try( dk <- qr.solve(A, b) , silent=silent)
		
		if(class(mycatch) == "try-error")
		{
			termcd <- 5
			if( length(strsplit(as.character(mycatch), "singul")[[1]]) == 2 )
				termcd <- 6
			break
		}
		
		if(global == "none")
		{
			xkp1 <- xk + dk
			inner.iter <- 0
		}
		
		if(global != "none")
		{
			normfk <- crossprod(fk)/2
			
			inner.iter <- 0
			
			slopek <- crossprod(Jacfk %*% dk, fk)
			
			minstep <- con$xtol / max( abs(dk) / pmax(xk, 1) )
			
			stepk <- 1
			
			if(global == "gline")
			while(stepk > minstep)
			{
				normfp <- crossprod(Phi(xk + stepk*dk, argfun=argfun))/2
				
				#traces in R console	
				if(con$trace >= 3)
				cat("i", inner.iter, "\tstep", stepk, "\tleft term", normfp, "\tright term\t", normfk + con$btol * stepk * slopek, "\n")			
				
				#cat("largest\t", max(abs(Phi(xk + stepk*dk, argfun=argfun))), "\n")	
				
				#check Armijo condition
				if(normfp <= normfk + con$btol * stepk * slopek)
				{
					break
				}
				
				inner.iter <- inner.iter + 1	
				stepk <- con$sigma * stepk
			}
			
			if(global == "qline")
			while(stepk > minstep)
			{
				normfp <- crossprod(Phi(xk + stepk*dk, argfun=argfun))/2
				
				#traces in R console	
				if(con$trace >= 3)
					cat("i", inner.iter, "\tstep", stepk, "\tleft term", normfp, "\tright term\t", normfk + con$btol * stepk * slopek, "\n")			
				
				
				#check Armijo condition
				if(normfp <= normfk + con$btol * stepk * slopek)
					break
				
				inner.iter <- inner.iter + 1	
				stepk <- - as.numeric( (stepk)^2 * slopek / 2 / (normfp -  normfk - slopek)	)	

			}
			#traces in R console	
			if(con$trace >= 3)
				cat("\n")	
			
			if(stepk <= minstep)
			{	
				termcd <- 3
				break
			}else if(normfp <= normfk + con$btol * stepk * slopek) 
			{	
				xkp1 <- xk + stepk*dk
			}else
			{
				stop("internal error in nseq.LM function.")
			}
		}
		

		
		xk_1 <- xk
		xk <- xkp1
		
		fk_1 <- fk	
		fk <- Phi(xk, argfun=argfun)
		Jacfk <- jacPhi(xk, argjac=argjac)
		
		if(any(is.nan(Jacfk)) || any(is.nan(fk)))
			termcd <- -10
			
		
		iter <- iter+1
		nfcnt <- nfcnt + 1 + inner.iter 
		njcnt <- njcnt + 1

		#traces in R console
		if(con$trace >= 1 && is.null(con$echofile))
			cat("**** k", iter, "\n")		
		if(con$trace >= 1 && is.null(con$echofile))
			cat("x_k", xk, "\n")
		if(con$trace >= 2 && is.null(con$echofile))
			cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n", "lambdak", lambdak, "\n")
		
		
		
		#termination criterion, see Schnabel algo A7.2.3
		if(max( abs( fk ) ) <= con$ftol)
			termcd <- 1
		if(iter >= con$maxit)
			termcd <- 4
		if(max( abs(xk - xk_1) / abs(xk) ) <= con$xtol)
			termcd <- 2
		
	}
	
	message <- NA
	if(termcd == 1)
		message <- "Function criterion near zero"
	if(termcd == 2)
		message <- "x-values within tolerance `xtol'"
	if(termcd == 3)
		message <- "No better point found (algorithm has stalled)"
	if(termcd == 4)
		message <- "Iteration limit exceeded"
	if(termcd == 5)
		message <- "Jacobian is too ill-conditioned"
	if(termcd == 6)
		message <- "Jacobian is singular"
	if(termcd == -10)
		message <- "Analytical Jacobian most likely incorrect"
	
    	
	
	list(x= as.vector(xk), fvec=fk, nfcnt=nfcnt, njcnt=njcnt, iter=iter, 
		 termcd=termcd, message=message)
	
}





nseq.LM.adapt <- function(xinit, Phi, jacPhi, argfun, argjac, control, global, silent=TRUE)	
{	
	global <- match.arg(global, c("none"))
	
	#default control parameters
	con <- list(ftol=1e-6, xtol=1e-5, maxit=100, trace=0, p0=1e-4, p1=1/4, p2=3/4, 
				echofile=NULL, mumin=1e-8, muinit=1e-2)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	
	
	xk_1 <- xinit
	xk <- xinit
	fk <- Phi(xinit, argfun=argfun)
	Jacfk <- jacPhi(xinit, argjac=argjac)
	muk <- con$muinit
	iter <- 0
	inner.iter <- 0
	nfcnt <- 0
	njcnt <- 0
	termcd <- 0
	
#traces in R console
	if(con$trace >= 1 && is.null(con$echofile))
		cat("**** k", iter, "\n")
	if(con$trace >= 1 && is.null(con$echofile))
		cat("x_k", xk, "\n")	
	if(con$trace >= 2 && is.null(con$echofile))
		cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n", 
			"lambdak", muk * sqrt( sum( fk^2 ) ), "muk", muk, "\n")
	
	
	while(termcd == 0 && iter < con$maxit) 
	{
		normfk <- sqrt( sum( fk^2 ) )
		A <- crossprod( Jacfk, Jacfk ) + muk * normfk * diag( length(xk) )
		b <- -crossprod( Jacfk, fk )
		
		mycatch <- try( dk <- qr.solve(A, b) , silent=silent )
		
		if(class(mycatch) == "try-error")
		{
			termcd <- 5
			if( length(strsplit(as.character(mycatch), "singul")[[1]]) == 2 )
				termcd <- 6
			break
		}
		
		rhok <- normfk^2 - sum(Phi(xk + dk, argfun=argfun)^2)
		rhok <- rhok / ( normfk^2 - sum((fk + Jacfk%*%dk)^2) )
		
		#sufficient decrease
		if(rhok > con$p0)
			xkp1 <- xk + dk
		else
			xkp1 <- xk 
		
		#update LM param
		if(rhok < con$p1)
			muk <- muk * 2
		if(rhok > con$p2)
			muk <- max(muk / 2, con$mumin)
		
		xk_1 <- xk
		xk <- xkp1
		
		fk_1 <- fk	
		fk <- Phi(xk, argfun=argfun)
		Jacfk <- jacPhi(xk, argjac=argjac)
		
		if(any(is.nan(Jacfk)) || any(is.nan(fk)))
		termcd <- -10
		
		
		iter <- iter+1
		nfcnt <- nfcnt + 2 
		njcnt <- njcnt + 1
		
		#traces in R console
		if(con$trace >= 1 && is.null(con$echofile))
			cat("**** k", iter, "\n")
		if(con$trace >= 1 && is.null(con$echofile))
			cat("x_k", xk, "\n")
		if(con$trace >= 2 && is.null(con$echofile))
			cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n", 
				"lambdak", muk * normfk, "muk", muk, "rhok", rhok, "\n")
		
		
		
		#termination criterion, see Schnabel algo A7.2.3
		if(max( abs( fk ) ) <= con$ftol)
			termcd <- 1
		if(iter >= con$maxit)
			termcd <- 4
		#NB: check x_k - x_k_1 within tolerance only if dk is accepted
		if(max( abs(xk - xk_1) / abs(xk) ) <= con$xtol && rhok > con$p0)
			termcd <- 2
		
	}
	
	message <- NA
	if(termcd == 1)
		message <- "Function criterion near zero"
	if(termcd == 2)
		message <- "x-values within tolerance `xtol'"
	if(termcd == 3)
		message <- "No better point found (algorithm has stalled)"
	if(termcd == 4)
		message <- "Iteration limit exceeded"
	if(termcd == 5)
		message <- "Jacobian is too ill-conditioned"
	if(termcd == 6)
		message <- "Jacobian is singular"
	if(termcd == -10)
		message <- "Analytical Jacobian most likely incorrect"
	
	
	
	list(x= as.vector(xk), fvec=fk, nfcnt=nfcnt, njcnt=njcnt, iter=iter, 
		 termcd=termcd, message=message)
	
}
