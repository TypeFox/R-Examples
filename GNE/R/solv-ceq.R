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
### constrained equations in GNE
###
###         R functions
### 

ceq <- function(xinit, dimx, dimlam, 
	Hfinal, jacHfinal, argfun, argjac,  
	method=c("PR", "AS"), global=c("gline", "qline", "pwldog", "none"), 
	xscalm=c("fixed", "auto"), 
	control=list(), silent=TRUE, ...)	
{
	#currently, the analytical Jacobian is mandatory
	
	method <- match.arg(method, c("PR", "AS"))
	global <- match.arg(global, c("gline", "qline", "pwldog", "none"))
	xscalm <- match.arg(xscalm, c("fixed", "auto"))
	if(method == "PR" && global == "pwldog")
		stop("potential reduction algo with trust region not yet implemented.")
	if(method == "AS" && global != "pwldog")
		stop("affine scaling algo with line search not yet implemented.")

	#default control parameters
	con <- list(ftol=1e-6, maxit=100, trace=0)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	if(method == "PR")
	test.try <- try( ceq.PR(xinit, dimx, dimlam, Hfinal, jacHfinal, 
		argfun, argjac, merit=psi.ce, gradmerit=gradpsi.ce, checkint=checkint.ce, 
		control= control, global=global, silent=silent), silent=silent)		
	if(method == "AS")
	test.try <- try( ceq.AS(xinit, dimx, dimlam, Hfinal, jacHfinal, 
		argfun, argjac, global=global, xscalm=xscalm, 
		control= control, silent=silent), silent=silent)		
		
	if(class(test.try) == "try-error")
		res <- list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					message=paste("Error in the non smooth problem:", test.try, "."), 
					fvec=NA)
	if(class(test.try) != "try-error")
		res <- list(par = test.try$x, value = sqrt(sum( test.try$fvec^2 )), 
					counts = c(phicnt = test.try$nfcnt, jaccnt = test.try$njcnt), 
					iter = test.try$njcnt, code = test.try$termcd, 
					message = test.try$message, fvec= test.try$fvec,
					xscalm = test.try$xscalm)
	res
}	


ceq.PR <- function(xinit, dimx, dimlam, Hfinal, jacHfinal, argfun, argjac,
	merit, gradmerit, checkint, control, global, silent=TRUE)	
{	
	global <- match.arg(global, c("gline", "qline", "none"))
	
	#default control parameters
	con <- list(ftol=1e-6, xtol=1e-5, btol=1e-2, maxit=100, trace=0, sigma=1/2, 
				echofile=NULL, delta=1, forcingpar=.1, zeta=length(xinit)/2)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	xk_1 <- xinit
	xk <- xinit
	fk <- Hfinal(xinit, argfun=argfun)
	sigmak <- con$forcingpar
	
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
		cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n")
	

	while(termcd == 0 && iter < con$maxit) 
	{
		
		dk <- ceq.PR.direction(xk, dimx, dimlam, Hfinal, jacHfinal, 
			argfun, argjac, sigmak, silent=silent)
		
		if(class(dk) == "try-error")
		{
			termcd <- 5
			if( length(strsplit(as.character(dk), "singul")[[1]]) == 2 )
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
			inner.iter <- 0
			minstep <- con$xtol / max( abs(dk) / pmax(xk, 1) )
			
			slopek <- crossprod( gradmerit(xk, dimx, dimlam,  
				Hfinal, jacHfinal, argfun, argjac, con$zeta), dk)
			stepk <- 1			
			
			if(global == "gline")
				LSres <- try( linesearch.geom.cond(xk, dk, slopek, con, merit= merit, 
					checkint=checkint, dimx=dimx, dimlam=dimlam, 
					Hfinal=Hfinal, argfun=argfun, zeta=con$zeta) )
			
			if(global == "qline")
				LSres <- try( linesearch.quad.cond(xk, dk, slopek, con, merit= merit, 
					checkint=checkint, dimx=dimx, dimlam=dimlam, 
					Hfinal=Hfinal, argfun=argfun, zeta=con$zeta) )

			if(class(LSres) == "try-error")
				stop("internal error in line search function.")
			
			if(LSres$stepk <= minstep)
			{	
				termcd <- 3
				break
			}else if(is.infinite(LSres$normfp))
			{
				termcd <- 7
				break
			}else if(LSres$normfp <= LSres$normfk + con$btol * LSres$stepk * slopek) 
			{	
				xkp1 <- xk + LSres$stepk*dk
			}else
			{
				stop("internal error in ceq.IP function.")
			}
		}
		
		xk_1 <- xk
		xk <- xkp1
		
		fk_1 <- fk	
		fk <- Hfinal(xk, argfun=argfun)		

		if(any(is.nan(fk)))
			termcd <- -10
			
		iter <- iter+1
		nfcnt <- nfcnt + 1 + ifelse(global != "none", LSres$inner.iter, 0)
		njcnt <- njcnt + 1

		#traces in R console
		if(con$trace >= 1 && is.null(con$echofile))
			cat("**** k", iter, "\n")		
		if(con$trace >= 1 && is.null(con$echofile))
			cat("x_k", xk, "\n")
		if(con$trace >= 2 && is.null(con$echofile))
			cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n")
		
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
	if(termcd == 7)
		message <- "No better point found in the interior of the domain"
	if(termcd == -10)
		message <- "Analytical Jacobian most likely incorrect"
	
    	
	
	list(x= as.vector(xk), fvec=fk, nfcnt=nfcnt, njcnt=njcnt, iter=iter, 
		 termcd=termcd, message=message, xscalm="fixed")
	
}


#z = (x, lam, w)
#with size (n, m, m)
ceq.PR.direction <- function(z, dimx, dimlam, Hfinal, jacHfinal, 
	argfun, argjac, sigma, silent=TRUE)
{
	n <- sum(dimx)
	m <- sum(dimlam)
	
	x <- z[1:n]
	lam <- z[(n+1):(n+m)]
	w <- z[(n+m+1):(n+2*m)]
		
	Hz <- Hfinal(z, argfun=argfun)
	Jz <- jacHfinal(z, argjac=argjac)
	
	a <- c( rep(0, n), rep(1, 2*m) )
	
	btot <- sigma * sum(a*Hz) / sum(a^2) * a - Hz
	bx <- btot[1:n]
	blam <- btot[(n+1):(n+m)]
	bw <- btot[(n+m+1):(n+2*m)]
	
	if(m > 0)
	{
		Ex <- Jz[1:n, (n+1):(n+m)] 
		Jacgx <- Jz[(n+1):(n+m), 1:n]
		
		mat4x <- Jz[1:n, 1:n] + Ex %*% diag(lam/w) %*% Jacgx
		vec4x <- bx - Ex %*% diag(1/w) %*% bw + Ex %*% diag(lam/w) %*% blam
		
		d4x <- try( qr.solve(mat4x, vec4x), silent=silent)
		if(class(d4x) != "try-error")
		{
			d4w <- blam - Jacgx %*% d4x
			d4lam <- diag(1/w) %*% (bw - diag(lam) %*% d4w)
			return(c(d4x, d4lam, d4w))
		}else
		{
			d4x.LU <- try( solve(mat4x, vec4x), silent=silent)
			if(class(d4x.LU) != "try-error")
			{
				d4w <- blam - Jacgx %*% d4x.LU
				d4lam <- diag(1/w) %*% (bw - diag(lam) %*% d4w)
				return(c(d4x.LU, d4lam, d4w))
			}			
	#		print(mat4x, digits=15)
	#		print(vec4x, digits=15)
	#		print(d4x)
	#		cat("\n\n")
			return(d4x)	
		}
	}else
	{
		mat4x <- Jz
		vec4x <- -Hz
		
		d4x <- try( qr.solve(mat4x, vec4x), silent=silent)
		if(class(d4x) != "try-error")
			return(d4x)
		else
		{
			d4x.LU <- try( solve(mat4x, vec4x), silent=silent)
			if(class(d4x.LU) != "try-error")
				return(d4x.LU)
			else
				return(d4x)
		}
	}
}

ceq.AS <- function(xinit, dimx, dimlam, Hfinal, jacHfinal, argfun, argjac,
	global, xscalm, control, silent=TRUE)	
{
	global <- match.arg(global, c("pwldog"))
	
	#default control parameters
	con <- list(ftol=1e-6, xtol=1e-5, maxit=100, trace=0, theta=0.99995, 
				echofile=NULL, radiusmin=1, reducmin=0.1, radiusmax=1e10, 
				radiusred=1/2, reducred=1/4, radiusexp=2, reducexp=3/4)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	n <- sum(dimx)
	m <- sum(dimlam)
	lower <- c(rep(-Inf, n), rep(0, 2*m))
	
	if(any(xinit < lower))
		stop("wrong initial point.")
	
	xk_1 <- xinit
	xk <- xinit
	fk <- Hfinal(xinit, argfun=argfun)
	Jfk <- jacHfinal(xinit, argjac=argjac)
	Jfkfk <- t(Jfk) %*% fk
#		c(t(Jfk[1:n, 1:n]) %*% fk[1:n] + t(Jfk[1:m+n, 1:n]) %*% fk[1:m+n],
#		   t(Jfk[1:n, 1:m+n]) %*% fk[1:n]
	
	iter <- 0
	inner.iter <- 0
	nfcnt <- 0
	njcnt <- 0
	termcd <- 0
	if(xscalm == "auto")
		scalvec <- scaling.AS(xinit, lower, n, m, Jfkfk)
	else
		scalvec <- rep(1, n+2*m) 
	deltak <- 1

	#traces in R console
	if(con$trace >= 1 && is.null(con$echofile))
		cat("**** k", iter, "\n")
	if(con$trace >= 1 && is.null(con$echofile))
		cat("x_k", xk, "\n")	
	if(con$trace >= 2 && is.null(con$echofile))
		cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n")
	
	while(termcd == 0 && iter < con$maxit) 
	{
		
		#Newton point
		pn <- try( qr.solve(Jfk, -fk), silent=silent)
		if(class(pn) == "try-error")
		{
			termcd <- 5
			if( length(strsplit(as.character(pn), "singul")[[1]]) == 2 )
			termcd <- 6
			break
		}
		
		#rescaled gradient of the merit function
		rgk <- Jfkfk / scalvec^2
		#Cauchy point (steepest descent) in the scaled setting
		pc <- - sum((Jfkfk/scalvec)^2) / sum( (Jfk %*% rgk)^2 ) * rgk		
#		if(sqrt(sum(pc^2)) > deltak) #truncated Cauchy point not needed
#			pc <- - deltak / sqrt(sum((Jfkfk/scalvec)^2)) * rgk
		
		ndir <- truncpoweldogleg(n+2*m, pn, sqrt(sum(pn^2)), pc, 
						sqrt(sum(pc^2)), deltak, lower, rep(Inf, n+2*m), 
						scalvec, xk, con)
		rhok <- (crossprod(Jfkfk, ndir$nextp) + crossprod(Jfk %*% ndir$nextp)/2) / (crossprod(Jfkfk, pc) + crossprod(Jfk %*% pc)/2)
		if(rhok < con$reducmin)
			ndir <- list(nextp=pc, type="Cauchy step")
		
		
		#compute the actual reduction
		fkpp <- Hfinal(xk + ndir$nextp, argfun=argfun)
		rhok <- (sum(fkpp^2) - sum(fk^2)) / (crossprod(Jfkfk, ndir$nextp) + crossprod(Jfk %*% ndir$nextp)/2)
		
		if(rhok > con$reducmin)
		{	
			xkp1 <- xk + ndir$nextp
			fkp1 <- fkpp
		}else
		{
			xkp1 <- xk
			fkp1 <- fk
		}
		
		#update radius
		if(rhok < con$reducred)
			deltak <- max(deltak * con$radiusred, con$radiusmin)
		else if(rhok > con$reducexp)
			deltak <- min(deltak * con$radiusexp, con$radiusmax)
		
		#prepare for next iteration
		xk_1 <- xk
		xk <- xkp1
		fk_1 <- fk	
		fk <- fkp1
		
		#Jacobian and scaling
		Jfk <- jacHfinal(xk, argjac=argjac)
		Jfkfk <- t(Jfk) %*% fk
		if(xscalm == "auto")
			scalvec <- scaling.AS(xk, lower, n, m, Jfkfk)
		else
			scalvec <- rep(1, n+2*m) 
		
		if(any(is.nan(fk)))
			termcd <- -10
		
		iter <- iter+1
		nfcnt <- nfcnt + 1
		njcnt <- njcnt + 1
			
		
		#traces in R console
		if(con$trace >= 1 && is.null(con$echofile))
			cat("**** k", iter, "\n")		
		if(con$trace >= 1 && is.null(con$echofile))
			cat("x_k", xk, "\n")
		if(con$trace >= 2 && is.null(con$echofile))
			cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n")
		if(con$trace >= 3 && is.null(con$echofile))
			cat(" dk -", ndir$type, ndir$nextp,  "\n")
		if(con$trace >= 3 && is.null(con$echofile) && xscalm == "auto")
			cat(" scaling", scalvec,  "\n\n")
			
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
		 termcd=termcd, message=message, xscalm=xscalm)
	
}

#lower bound only
scaling.AS <- function(x, l, n, m, Jff)
{
	xscalv <- rep(1, n+2*m)
	xscalv[Jff >= 0 & is.finite(l)] <- (x - l)[Jff >= 0 & is.finite(l)]
	1/sqrt(abs(xscalv))
}

