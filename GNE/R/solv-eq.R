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
### equation solver in GNE
###
###         R functions
### 

eqsolve <- function(xinit, f, jac,
    method=c("Newton", "Levenberg-Marquardt", "Broyden"),
	global=c("line search", "none"), control=list())
{
	method <- match.arg(method, c("Newton", "Levenberg-Marquardt","Broyden"))
	global <- match.arg(global, c("line search", "none"))

	
	#default control parameters
	con <- list(tol=1e-6, maxit=100, echo=0, sigma=1/2, btol=1e-2, 
		echofile=NULL, echograph="NULL", delta=2, initlnsrch=0, minstep=1e-5)
	namc <- names(con)
	con[namc <- names(control)] <- control
 
	
	xk <- xinit
	fk <- f(xinit)
	if(method == "Broyden")	
	{
		# Wk <- sqrt(sum(fk^2)) * diag(length(xk))
		Bk <- jac(xinit)
	}	
	if(method != "Broyden")		
		Jacfk <- jac(xinit)
	iter <- 0
	
	if(class(con$echo) == "character")
		stop("wrong echo parameter.")
	if(!is.numeric(con$tol) || !is.numeric(con$maxit) || !is.numeric(con$sigma) || !is.numeric(con$btol))	
		stop("wrong control parameters.")
		

	if(!is.null(con$echofile))
	{
		dirbefore <- getwd()
		
		dirname <- sub(":", "m", sub(":", "h", Sys.time()))
		dirname <- paste(dirname, " - example", sep="")
		dir.create(dirname, showWarnings=FALSE)
		setwd(paste("./", dirname, sep=""))
	}else
	{
		dirname <- NULL
	}
	
	#traces in R console
	if(con$echo >= 1 && is.null(con$echofile))
	{
		cat("**** k", iter, "\n")
		cat("x_k", xk, "\n")
	}	

	#traces in text file
	if(con$echo >= 2 && !is.null(con$echofile))
		cat("**** k", iter,  "\n", file=con$echofile, append=TRUE)	
	if(con$echo >= 1 && !is.null(con$echofile))
		cat("x_k", xk, "\n", file=con$echofile, append=TRUE)	
				
	
		
	
	while(sqrt( sum( fk^2 ) ) > con$tol && iter < con$maxit) 
	{
		if(method == "Newton")
			xkp1 <- try( NewtonNext(xk, fk, Jacfk, TRUE) )
		if(method == "Levenberg-Marquardt")
			xkp1 <- try( LevenMarqNext(xk, fk, Jacfk, TRUE, con$delta) )
		if(method == "Broyden")	
			xkp1 <- try( QuasiNewtonNext(xk, fk, Bk, silent=TRUE, inv=FALSE)	)
			
		if(class(xkp1) == "try-error")
			break
	
		if(con$echo >= 1 && is.null(con$echofile))
			cat("**** k", iter+1, "\n")
		if(con$echo >= 2 && !is.null(con$echofile))
			cat("**** k", iter+1,  "\n", file=con$echofile, append=TRUE)

		dk <- xkp1 - xk	

		if(con$echograph == "line")	
		{
			if(!is.null(con$echofile))
				png(paste("line search - ", iter, ".png", sep=""), 1600, 800)		

			par(mfrow=c(1, 2))
			graphlinesearch(f, xk, dk, iter)
			graphContourOneIter(f, xk, dk, iter, delta=.2)
			if(!is.null(con$echofile))			
				dev.off()		
		}		
	
		if(global == "line search")
		{
			dk <- xkp1 - xk
			normfk <- crossprod(fk)/2

			i <- con$initlnsrch
			
			if(method == "Broyden")	
				slopek <- crossprod(Bk %*% dk, fk)
		
			if(method != "Broyden")	
				slopek <- crossprod(Jacfk %*% dk, fk)
				
			

			while(i < 32)
			{
				stepk <- con$sigma^i
				normfp <- crossprod(f(xk + stepk*dk))/2
				
				#traces in R console	
				if(con$echo >= 3 && is.null(con$echofile))
					cat("i", i, "\tlambda", stepk, "\tright term\t", normfk + con$btol * stepk * slopek, "\tleft term", normfp, "\n")			
				#traces in text file	
				if(con$echo >= 3 && !is.null(con$echofile))
					cat("i", i, "\tlambda", stepk, "\txk + stepk*dk\t", xk + stepk*dk, "\n", file=con$echofile, append=TRUE)			
					
				cat("largest\t", max(abs(f(xk + stepk*dk))), "\n")	

				#check Armijo condition
				if(normfp <= normfk + con$btol * stepk * slopek)
				{
					cat("Armijo satisfied\n")
					break
				}
				
				cat("\tnorm stepk*dk\t", sqrt(sum(stepk*dk^2)), "\n\n")
				i <- i+1	
					
			}
			xkp1 <- xk + stepk*dk

		}
		xk_1 <- xk
		xk <- xkp1
		
		#control step size
		if(sqrt( sum( (xk - xk_1)^2 ) ) <= con$minstep)
		{
			cat("too small step\n")
			break
		}
			
		fk_1 <- fk	
		fk <- f(xk)
		if(method == "Broyden")	
		{
			sk <- xk - xk_1
			yk <- fk - fk_1
			# Wkyk <- Wk %*% yk
			# Wk <- Wk + (sk - Wkyk) %*% t(Wkyk)  / as.numeric(crossprod(sk, Wkyk))
			
			Bk <- Bk + (yk - Bk %*% sk) %*% t(sk) / as.numeric(crossprod(sk))
			

			if(any(is.nan(Bk)))
				break
		}	
		if(method != "Broyden")	
		{	
			Jacfk <- jac(xk)
			if(any(is.nan(Jacfk)))
				break
		}	
		
		iter <- iter+1

		#traces in R console
		if(con$echo >= 1 && is.null(con$echofile))
			cat("x_k", xk, "\n")
		if(con$echo >= 2 && is.null(con$echofile))
			cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n")
		if(con$echo >= 4 && is.null(con$echofile))
		{
			if(method != "Broyden")	
			{
				cat(" Jac Phi(x_k)\n")
				print(Jacfk)
			}	
			if(method == "Broyden")	
			{
				cat(" B_k\n")
				print(Bk)
			}
		}
		if(con$echo >= 1 && is.null(con$echofile))
			cat("\n")


		#traces in text file
		if(con$echo >= 1 && !is.null(con$echofile))
			cat("x_k", xk, "\n", file=con$echofile, append=TRUE)			
		if(con$echo >= 2 && !is.null(con$echofile))
			cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n", file=con$echofile, append=TRUE)			
		if(con$echo >= 4&& !is.null(con$echofile))
		{
			if(method != "Broyden")	
			{
				cat(" Jac Phi(x_k)\n", file=con$echofile, append=TRUE)			
				capture.output(print(Jacfk), file=con$echofile, append=TRUE)		
			}
			if(method == "Broyden")	
			{
				cat(" B_k\n", file=con$echofile, append=TRUE)			
				capture.output(print(Bk), file=con$echofile, append=TRUE)		
			}		
		}
		if(con$echo >= 2 && !is.null(con$echofile))
			cat("\n", file=con$echofile, append=TRUE)




	}
	if(!is.null(con$echofile))
		setwd(dirbefore)
		
	if(method == "Broyden")	
		exitcode <- (iter >= con$maxit)*1 + (sqrt( sum( fk^2 ) ) > con$tol)*10 + any(is.nan(Bk))*100	
	if(method != "Broyden")	
		exitcode <- (iter >= con$maxit)*1 + (sqrt( sum( fk^2 ) ) > con$tol)*10 + any(is.nan(Jacfk))*100	
	
	
	list(par= as.vector(xk), value=sqrt( sum( fk^2 ) ), counts=iter, iter=iter, 
		 code= exitcode, message=NULL, dirname=dirname, echofile=con$echofile)

}

graphlinesearch <- function(f, xk, dk, k, tbound=c(0, 2))
{
	phi <- function(t) sapply(t, function(u) crossprod(f(xk + u* dk)) /2)
	times <- seq(tbound[1], tbound[2], length=201)
	yval <- phi(times)
	

	
	plot(times, yval, type="l", main=paste("phi_k(t) (k=", k, ")", sep=""), ylim=c(min(yval), yval[1]), xlab="t", ylab="phi_k(t)")
	abline(h=phi(0), lty=2, col="grey")
	legend("topright", legend=paste("xk", as.vector(xk)))
	
	
	
}	
	
graphContourOneIter <- function(f, xk, dk, k, delta=.2, tols=c(1e-1, 1e-3, 1e-6), 
	cex=2, line=TRUE, tbound=c(0, 1/2, 1, 2))
{	
	if(length(xk) == 2)
	{
		xkp1 <- xk + dk
		xlim <- cbind( min = (1-delta)*pmin(xk, xkp1)*(pmin(xk, xkp1) >0) + (1+delta)*pmin(xk, xkp1)*(pmin(xk, xkp1) <0) - delta*(pmin(xk, xkp1) == 0),
				max = (1+delta)*pmax(xk, xkp1)*(pmax(xk, xkp1) >0) + (1-delta)*pmax(xk, xkp1)*(pmax(xk, xkp1) <0) + delta*(pmin(xk, xkp1) == 0) 	)	
				
		x <- seq(xlim[1, "min"], xlim[1, "max"], length=101)		
		y <- seq(xlim[2, "min"], xlim[2, "max"], length=101)				
	
		normf <- function(x, y)
		{
			xy <- cbind(x, y)
			sapply(1:NROW(xy), function(i) crossprod(f(xy[i, ])/2) )
		}
		z <- outer(x, y, normf )
		colnames(z) <- paste("y", y)
		rownames(z) <- paste("x", x)
		
		
		
		contour(x, y, z, nlevels=20, main="Contour of norm(F(x)), with 'zeros'", cex=cex, col="grey50")

		idzero <- which(abs(z) < tols[1], arr.ind=TRUE)
		points(cbind( x[idzero[, "row"]], y[idzero[, "col"]] ), col="yellow", pch=".", cex=1.5)

		idzero <- which(abs(z) < tols[2], arr.ind=TRUE)
		points(cbind( x[idzero[, "row"]], y[idzero[, "col"]] ), col="orange", pch=20)

		idzero <- which(abs(z) < tols[3], arr.ind=TRUE)
		points(cbind( x[idzero[, "row"]], y[idzero[, "col"]] ), col="red", pch=19)

		if(line)
		{
			trials <- function(t) sapply(t, function(u) xk + u* dk)
			

			times <- c(tbound[1], tbound[2], tbound[3], tbound[4])
			vals <- t(trials(times))
	
			lines(vals[3:4,1], vals[3:4,2], col="black", lty="dotted")
			lines(vals[2:3,1], vals[2:3,2], col="black", lty="dotdash")
			lines(vals[1:2,1], vals[1:2,2], col="black", lty="solid")
			points(vals[1,1], vals[1,2], pch=4)	
			
			legend("topright", legend=c("t<=1/2", "t<=1", "t<=2"), lty=c("dotted","dotdash","solid"), col="black")
		}
		
	}
	
}


graphContourMultIter <- function(f, xks, xlim, 
	option, nbgrid=101, tols=c(1e-1, 1e-3, 1e-6), 
	cex=1, lwd=1, itercol=c("blue", "green"), zerocol=c("yellow", "orange","red"), 
	title=NULL, posleg=NULL)
{	
	if(missing(xlim))
		xlim <- apply(xks, 2, range)	
	if(is.null(colnames(xlim)))	
		colnames(xlim) <- c("min", "max")
	if(any( colnames(xlim) != c("min", "max") ) )
		colnames(xlim) <- c("min", "max")
	
	option <- match.arg(option, c("none", "zeros", "colored", "numbered"), several.ok=TRUE)

	if(!is.null(title))
	{
		mytitle <- paste("Contour of ||F(x)||", title)
		
	}else if(length(option) == 1)
	{
		mytitle <- switch(option, 
		none= "Contour of norm(F(x))",
		zeros = "Contour of norm(F(x)), with 'zeros'",	
		colored = "Contour of norm(F(x)), with colored iterates",
		numbered = "Contour of norm(F(x)), with numbered iterates")		
	}else
		mytitle <- "Contour of norm(F(x)), with iterates and zeros"	
	
	
						
		x <- seq(xlim[1, "min"], xlim[1, "max"], length= nbgrid)		
		y <- seq(xlim[2, "min"], xlim[2, "max"], length= nbgrid)				
	
		normf <- function(x, y)
		{
			xy <- cbind(x, y)
			sapply(1:NROW(xy), function(i) crossprod(f(xy[i, ])/2) )
		}
		z <- outer(x, y, normf )
		colnames(z) <- paste("y", y)
		rownames(z) <- paste("x", x)
		

		
		contour(x, y, z, nlevels=20, main=mytitle, cex=cex, col="grey50", xlab="x_1", ylab="x_2")

	if("zeros" %in% option)		
	{
		idzero <- which(abs(z) < tols[1], arr.ind=TRUE)
		points(cbind( x[idzero[, "row"]], y[idzero[, "col"]] ), col=zerocol[1], pch=".", cex=cex)

		idzero <- which(abs(z) < tols[2], arr.ind=TRUE)
		points(cbind( x[idzero[, "row"]], y[idzero[, "col"]] ), col=zerocol[2], pch=20)

		idzero <- which(abs(z) < tols[3], arr.ind=TRUE)
		points(cbind( x[idzero[, "row"]], y[idzero[, "col"]] ), col=zerocol[3], pch=19)
		
		if(!is.null(posleg))
			legend(posleg, legend=c(paste("<", rev(tols), sep=""), paste(">", tols[1], sep="")), 
				   fill=c(rev(zerocol), "white"), bg ="white")

	}

	if("colored" %in% option)
	{
		mycol <- colorRampPalette( itercol )(NROW(xks)) 

		points(xks[,1], xks[,2], col=mycol, pch=3, cex=cex, lwd=lwd)		
	}
	if("numbered" %in% option)
	{
		n <- NROW(xks)
		mycol <- colorRampPalette( itercol )(n) 

		text(xks[,1], xks[,2], 0:(n-1), col=mycol, cex=cex, lwd=lwd)		
	}
}


