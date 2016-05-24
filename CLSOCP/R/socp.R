######################################################################
#Copyright Jason Rudy & Faramarz Valafar 2009-2010

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
######################################################################

socp <-
function(A, b, c, kvec, type = rep('q',length(kvec)), use_sparse=TRUE, gamma_fac=.95, delta0 = .75, sigma0 = .25, mu0 = 0.01, zero_tol = 1e-6, max_iter = 100, min_progress = zero_tol/100){
	#A is m x sum(kvec)
	#b is a sum(kvec) vector
	#kvec = c(k1, ..., kn)
	
	#Set return code
	code = -1
	#Codes are: 
	#0 - Full success, 
	#1 - Singularity or other error occured,
	#2 - Lack of progress,
	#3 - Maximum number of iterations reached
	

	#Change linear constraints to SOC
	converted <- convert_l(kvec,type)
	kvec <- converted$kvec
	type <- converted$type
	
	
	if(use_sparse){
		A <- as(A,"sparseMatrix")
		cat("Using sparse constraint and jacobian matrices.\n")

	}
	
	#Dimensions
	n <- length(kvec)
	m <- dim(A)[1]
	k <- dim(A)[2]
	
	#Create the SOC index list
	inds <- make_inds(kvec, type)
	
	
	#Step 0
	delta <- delta0
	sigma <- sigma0
	mu <- rep(mu0,n)
	zbar <- c(rep(0,k + m), mu)
	x <- jordan_identity(kvec)
	y <- rep(0, m)
	s <- c - t(A)%*%y
	H <- CL_H(x, s, mu, A, b, inds)
	nrm_H <- two_norm(H)
	oldnrm_H <- nrm_H + 10*min_progress
	gamma <- gamma_fac*min(1,1/nrm_H)
	iter <- 0
	lklast <- 0
	
	while(TRUE){
		
		#Step 1
		if (nrm_H < zero_tol){
			code <- 0
			cat("Solution achieved within tolerance for socp.\n")
			break
		}else if(oldnrm_H - nrm_H < min_progress){
			code <- 2
			cat("Minimum progress not achieved for socp.\n")
			break
		}else if(iter >= max_iter){
			code <- 3
			cat("Maximum number of iterations reached for socp.\n")
			break
		}else{
			rho <- CL_rho(nrm_H, gamma)
		}
		
		#Step 2
		#print(rcond(CL_grad_H(x, y, s, mu, A, c, inds), useInve = TRUE))
		#First try Newton method
		delta_z <- try(solve(CL_grad_H(x, y, s, mu, A, c, inds), rho*zbar - H, tol=.Machine$double.xmin))
		if (class(delta_z)=="try-error"){
			
			code <- 1
			cat("Singularity or other error occured for socp.  Solution could be inaccurate.\n")
			break

		}
		
		#Steps 3 and 4
#		lk <- -1
#		done <- FALSE
#		while(!done){
#			lk <- lk +1
#			newx <- x + (delta^lk)*delta_z[1:k]
#			newy <- y + (delta^lk)*delta_z[(k+1):(k+m)]
#			newmu <- mu + (delta^lk)*delta_z[(k+m+1):(k+m+n)]
#			news <- c - t(A)%*%newy
#			newH <- CL_H(newx, news, newmu, A, b, inds)
#			newnrm_H <- two_norm(newH)
#			done <- newnrm_H <= (1-sigma*(1-gamma*mu0)*(delta^lk))*nrm_H
#		}
#		print("old lk")
#		print(lk)
		
		
		
#		Optimized based on the idea that lk from this iteration should be close to lk form the previous iteration
		lk <- lklast
		signswitch <- FALSE
		newx <- x + (delta^lk)*delta_z[1:k]
		newy <- y + (delta^lk)*delta_z[(k+1):(k+m)]
		newmu <- mu + (delta^lk)*delta_z[(k+m+1):(k+m+n)]
		news <- c - t(A)%*%newy
		newH <- CL_H(newx, news, newmu, A, b, inds)
		newnrm_H <- two_norm(newH)
		oldsign <- newnrm_H <= (1-sigma*(1-gamma*mu0)*(delta^lk))*nrm_H
		
		done <- FALSE
		while(!done){
	
			if(oldsign){
				if (lk == 0){
					break
				}
				lk <- lk - 1
				
			}else{
				lk <- lk + 1
			}
			
			newx <- x + (delta^lk)*delta_z[1:k]
			newy <- y + (delta^lk)*delta_z[(k+1):(k+m)]
			newmu <- mu + (delta^lk)*delta_z[(k+m+1):(k+m+n)]
			news <- c - t(A)%*%newy
			newH <- CL_H(newx, news, newmu, A, b, inds)
			newnrm_H <- two_norm(newH)
			sign <- newnrm_H <= (1-sigma*(1-gamma*mu0)*(delta^lk))*nrm_H
			
			if(sign & !oldsign){
				done <- TRUE
			}
			oldsign <- sign
			
		}
		
		
		lklast <- max(lk - 1,0)
		#cat("lk =",as.character(lk),"\n")
		
	
		x <- newx
		y <- newy
		mu <- newmu
		s <- news
		H <- newH
		oldnrm_H <- nrm_H
		nrm_H <- newnrm_H
		iter <- iter + 1
		cat("Iteration ", iter, " complete.  ||H|| = ",nrm_H,".\n",sep="")
	
		
	}
	
	
	return(list(x=as.vector(x), y=as.vector(y), s=as.vector(s), obj=sum(c*x), code=code, mu=mu, iter=iter))
	
}

