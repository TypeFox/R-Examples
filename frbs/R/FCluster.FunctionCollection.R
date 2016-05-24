#############################################################################
#
#  This file is a part of the R package "frbs".
#
#  Author: Lala Septem Riza
#  Co-author: Christoph Bergmeir
#  Supervisors: Francisco Herrera Triguero and Jose Manuel Benitez
#  Copyright (c) DiCITS Lab, Sci2s group, DECSAI, University of Granada.
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' This function is a part of the DENFIS method to generate cluster centers. 
#'
#' @title Evolving Clustering Method
#'
#' @param data.train a matrix (\eqn{m \times n}) of data for training, where \eqn{m} is the number of instances and 
#' \eqn{n} is the number of variables where the last column is the output variable.
#' @param Dthr the threshold value for the evolving clustering method (ECM), between 0 and 1.
#' @seealso \code{\link{DENFIS}} and \code{\link{DENFIS.eng}}
#' @return a matrix of cluster centers
#' @export
ECM <- function(data.train, Dthr){

Cc1 <- data.train[1, ,drop = FALSE]

Ru1 <- 0
nrow.dt <- nrow(data.train)
ncol.dt <- ncol(data.train)
Cc.j <- Cc1
Ru.j <- matrix(Ru1)

D.ij <- matrix()
temp <- matrix()
temp <- c(1)
n.cc <- 1
D.ij <- matrix()

num.cls <- 1
num.mem <- 1
for (i in 2 : nrow.dt){
	for (j in 1 : nrow(Cc.j)){
		norm.euc <- dist(rbind(data.train[i,], Cc.j[j, ]))
		D.ij[j] <- (norm.euc)/sqrt(ncol.dt)
	}
	
	D.ij <- matrix(D.ij)
	indx <- which.min(D.ij)	
	
	if (any(D.ij[indx] <= Ru.j)){
		#do nothing
	}
	else {
		S.ij <- Ru.j + D.ij
		indx.s <- which.min(S.ij)
		conts <- (2 * Dthr)
		if (S.ij[indx.s] > conts){
			Cc.j <- rbind(Cc.j, data.train[i, ])
			Ru.j <- rbind(Ru.j, 0)
		} else {		
			Ru.j.temp <- 0.5 * S.ij[indx.s]	
			
			if (Ru.j.temp > Dthr){
				Cc.j <- rbind(Cc.j, data.train[i, ])
				Ru.j <- rbind(Ru.j, 0)
			} else {
				temp <- data.train[i, ] - Cc.j[indx.s, ]
				d.temp <- sqrt(sum(temp^2))
				ratio <- abs((d.temp - Ru.j.temp))/d.temp
				new.vec <- ratio * temp
				new.cls <- Cc.j[indx.s, ] + new.vec
				Cc.j[indx.s, ] <- new.cls
				Ru.j[indx.s] <- Ru.j.temp
			}
		}
	}
}

res <- Cc.j

return(res)
}

# This function is to calculate the potential of data point that will be used on 
# subtractive clustering (SBC) method. 
#
# @title The potential of data
# @param dt.norm a matrix (n x m) of the normalized data.
# @param alpha a parameter of distance  
# @seealso \code{\link{SBC}} and \code{\link{SBC.test}}
# @return the potential on each data
# @export
potential <- function(dt.norm, alpha){
	counter <- 0
	n <- nrow(dt.norm)
	potential <- matrix(nrow = n)
	
	for (i in 1:n){
		temp.x <- dt.norm[i,]
		for (j in 1:n){
			euc <- dist(rbind(temp.x, dt.norm[j,]))
			cum <- exp(-alpha * (euc)^2)
			counter <- cum + counter
		}
		potential[i] <- counter
		counter <- 0
	}	
	return(potential)
}

# This function is to revise of the potential. 
#
# @title The revision of the potential
# @param P.star a value of potential of cluster center
# @param x.star a vector of cluster center
# @param dt.norm a matrix(n x m) of the normalized data.
# @param Beta a parameter of distance  
# @seealso \code{\link{SBC}} and \code{\link{SBC.test}}
# @return the revision of potential on each data
# @export
revise.pot <- function(P.star, x.star, dt.norm, Beta){
	counter <- 0
	n <- nrow(dt.norm)
	potential <- matrix(nrow = n)
	
	temp.x <- x.star
	for (i in 1:n){
		euc <- dist(rbind(temp.x, dt.norm[i,]))
		potential[i] <- exp(-Beta * (euc)^2)		
	}
	
	rev.pot <- P.star * potential 
	return(rev.pot)
}

# This function is to get stopping criteria on SBC
#
# @title The stopping criteria on Subtractive Clustering (SBC)
# @param r.a a positive constant which is effectively the radius defining a neighborhood. (default = 0.5)
# @param eps.high a parameter which is upper threshold value. (default = 0.5)
# @param eps.low a parameter which is lower threshold value. (default = 0.15)
# @param PP.star a value of potential of cluster center on k-th iteration
# @param P.star a value of potential of cluster center on 1-st iteration
# @param cluster.ctr a matrix of cluster center
# seealso \code{\link{SBC.test}} and \code{\link{SBC}}
# @return stopping criteria
# @export
stop.criteria <- function(r.a=0.5, eps.high = 0.5, eps.low = 0.15, PP.star, P.star,  cluster.ctr){
	temp.c.ctr <- na.omit(cluster.ctr)
	
	if (PP.star > (eps.high * P.star)){
		st.cr <- c(1)
	}
	else if (PP.star < (eps.low * P.star)){
		st.cr <- c(2)
	} else {
		x.k <- temp.c.ctr[nrow(temp.c.ctr), ]
		for (i in 1: (nrow(temp.c.ctr) - 1)){
			dmin <- dist(rbind(x.k, temp.c.ctr[i, ]))
			crt <- ((dmin/r.a) + (PP.star/P.star))
			if (crt >= 1){
				st.cr <- c(1)
			} else {
				st.cr <- c(3)
			}
		}
	}
	return(st.cr)
}

## it is used to calculate degree of membership function
# @param data.tst a matrix of testing data
# @param cluster.cls a matrix of cluster centers
# @param Dthr the threshold value for the evolving clustering method (ECM), between 0 and 1. 
# @param d a parameter for the width of the triangular membership function.
calc.degree.MF <- function(data.tst, cluster.cls, d, Dthr){
num.dt <- nrow(data.tst)
num.cls <- nrow(cluster.cls)
cluster.c <- cluster.cls[, -ncol(cluster.cls), drop = FALSE]

miu.rule <- matrix(nrow = num.dt, ncol = num.cls)
for (i in 1 : num.dt){
	for (j in 1 : num.cls){
		a <- cluster.c[j, ] - d * Dthr
		cc <- cluster.c[j, ] + d * Dthr
		left <- (data.tst[i, ] - a)/
		        (cluster.c[j, ] - a)
		right <- (cc - data.tst[i, ])/
		         (cc - cluster.c[j, ])		
		temp <- prod(pmax(pmin(left, right), 0))
		miu.rule[i, j] <- temp
	}	
}

return(miu.rule)
}