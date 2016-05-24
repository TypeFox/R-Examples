## Functions for performing Q-ball reconstruction from the  dti-package
## Authors: Karsten Tabelow, Joerg Polzehl 
## http://www.jstatsoft.org/v44/i12/
## Copyright (C) 2005-2012 Weierstrass
## Institute for Applied Analysis and Stochastics.
## URL: http://www.wias-berlin.de/projects/matheon_a3

sphcoord <-
function(ccoord)
{
#
#  transform cartesian into sherical coordinates
#
  ccoord <- ccoord/sqrt(sum(ccoord^2))
  phi <- atan2(ccoord[2],ccoord[1])
  theta <- atan2(sqrt(ccoord[2]^2+ccoord[1]^2),ccoord[3])
  c(theta,phi)
}

design.spheven <-
function(order,gradients,lambda=NULL)
{
#
#  compute design matrix for Q-ball
#
  if(is.null(lambda)) lambda <- 0.006
  order <- as.integer(max(0,order))
  if(order%%2==1){
    warning("maximum order needs to be even, increase order by one")
    order <- order+1
  } 
  # calculate spherical angles theta and phi corresponding to the gradients
  n <- dim(gradients)[2]
  theta <- phi <- numeric(n)
  for( i in 1:n){
    angles <- sphcoord(gradients[,i])
    theta[i] <- angles[1]
    phi[i] <-  angles[2]
  }
  # values of SH on specified spherical angles
  sphharmonics <- getsphericalharmonicseven(order,theta,phi)
  # Laplace-Beltrami-Regularization term
  lord <- rep(seq(0,order,2),2*seq(0,order,2)+1)
  L <- lambda*diag(lord^2*(lord+1)^2)
  # transformation matrix for SH coefficients
  ttt <- solve(sphharmonics%*%t(sphharmonics)+L)%*%sphharmonics
  # results
  list(design = sphharmonics,
       matrix = ttt,
       theta = theta,
       phi = phi)
}

plzero <-
function(order)
{
  l <- seq(2,order,2)
  pl <- l
  for(i in 1:length(l)) pl[i] <- (-1)^(l[i]/2)*prod(seq(1,(l[i]-1),2))/prod(seq(2,l[i],2))
  2*pi*diag(rep(c(1,pl),2*seq(0,order,2)+1))
}

getsphericalharmonicseven <-
function(order,theta,phi)
{
#
#   compute spherical harmonics
#
	# require(gsl,quietly = TRUE,warn.conflicts = FALSE)
	order <- as.integer(max(0,order))
	if(order%%2==1){
		warning("maximum order needs to be even, increase order by one")
		order <- order+1
	} 
	if(length(theta)!=length(phi)) stop("need same length of theta and phi")
	kseq <- seq(0,order,2)
	n <- length(phi)
	values <- matrix(0,(order+1)*(order+2)/2,n)
	for(k in kseq){
		mseq <- seq(-k,k,1)
		for(m in mseq){
			ind <- (k^2+k+2)/2+m
			z <- gsl::legendre_sphPlm(k,abs(m),cos(theta))
			if(m < 0){
				z <- sqrt(2)*z*cos(m*phi)
			} 
			if(m > 0){
				z <- sqrt(2)*z*sin(m*phi)
			}
			values[ind,] <- z
		}
	}
	# detach(package:gsl)
	values
}

getsphericalharmonicsall <-
function(order,theta,phi)
{
#
#   compute spherical harmonics
#
  order <- as.integer(max(0,order))
  if(length(theta)!=length(phi)) stop("need same length of theta and phi")
  kseq <- 0:order
  n <- length(phi)
  values <- matrix(0,(order+1)^2,n)
  l <- 1
  for(k in kseq){
    mseq <- (-k):k
    for(m in mseq){
      z <- gsl::legendre_sphPlm(k,abs(m),cos(theta))
      if(m < 0){
        z <- sqrt(2)*z*cos(m*phi)
      } 
      if(m > 0){
        z <- sqrt(2)*z*sin(m*phi)
      }
      values[l,] <- z
      l <- l+1
    }
  }
  #detach(package:gsl)
  values
}


#-----------------------------
##
## Aganj data transformation 
## Regularization following Aganj et al. (2010) delta=1e-3
## 
datatrans <-
function(si, s0=1)
{
	# cat("Data transformation started ",date(),"\n")
	dim(s0) <- dim(si) <- NULL
	si <- si/s0
	si[is.na(si)] <- 0
	si[si>=1] <- 1-.Machine$double.neg.eps
	si <- log( -log(si))
	si[is.na(si)] <- 0
	si[(si == Inf)] <- 0
	si[(si == -Inf)] <- 0
	# dim(si) <- c(prod(ddim),ngrad0)
##	si <- t(si)
	# cat("Data transformation completed ",date(),"\n")
	return(si)
}
#-----------------------------

