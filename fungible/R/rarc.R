################################################################## 
# Sept 21 2009                                                   #
# rarc:  rotate between two poionts on the surface on            #
# an n-dimensional ellipsoid                                     #
#                                                                #
# Rxx     : n-dimensional Correlation matrix                     #
# b1, b2  : scalars denoting regression vectors parallel to      #
#           eigenvectors v1 and v2  OR                           #
#           regression vectors on the n-dimensional ellipsoid    #
#           defined by b'Rxx b = Rsq                             #
# Rsq     : scalar defining the relative size of the             #
#           n-dimensional ellipsoid                              #
# Npoints : number of points in elliptical arc                   #               
##################################################################


rarc <- function(Rxx,Rsq,b1,b2,Npoints) {
		

# if b1 and b2 are scalars then choose scaled eigenvectors 
# v[b1] and v[b2] as the start and end vectors 
	
 if(length(b1)==1){
 	  cat(paste("\nb1 and b2 parallel to eigenvectors v",
 	  b1, " and v",b2,sep=""))
	  V <- eigen(Rxx)$vectors
	  b1 <- V[,b1]; if(sum(b1)< 0) b1<- -1*b1
	  b2 <- V[,b2]; if(sum(b2)< 0) b2<- -1*b2
 	}  


# if b1 and b2 are vectors but b'R b != Rsq: 
# scale to appropriate length   
    if( !identical( t(b1) %*% Rxx %*% b1, Rsq)) {
    	b1 <- b1 * sqrt(Rsq/as.numeric(t(b1)%*%Rxx%*%b1)) 
    }	     
     if( !identical( t(b2) %*% Rxx %*% b2, Rsq)) {
     	b2 <- b2 * sqrt(Rsq/as.numeric(t(b2)%*%Rxx%*%b2)) 
    } 	  
  	   
   
#~~~~~~FUNCTION DEFINITIONS~~~~~~~~~~~~~~~~#

# convert degrees to radians
  d2r <- function(deg) pi/180 * deg
	
# convert radians to degrees
  r2d <- function(r) r*180/pi
	
# rotation matrix
  Trot <- function(Rxx,deg,rot1,rot2){
		rot <- diag(ncol(Rxx))
		rot[rot1,rot1] <-  cos(d2r(deg))
		rot[rot1,rot2] <-  sin(d2r(deg))
		rot[rot2,rot1] <- -sin(d2r(deg))
		rot[rot2,rot2] <-  cos(d2r(deg))
		rot
	}
	
# scale vector to unit length
	norm <- function(x) as.matrix(x/as.numeric(sqrt(t(x) %*% x)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	
# Begin Main Function	
		
	numX <- ncol(Rxx)		
	
# norm the regression vectors
	xhat <- norm(b1)
	yhat <- norm(b2)	
 
# Create rotation matrix M1 using qr() and qr.Q() functions
	tmp <- matrix(c(xhat,rnorm(numX*(numX-1))),numX,numX)
	M1 <- t(qr.Q(qr(tmp)))
	if(sign(M1[1,1])!=sign(xhat[1]))  M1[1,] <- -1*M1[1,]
 
# rotate xhat to e1
  z <- M1 %*% xhat
	
# create v = M1%*%yhat and norm
	v <- norm(M1 %*% yhat)	
  		 
# build rotation matrix to create w
	old.rot <- diag(numX)
	old.v <- v
	for(i in 1:(numX-2)) {
		theta <- atan(old.v[i+2]/old.v[2])
		M2 <- diag(numX)
		M2[2,2] <- M2[i+2,i+2] <- cos(theta)
		M2[2,i+2] <- sin(theta); M2[i+2,2] <- -sin(theta)
		old.v <- M2 %*% old.v
		M2 <- M2 %*% old.rot
		old.rot <- M2	
	}
	
# create w  
	w <- M2 %*% v
	
# Rotate z to w along shortest path
  sgnDegrees <- 1  
  if( (sign(w[1]) ==  1) & (sign(w[2])==1)) sgnDegrees <- -1
  if( (sign(w[1]) == -1) & (sign(w[2])==1)) sgnDegrees <- -1
		
# rotate z towards w
	degrees <- r2d(acos( t(w) %*% z) ) 
	cat("\nb1 and b2 are separated by ",degrees, " degrees\n")
	zstar <- matrix(0,numX,Npoints)
	for(i in 1:Npoints) {
		zstar[,i] <- Trot(Rxx, i*( sgnDegrees* degrees/Npoints),1,2)%*%z
	}	
	zstar <- cbind(z,zstar)
	
# create vstar
	vstar <- solve(M1) %*% solve(M2) %*% zstar
	
# scale vstar to terminate on ellipsoid
	c1 <- apply(vstar,2,function(x) t(x) %*% Rxx %*% x)
	b <- matrix(0,numX,(Npoints+1))
	for(i in 1:(Npoints+1)) b[,i] <- sqrt(Rsq/c1[i]) * vstar[,i]

   return(b)
	
}  # End Function ellipsoid.arc	

#=======================================================================#

# ##-------------------##
# ## EXAMPLE
# ## GRE/GPA Data
# ##-------------------##
# R <- Rxx <- matrix(c(1.00, .56, .77,
#                      .56, 1.00, .73,
#                      .77, .73, 1.00),3,3)
#                  
# # GPA validity correlations                 
# rxy <- c(.39, .34, .38)
# b <- solve(Rxx)%*%rxy
# 
# Rsq <- t(b) %*% Rxx %*% b
# N <- 200       
#                   
# b <- rarc(Rxx=R,Rsq,b1=1,b2=3,Npoints=N) 
# 
# #compute validity vectors
#   r <- Rxx %*% b
# 	N <- N+1
# 	Rsq.r <- Rsq.unit <- rep(0,N)
# 	for(i in 1:N){
# 		# performance of unit weights
# 		Rsq.unit[i] <- (t(sign(r[,i])) %*% r[,i])^2 /
# 		               (t(sign(r[,i])) %*% R %*% sign(r[,i]))
# 		# performance of correlation weights               
#      	Rsq.r[i] <- (t(r[,i]) %*% r[,i])^2 /(t(r[,i]) %*% R %*% r[,i])	}
# 
# cat("\nAverage relative performance of unit weights across elliptical arc:",
# 	    round(mean(Rsq.unit)/Rsq,3) )     
# cat("\n\nAverage relative performance of r weights across elliptical arc:",
# 	    round(mean(Rsq.r)/Rsq,3) ) 
# 
# postscript(file="~/Desktop/FigArcEx.eps", width=8, height=11)
# plot(seq(0,90,length=N),Rsq.r, typ="l", 
#          ylim=c(0,.20),
#          xlim=c(0,95),
#          lwd=3,
#          ylab=expression(R^2),
#          xlab=expression(paste("Degrees from ",b[1]," in the direction of ",b[2])),
#          cex.lab=1.25,lab=c(10,5,5))
# points(seq(0,90,length=N), Rsq.unit, 
#          type="l", 
#          lty=2,lwd=3)
# legend(x=0,y=.12,
#        legend=c("r weights", "unit weights"), 
#        lty=c(1,2),
#        lwd=c(4,3),
#        cex=1.5) 
# dev.off()               
# 
# 



