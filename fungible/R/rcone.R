##############################################
# Sept 26, 2009                              #
# rcone                                      #
# generate regression vectors with constant  #
# cosine with a target vector                #
##############################################


rcone <- function(R,Rsq,b,axis1,axis2,deg,Npoints=360) {
	
	# R			- Correlation matrix
	# Rsq		   - Regression Coefficient of Determination: corr(y, yhat)^2
	# b			- target vector 
	# axis1		- 1st axis of the rotation plane
	# axis2		- 2nd axis of the rotation plane
	# deg      - all vectors, b.i, will be "deg" degress from b
	# Npoints	- number of rotation vectors; default is 360
	
# convert degrees to radians
  d2r <- function(deg) pi/180 * deg
	
# convert radians to degrees
  r2d <- function(r) r*180/pi

# scale vector to unit length
	norm <- function(x) as.matrix(x/as.numeric(sqrt(t(x) %*% x)))
	
# cosine btwn two vectors	
   vec.cos <- function(x,y) t(norm(x)) %*% norm(y)
   	
# create rotation matrix
    Trot <- function(Rxx,deg){
	        	rot <- diag(ncol(Rxx))
				rot[2,2] <-  cos(d2r(deg))
				rot[2,3] <- -sin(d2r(deg))
				rot[3,2] <-  sin(d2r(deg))
				rot[3,3] <-  cos(d2r(deg))
				rot
			}
			
# Number of predictors
	Nvar <- length(b)

# norm b to unit length
   unit.b <- norm(b)

# B^-1 will rotate [b, axis1, axis2] to [e1,e2,e3]   
 if(Nvar ==3) B <- cbind(unit.b,axis1,axis2)

# fill out orthogonal basis for b 
 if(Nvar > 3){  
      GenB <- function(b, axis1,axis2){
        n <- length(b)
        B <-matrix(0,n,n)
        B[,1:3]<-cbind(b,axis1,axis2)

        for(i in 4:n){
           B[,i] <- resid(lm(rnorm(n)~-1+B[,1:(i-1)]))
        }

# norm all columns in the basis        
        d <- diag(1/sqrt(diag(crossprod(B))))
        B <- B%*%d
        B
     }#end GenB

    B <- GenB(unit.b,axis1,axis2)
	 }
	 
# create an arbitrary vector w that makes an angle
# of "deg" degrees with e1
	w <- c(cos(d2r(deg)),sin(d2r(deg)),rep(0,(Nvar-2)))

# Rotate w around e1 in the e2,e3 plane
	w.i <- matrix(0,Nvar,Npoints)
	for(i in 1:Npoints){
		 w.i[,i] <- Trot(R,(360/Npoints)*i) %*% w
	}	 
	
# Rotate w.i to original basis
# b.u = unscaled b.i 
  b.u <- B %*% w.i
	
# Scale b.u to terminate on the ellispoid
	b.i <- apply(b.u,2,function(x) x*sqrt(Rsq/as.numeric(t(x)%*%R%*%x)))    

	return(b.i)
}  #END FUNCTION rcone

#=======================================================================#

# source("WAIS3.DAT")
# 
# Npoints <- 1000
# Rsq <- .40
# NumDeg <- 20
# V <- eigen(R)$vectors
# bvec <- c(1,2,6,8)   #8
# ABCD <- c("A","B","C","D")
# 
# #postscript(file="~/Desktop/FigCone.eps", width=8, height=11)
# par(mfrow=c(2,2))
# for(plot.num in 1:4){
# 	# create b parallel to v[plot.num]
#    # rotate in the 7 - 14 plane
#    b <- V[,bvec[plot.num]]
#    bsq <- t(b) %*% R %*% b 
#    b <- b * sqrt(Rsq/bsq)                
#    b.i <- rcone(R,Rsq,b, V[,7],V[,14], deg=NumDeg,Npoints)
# 
#    #compute validity vectors
#    r <- R %*% b.i
#   	Rsq.r <- Rsq.unit <- rep(0,Npoints)
# 	for(i in 1:Npoints){
# 		# performance of unit weights
# 		Rsq.unit[i] <- (t(sign(r[,i])) %*% r[,i])^2 /
# 		               (t(sign(r[,i])) %*% R %*% sign(r[,i]))
# 		# performance of correlation weights               
#      	Rsq.r[i] <- (t(r[,i]) %*% r[,i])^2 /(t(r[,i]) %*% R %*% r[,i])	}
# 	
#    plot(seq(0,360,length=Npoints),Rsq.r, typ="l", 
#          ylim=c(0,Rsq),
#          xlim=c(0,360),
#          lwd=3,
#          ylab=expression(R^2),
#          xlab="Degrees around the cone ",
#          cex.lab=1.25,lab=c(10,5,5),
#          main=paste("b  ||  v",bvec[plot.num],sep=""))
#    points(seq(0,360,length=Npoints), Rsq.unit, 
#          type="l", 
#          lty=2,lwd=3)
#    mtext(text=ABCD[plot.num],at=5,cex=1.5,line=1.5)      
#    if(plot.num==1){
#    	legend(x=10,y=.30,
#    	      legend=c("unit weights", "r-weights"),
#    	      lty=c(1,2),lwd=3,cex=1.3)
#    }	      
#  } # End create plots
# 
# #dev.off()       

