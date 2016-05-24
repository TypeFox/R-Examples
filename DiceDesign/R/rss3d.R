rss3d<-function(design, lower, upper, gof.test.type="greenwood", gof.test.stat=NULL, transform=NULL, n.angle=60, graphics=1, trace=TRUE){


design <- as.matrix(design)
n <- dim(design)[1]
d <- dim(design)[2]

#if (ncol(design)!=3) stop("'design' must have 3 columns")
#if ((length(lower)!=3) & (length(upper)!=3)) stop("lower and upper must be 3-dimensional vectors")
#dimension <- 3

	# some arguments checks

if (d < 3) stop("rss3d is for d-dimensional designs with d >= 3")	
if (!is.element(gof.test.type, c("greenwood", "qm", "ks", "V", "cvm"))) stop("The goodness-of-fit test must be one of: Greenwood, Quesenberry-Miller, Kolmogorov-Smirnov, V = (D+) +  (D-), or Cramer-Von Mises")

if ((length(lower)!=d) & (length(upper)!=d)) stop("lower and upper must be d-dimensional vectors")


	# domain transformation to [-1,1]^d 

for (j in 1:d) {
	design.col <- design[,j]
	design.col.min <- min(design.col)
	if (design.col.min < lower[j]) stop('the minimum of design values is not compatible with lower')
	design.col.max <- max(design.col)
	if (design.col.max > upper[j]) stop('the maximum of design values is not compatible with upper')
	design.col <- 2*((design.col - lower[j])/(upper[j] - lower[j]) - 0.5)
	design[,j] <- design.col
}


	# angles definition

theta.degree <- seq(from=0, to=180, length=n.angle+1)
theta.degree <- theta.degree[1:n.angle]
theta.degree <- as.matrix(theta.degree)
theta <- theta.degree*2*pi/360
n.theta <- length(theta)        
cos.theta <- cos(theta); sin.theta <- sin(theta)

phi.degree <- seq(from=-90, to=90, length=n.angle+10)
phi.degree <- phi.degree[1:n.angle]
phi.degree <- as.matrix(phi.degree)
phi <- phi.degree*2*pi/360
n.phi <- length(phi)
cos.phi <- cos(phi); sin.phi <- sin(phi)


	# Loops on theta and phi to build the matrix of statistic values, one of which is implemented in C

ax <- cos.theta%*%t(cos.phi)
ay <- sin.theta%*%t(cos.phi)
az <- rep(1, n.theta)%*%t(sin.phi)

print.out <- data.frame(global.stat = rep(NA, (d*(d-1)*(d-2))%/%6))

	# loop over all triplets of dimensions

anglewise.stat.max <- anglewise.stat <- matrix(0, n.theta, n.phi)
global.stat.array <- array(NA, c(d, d, d))
global.stat.max <- 0
meter <- 0

if (trace) {
	cat("\n3D Radial Scanning Statistic (RSS) with ", toupper(gof.test.type), " statistic\n", sep="")
	cat("Discretization step (in degree) : ", 180/n.angle, sep="")
	cat("\n\nMaximum of RS statistic values (global statistic) per triplet of dimensions")
}

for (i1 in 1:(d-2)) {
	x <- design[, i1]

	for (i2 in ((i1+1):(d-1))) {
		y <- design[, i2]

		for (i3 in ((i2+1):d)) {
			z <- design[, i3]


			for (j in 1:n.phi){                  

					# 1st step : compute the matrix of the F(projected points onto Vect(cos.theta, sin.theta))
							
				out <- .C("C_rss3Dloop", as.double(x), as.double(y), as.double(z), as.double(ax), as.double(ay), as.double(az), as.integer(n), as.integer(n.theta), as.integer(n.phi), as.integer(j-1), ans=double(n * n.theta), PACKAGE="DiceDesign") 

				F.projections <- matrix(out$ans, n, n.theta)

	
					# 2nd step: test statistic computations

				for (k in 1:n.theta) {
					anglewise.stat[k,j] <- unif.test.statistic(x=F.projections[,k], type=gof.test.type, transform=transform)
				}	# end loop over theta (k)

			} # end loop over phi (j)
	
				# compute the worst value over all angles and store it
			global.stat <- max(anglewise.stat)      
			global.stat.array[i1,i2,i3] <- global.stat
		
			if (global.stat > global.stat.max) {
				global.stat.max <- global.stat
				triplet.worst <- c(i1,i2,i3)
				anglewise.stat.max <- anglewise.stat
			} 
			
			meter <- meter + 1
			print.out[meter, 1] <- global.stat
			name.current <- paste("(", i1, ",", i2, ",", i3, ")", sep="")
			row.names(print.out)[meter] <- name.current
			if (trace) cat("\n", name.current, " ", global.stat, sep="")
	
		} # end loop i3
	} # end loop i2
} # end loop i1

if (trace) cat("\n\n")
	
	# threshold at significance level 5%
	
if (is.null(gof.test.stat)) {
	gof.test.stat <- unif.test.quantile(type=gof.test.type, n=n, alpha=0.05)
}

	
	# 3D graph with package rgl ##

if (graphics > 0) {
	
  if (requireNamespace("rgl", quietly = TRUE)) {
    rgl::open3d()
    design <- design[, triplet.worst]
    x <- design[, 1]
    y <- design[, 2]
    z <- design[, 3]
    
    index.max <- which.max(anglewise.stat.max)
    ax.max <- ax[index.max]
    ay.max <- ay[index.max]
    az.max <- az[index.max] 
    
    phi.max <- theta.max <- NA
    for (j in 1:n.phi) {
      for (k in 1:n.theta) {
        if (abs(anglewise.stat.max[k,j]-global.stat.max)<1e-10) {
          phi.max <- phi[j]
          theta.max <- theta[k]
        } 			
      }
    }
    
    dir.max <- c(ax.max, ay.max, az.max)
    projections <- as.matrix(design) %*% dir.max
    
    rgl::plot3d(x, y, z, size=5, col="blue", xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1), xlab="", ylab="", zlab="")
    
    a.max <- max(abs(dir.max))
    
    rgl::plot3d(c(-1,1)*ax.max/a.max, c(-1,1)*ay.max/a.max, c(-1,1)*az.max/a.max, col="red", type="l", size=2, add=TRUE)
    
    
    for (i in 1:n) {
      h <- projections[i]*dir.max - design[i,]
      if (max(abs(h)) > 1e-8) {
        lambda <- min(min((sign(h) - design[i,])/h), 1)
        rgl::plot3d(c(design[i,1], design[i,1] + lambda*h[1]), c(design[i,2], design[i,2] + lambda*h[2]), c(design[i,3], design[i,3] + lambda*h[3]), type="l", col="red", add=TRUE)
      }
      if (max(abs(projections[i]*dir.max))<=1) rgl::plot3d(projections[i]*dir.max[1], projections[i]*dir.max[2], projections[i]*dir.max[3], pch=20, col="red", size=5, add=TRUE)
    }
    par(mfrow=c(1,1))
    
  } else {
   print("Error : the package rgl is not installed")
  }
	
	
} # end of conditional block: "if graphics > 0"

return(list(stat=anglewise.stat.max, angle=data.frame(theta=theta, phi=phi), global.stat=global.stat.array, print.out=print.out, gof.test.stat=gof.test.stat, worst.case=triplet.worst, worst.dir=dir.max))

} 