rss2d <- function(design, lower, upper, gof.test.type="greenwood", gof.test.stat=NULL, transform=NULL, n.angle=360, graphics=1, trace=TRUE, lines.lwd=1, lines.lty="dotted", ...)  {

design <- as.matrix(design)
n <- dim(design)[1]
d <- dim(design)[2]


	# some arguments checks
	
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


	# compute the subdivision of the (half) circle in cartesian coordinates

theta.degree <- seq(from=0, to=180, length=n.angle+1)
theta.degree <- theta.degree[1:n.angle]
theta.degree <- as.matrix(theta.degree)
theta <- theta.degree*2*pi/360
cos.theta <- cos(theta)
sin.theta <- sin(theta)
n.theta <- length(theta)
subdiv.halfcircle <- cbind(cos.theta, sin.theta)


	# loop over all pairs of dimensions

global.stat <- matrix(NA, d, d)
global.stat.max <- 0

if (trace) {
	cat("\n2D Radial Scanning Statistic (RSS) with ", toupper(gof.test.type), " statistic\n", sep="")
	cat("Discretization step (in degree) : ", 180/n.angle, sep="")
	cat("\n\nMaximum of RS statistic values (global statistic) per pair of dimensions")
}


print.out <- data.frame(global.stat = rep(NA, (d*(d-1))%/%2))
meter <- 0

for (i in 1:(d-1)) {
	x <- design[,i]

	for (j in ((i+1):d)) {
		y <- design[,j]
		
			# compute anglewise statistic
		
				
				# 1st step : compute the matrix of the F(projected points onto Vect(cos.theta, sin.theta))
		out <- .C("C_rss2Dloop", as.double(x), as.double(y), as.double(cos.theta), as.double(sin.theta), as.integer(n), as.integer(n.theta), ans=double(n * n.theta), PACKAGE="DiceDesign") 
		F.projections <- matrix(out$ans, n, n.theta)

				# 2nd step : for each angle, compute the statistic 
				# In future version, should be computed inside the C loop
		anglewise.stat.ij <- matrix(NA, n.theta, 1)
		for (angle in 1:n.theta) {
			anglewise.stat.ij[angle] <- unif.test.statistic(x=F.projections[,angle], type=gof.test.type, transform=transform)
		}
			

			# compute the worst value over all angles and store it
		global.stat.ij <- max(anglewise.stat.ij)      
		global.stat[i,j] <- global.stat[j,i] <- global.stat.ij
		
		if (global.stat.ij > global.stat.max) {
			global.stat.max <- global.stat.ij
			pair.worst <- c(i,j)
			anglewise.stat <- anglewise.stat.ij
		} 
		
		meter <- meter + 1
		print.out[meter, 1] <- global.stat.ij
		name.current <- paste("(", i, ",", j, ")", sep="")
		row.names(print.out)[meter] <- name.current
		if (trace) cat("\n", name.current, " ", global.stat.ij, sep="")
		
	} # end loop j
} # end loop i
	
if (trace) cat("\n\n")

	
	# rss curve
rss.curve.x <- anglewise.stat*subdiv.halfcircle[,1]
rss.curve.y <- anglewise.stat*subdiv.halfcircle[,2]


	# statistic upper tail percentage points
	# see D'Agostino and Stephens "Goodness-of-fit techniques", 1986

if (is.null(gof.test.stat)) {
	gof.test.stat <- unif.test.quantile(type=gof.test.type, n=n, alpha=0.05)
}


	# --------------------
	# graphics if required
	# --------------------
	
if (graphics>=0) {

	design.names <- names(as.data.frame(design))
	
	design <- design[ , pair.worst] 
	design.names <- design.names[pair.worst]
	
	anglewise.stat.max <- max(anglewise.stat)
	index.max <- which.max(anglewise.stat)
 	cos.theta.max <- subdiv.halfcircle[index.max, 1]
 	sin.theta.max <- subdiv.halfcircle[index.max, 2]
 	dir.max <- c(cos.theta.max, sin.theta.max)
	
	if (is.element(graphics, c(0,1))) {
		par(mfrow=c(1,2+graphics))
		plot(design, xlim=c(-1,1), ylim=c(-1,1), 
		             xlab=design.names[1], ylab=design.names[2])
	}

	
		# draw the rss curve
	if (graphics>0) {
		rx <- c(rss.curve.x, -rss.curve.x, rss.curve.x[1])
		ry <- c(rss.curve.y, -rss.curve.y, rss.curve.y[1])
		graph.size <- max(abs((anglewise.stat.max)*dir.max), gof.test.stat)*1.05
		plot(rx, ry, xlim=c(-graph.size,graph.size), ylim=c(-graph.size, graph.size), 
									 xlab="", ylab="", ...)

			# draw the circle with radius equal to the threshold at significance level 5%
		theta_aux <- seq(from=0, to=2*pi+0.1,by=0.1)
		lines(gof.test.stat*cos(theta_aux), gof.test.stat*sin(theta_aux))

			# draw the coordinate axis in dotted lines
		abline(h=0, v=0, lty=lines.lty, col="black", lwd=lines.lwd)
	}
 	
 	if (is.element(graphics, c(0,1))) {
 		plot(design, xlim=c(-1,1), ylim=c(-1,1), 
		             xlab=design.names[1], ylab=design.names[2])
 		projections <- design %*% dir.max
		points(projections*dir.max[1], projections*dir.max[2], pch=20, col="red")
		if (cos.theta.max==0) {
			lines(c(0,0), c(-1,1), col="red")
		} else lines(c(-1,1), c(-1,1)*sin.theta.max/cos.theta.max, col="red")
		for (i in 1:n) lines(c(design[i,1], projections[i]*cos.theta.max), c(design[i,2], projections[i]*sin.theta.max), lty=lines.lty, lwd=lines.lwd)
 	}
	
	par(mfrow=c(1,1))
}

return(list(global.stat=global.stat, worst.case=pair.worst, worst.dir=dir.max, stat=as.numeric(anglewise.stat), angle=as.numeric(theta), curve=cbind(c(rss.curve.x,-rss.curve.x), c(rss.curve.y,-rss.curve.y)),  gof.test.stat=gof.test.stat))
}