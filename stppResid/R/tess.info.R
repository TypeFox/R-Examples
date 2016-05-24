tess.info <- function(X, cifunction, theta = NULL, areas, tl, n.def = 100)
{
	data <- data.frame(cbind(X$x, X$y, X$t))
	names(data) <- c("x", "y", "t")
	if(is.null(theta)) {
		fun2 <- function(z)
		{
			Y <- rbind(data, z)
			Y <- Y[which(Y$t <= z[3]),]
			Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
			tail(cifunction(Y),1)
		}
	}
	if(!is.null(theta)) {
		fun2 <- function(z)
		{
			Y <- rbind(data, z)
			Y <- Y[which(Y$t <= z[3]),]
			Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
			tail(cifunction(Y, theta),1)
		}
	}
	all <- function(x)
	{
		div <- length(x)/3
		xt <- x[1:div]
		yt <- x[(div+1):(2*div)]
		tt <- x[(2*div+1):(3*div)]
		new.pts <- data.frame(cbind(xt, yt, tt))
		new.pts <- new.pts[order(new.pts[,3]),]
		lamb2 <- apply(new.pts, 1, fun2)
		return(lamb2)
	}
	new.points <- function(tl) {
		poly.temp <- as.points(tl$x, tl$y)
		pts <- csr(poly.temp, n.def)
		pts2 <- c(pts[,1], pts[,2])
		return(pts2)
	}
	vol <- areas * diff(X$tcoord)
	n <- 0
	lamb2 <- c()
	temp.tl <- tl
	place <- 1:length(temp.tl)
	lamb.ave.final <- rep(0, length(place))
	lamb.sd.final <- rep(0, length(place))
	n.final <- rep(0, length(place))
	n = n + n.def
	new.pts <- lapply(temp.tl, new.points)
	xyp <- data.frame(new.pts)
	xyp <- data.matrix(xyp)
	tp.t <- runif(n.def*length(temp.tl), X$tcoord[1], X$tcoord[2])
	tp <- matrix(tp.t, ncol = length(temp.tl))
	new.p <- rbind(xyp, tp)
	lamb2 <- rbind(lamb2, apply(new.p, 2, all))
	lamb2.ave <- apply(lamb2, 2, mean)
	lamb2.sd <- apply(lamb2, 2, sd)
	while(length(place) > 0) {
		n = n + n.def
		cat("number of bins left: ", length(temp.tl), "\n")
		new.pts <- lapply(temp.tl, new.points)
		xyp <- data.frame(new.pts)
		xyp <- data.matrix(xyp)
		tp.t <- runif(n.def*length(temp.tl), X$tcoord[1], X$tcoord[2])
		tp <- matrix(tp.t, ncol = length(temp.tl))
		new.p <- rbind(xyp, tp)
		if(length(temp.tl) > 1) {
			lamb2 <- rbind(lamb2, apply(new.p, 2, all))
			lamb2.ave <- apply(lamb2, 2, mean)
			lamb2.sd <- apply(lamb2, 2, sd)
		} else {
			lamb2 <- c(lamb2, apply(new.p, 2, all))
			lamb2.ave <- mean(lamb2)
			lamb2.sd <- sd(lamb2)
		}
		error <- lamb2.ave/100
		std.error <- lamb2.sd/sqrt(n)
		if(any(std.error < error)) {
			lose <- which(std.error < error)
			replace <- place[lose]
			lamb.ave.final[replace] <- lamb2.ave[lose]
			lamb.sd.final[replace] <- lamb2.sd[lose]
			n.final[replace] <- n
			temp.tl <- temp.tl[-lose]
			if(sum(std.error < error) < length(error)) {
				lamb2 <- lamb2[ , -lose]
			} 
			place <- place[-lose]
		} 
	}
	int.approx <- lamb.ave.final * vol	
	tess <- list(n = n.final, integral = int.approx, mean.lambda = lamb.ave.final, sd.lambda = lamb.sd.final)
	class(tess) <- "tess.info"
	return(tess)
}