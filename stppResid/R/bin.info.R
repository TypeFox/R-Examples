bin.info <- function(X, cifunction, theta = NULL, gf, type = c(1,2), n.def = 100)
{
	xpoints <- function(x) 
	{
		y <- runif(n.def, x[1], x[2])
		return(y)	
	}
	ypoints <- function(x)
	{
		y <- runif(n.def, x[3], x[4])	
		return(y)
	}
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
	vol <- diff(X$xcoord) * diff(X$ycoord) * diff(X$tcoord) / nrow(gf$grid.full)
	data <- data.frame(cbind(X$x, X$y, X$t))
	names(data) <- c("x", "y", "t")
	n <- 0
	lamb2 <- c()
	temp.grid <- gf$grid.full
	place <- 1:nrow(gf$grid.full)
	lamb.ave.final <- rep(0, length(place))
	lamb.sd.final <- rep(0, length(place))
	n.final <- rep(0, length(place))
	n = n + n.def
	xp <- apply(temp.grid, 1, xpoints)
	yp <- apply(temp.grid, 1, ypoints)
	tp.t <- runif(n.def*nrow(temp.grid), X$tcoord[1], X$tcoord[2])
	tp <- matrix(tp.t, ncol = nrow(temp.grid))
	new.p <- rbind(xp, yp, tp)
	lamb2 <- rbind(lamb2, apply(new.p, 2, all))
	if(type == 1) {
		lamb2.ave <- apply(lamb2, 2, mean)
		lamb2.sd <- apply(lamb2, 2, sd)
	} else {
		lamb2.ave <- apply(sqrt(lamb2), 2, mean)
		lamb2.sd <- apply(sqrt(lamb2), 2, sd)
	}
	while(length(place) > 0) {
		n = n + n.def
		cat("number of bins left: ", nrow(temp.grid), "\n")
		xp <- apply(temp.grid, 1, xpoints)
		yp <- apply(temp.grid, 1, ypoints)
		tp.t <- runif(n.def*nrow(temp.grid), X$tcoord[1], X$tcoord[2])
		tp <- matrix(tp.t, ncol = nrow(temp.grid))
		new.p <- rbind(xp, yp, tp)
		if(nrow(temp.grid) > 1) {
			lamb2 <- rbind(lamb2, apply(new.p, 2, all))
			if(type == 1) {
				lamb2.ave <- apply(lamb2, 2, mean)
				lamb2.sd <- apply(lamb2, 2, sd)
			} else {
				lamb2.ave <- apply(sqrt(lamb2), 2, mean)
				lamb2.sd <- apply(sqrt(lamb2), 2, sd)
			}
		} else {
			lamb2 <- c(lamb2, apply(new.p, 2, all))
			if(type == 1) {
				lamb2.ave <- mean(lamb2)
				lamb2.sd <- sd(lamb2)
			} else {
				lamb2.ave <- mean(sqrt(lamb2))
				lamb2.sd <- sd(sqrt(lamb2))
			}
		}
		error <- lamb2.ave/100
		std.error <- lamb2.sd/sqrt(n)
		if(any(std.error < error)) {
			lose <- which(std.error < error)
			replace <- place[lose]
			lamb.ave.final[replace] <- lamb2.ave[lose]
			lamb.sd.final[replace] <- lamb2.sd[lose]
			n.final[replace] <- n
			temp.grid <- temp.grid[-lose, ]
			if(sum(std.error < error) < length(error)) {
				lamb2 <- lamb2[ , -lose]
			} 
			place <- place[-lose]
		} 
	}
	int.approx <- lamb.ave.final * vol
	bins <- list(n = n.final, integral = int.approx, mean.lambda = lamb.ave.final, sd.lambda = lamb.sd.final)
	class(bins) <- "bin.info"
	return(bins)
}
