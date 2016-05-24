tessdev.info <- function(X, cifunction1, cifunction2, theta1 = NULL, theta2 = NULL, areas, tl)
{
	data <- data.frame(cbind(X$x, X$y, X$t))
	names(data) <- c("x", "y", "t")
	if(is.null(theta1)) {
		lamb1 <- cifunction1(X)
		fun2.1 <- function(z)
		{
			Y <- rbind(data, z)
			Y <- Y[which(Y$t <= z[3]),]
			Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
			tail(cifunction1(Y),1)
		}
	}
	if(!is.null(theta1)) {
		lamb1 <- cifunction1(X, theta = theta1)
		fun2.1 <- function(z)
		{
			Y <- rbind(data, z)
			Y <- Y[which(Y$t <= z[3]),]
			Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
			tail(cifunction1(Y, theta1),1)
		}
	}
	if(is.null(theta2)) {
		lamb2 <- cifunction2(X)
		fun2.2 <- function(z)
		{
			Y <- rbind(data, z)
			Y <- Y[which(Y$t <= z[3]),]
			Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
			tail(cifunction2(Y),1)
		}
	}
	if(!is.null(theta2)) {
		lamb2 <- cifunction2(X, theta = theta2)
		fun2.2 <- function(z)
		{
			Y <- rbind(data, z)
			Y <- Y[which(Y$t <= z[3]),]
			Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
			tail(cifunction2(Y, theta2),1)
		}
	}
	all <- function(tl)
	{
		cells <- as.points(tl$x, tl$y)
		place <- inpip(all.p[ ,1:2], cells)
		lamb1.ave <- mean(lamb.t1[place])
		lamb1.sd <- sd(lamb.t1[place])
		lamb1.sd[which(is.na(lamb1.sd))] <- lamb1.ave[which(is.na(lamb1.sd))]
		lamb2.ave <- mean(lamb.t2[place])
		lamb2.sd <- sd(lamb.t2[place])
		lamb2.sd[which(is.na(lamb2.sd))] <- lamb2.ave[which(is.na(lamb2.sd))]
		n2 <- length(place)
		return(c(lamb1.ave, lamb1.sd, lamb2.ave, lamb2.sd, n2))
	}
	vol <- areas * diff(X$tcoord)
	e <- 0
	n <- 0
	lamb.t1 <- lamb1
	lamb.t2 <- lamb2
	all.p <- data
	while(e < (2*length(length(X$x)))) {
		n <- n + 10000
		cat("total number of points used to estimate integrals: ", n, "\n")
		xpoints <- runif(10000, X$xcoord[1], X$xcoord[2])
		ypoints	<- runif(10000, X$ycoord[1], X$ycoord[2])
		tpoints <- runif(10000, X$tcoord[1], X$tcoord[2])
		new.p <- data.frame(cbind(xpoints, ypoints, tpoints))
		names(new.p) <- c("x", "y", "t")
		lamb.temp1 <- apply(new.p, 1, fun2.1)
		lamb.temp2 <- apply(new.p, 1, fun2.2)
		lamb.t1 <- c(lamb.t1, lamb.temp1)
		lamb.t2 <- c(lamb.t2, lamb.temp2)
		all.p <- rbind(all.p, new.p)
		info <- lapply(tl, all)
		ave1 <- unlist(lapply(info, function(w){w[1]}))
		ave2 <- unlist(lapply(info, function(w){w[3]}))
		sd1 <- unlist(lapply(info, function(w){w[2]}))
		sd2 <- unlist(lapply(info, function(w){w[4]}))
		n3 <- unlist(lapply(info, function(w){w[3]}))
		e <- sum(sd1/sqrt(n3) < ave1/100) + sum(sd2/sqrt(n3) < ave2/100)
	}
	int1.approx <- ave1 * vol
	int2.approx <- ave2 * vol	
	tess <- list(n = n, integral1 = int1.approx, integral2 = int2.approx, mean.lambda1 = ave1, mean.lambda2 = ave2, sd.lambda1 = sd1, sd.lambda2 = sd2)
	class(tess) <- "tess.info"
	return(tess)
}