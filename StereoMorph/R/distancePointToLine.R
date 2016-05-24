distancePointToLine <- function(p, l1, l2 = NULL){

	# CONVERT SINGLE POINTS TO ONE-ROW MATRIX
	if(!is.matrix(p)) p <- matrix(p, nrow=1, ncol=length(p))

	# GET LINE PARAMETERS IF L1 IS LIST
	if(is.list(l1)){
		# GET LINE POINTS IF IN l1 LIST
		if('l1' %in% names(l1)){l2 <- l1$l2;l1 <- l1$l1}
	
		# GET LINE POINTS IF m, b
		if('m' %in% names(l1)){l2 <- c(1, l1$m + l1$b);l1 <- c(0, l1$b)}
	
		# GET LINE POINTS IF a, b, c
		if('a' %in% names(l1)){l2 <- c(1, -(l1$a + l1$c)/l1$b);l1 <- c(0, -l1$c/l1$b)}
	}

	# IF INPUT IS 2D POINT(S) AND l1, l2, ADD ZERO THIRD DIMENSION
	if(ncol(p) == 2){p <- cbind(p, rep(0, nrow(p)));l1 <- c(l1, 0);l2 <- c(l2, 0)}

	# FIND SHORTEST DISTANCE BETWEEN POINT(S) AND LINE
	d <- rep(NA, nrow(p))
	for(i in 1:nrow(p)) d[i] <- sqrt(sum(cprod_SM(p[i, ] - l1, p[i, ] - l2)^2)) / sqrt(sum((l2 - l1)^2))

	d
}