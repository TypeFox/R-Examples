Rclusterpp.setThreads <- function(threads=1) {
	threads <- ifelse(is.null(threads), .Call("rclusterpp_get_num_procs", PACKAGE="Rclusterpp"), threads)
	invisible(.Call("rclusterpp_set_num_threads", threads=as.integer(threads), NAOK=FALSE, PACKAGE="Rclusterpp"))
}
	

Rclusterpp.linkageKinds <- function() {
	.Call("linkage_kinds", PACKAGE="Rclusterpp")
}

Rclusterpp.distanceKinds <- function() {
	.Call("distance_kinds", PACKAGE="Rclusterpp")
}

Rclusterpp.hclust <- function(x, method="ward", members=NULL, distance="euclidean", p=2) {
	METHODS <- Rclusterpp.linkageKinds()
	method  <- pmatch(method, METHODS)
	if (is.na(method))
		stop("Invalid clustering method")
  if (method == -1) 
    stop("Ambiguous clustering method")

	if (class(x) == "dist") {
		dist.method = attributes(x)$method
		labels      = attributes(x)$Labels

		hcl <- .Call("hclust_from_distance", 
								 data = as.double(x),
								 size = as.integer(attributes(x)$Size),
								 link = as.integer(method), 
								 DUP = FALSE, NAOK = FALSE, PACKAGE = "Rclusterpp" )
	
		hcl$labels      = labels 
		hcl$method      = METHODS[method]
		hcl$call        = match.call()
		hcl$dist.method = dist.method 
		class(hcl) <- "hclust"
	
		return(hcl)
	} else {
		if (!is.null(members)) {
			stop("members must be null when clustering from data")
		}

		DISTANCES <- Rclusterpp.distanceKinds()
		distance  <- pmatch(distance, DISTANCES)
		if (is.na(distance))
			stop("Invalid distance metric")
		if (method == -1)
			stop("Ambiguous distance metric")

		if (METHODS[method] == "ward" && DISTANCES[distance] != "euclidean") {
			warning("Distance method is forced to (squared) 'euclidean' distance for Ward's method")
			distance <- which(DISTANCES == "euclidean")[1]
		}
	
		N <- nrow(x <- as.matrix(x))
		hcl <- .Call("hclust_from_data", 
		             data = x,
								 link = as.integer(method), 
								 dist = as.integer(distance),
								 p    = as.numeric(p),
								 DUP = FALSE, NAOK = FALSE, PACKAGE = "Rclusterpp" )
		
		hcl$labels = row.names(x)
		hcl$method = METHODS[method]
		hcl$call   = match.call()
		hcl$dist.method = DISTANCES[distance]
		class(hcl) <- "hclust"
		
		return(hcl)
	}
}

