pointsAtEvenSpacing <- function(x, n){

	# IF X IS LIST, TAKE ONLY FIRST ELEMENT
	if(is.list(x)) x <- x[[1]]

	# GET ORIGINAL ROW NAMES
	row_names <- rownames(x)

	# REMOVE NA VALUES
	x <- x[rowSums(is.na(x)) == 0, ]

	# ERROR IF THE NUMBER OF NEW POINTS EXCEEDS THE NUMBER OF INPUT POINTS
	#if(n > nrow(x)) stop(paste0("The number of new points (", n, ") exceeds the number of non-NA input points (", nrow(x), ")."))

	# EMPTY MATRIX FOR NEW POINTS
	r <- matrix(NA, nrow=n, ncol=ncol(x))

	# SET STARTING AND END POINTS
	r[1, ] <- x[1, ]
	r[nrow(r), ] <- x[nrow(x), ]

	# FIND INTERPOINT DISTANCES IN X
	if(nrow(x) == 2){
		d <- sqrt(sum((x[2, ] - x[1, ])^2))
	}else{
		d <- sqrt(rowSums((x[2:nrow(x), ] - x[1:(nrow(x)-1), ])^2))
	}

	# FIND CUMULATIVE DISTANCE ALONG POINTS
	cd <- c(0, cumsum(d))

	# GET TOTAL DISTANCE
	td <- cd[length(cd)]

	# FIND SPACING BETWEEN POINTS
	s <- td / (n-1)
	
	# SET STARTING POINT FOR SEARCH
	under_pt <- 2

	for(i in 1:(n-2)){

		# IF UNDER POINT GOES DOWN TO 1, SET TO 2
		if(under_pt == 1) under_pt <- 2

		# FIND CUMULATIVE DISTANCES LESS OR EQUAL TO NEXT INTERVAL LIMIT - ONLY SEARCH IN POSSIBLE RANGE
		le_int <- which(cd[(under_pt-1):length(cd)] <= s*i) + under_pt - 2

		# GET LAST INDEX - AT OR JUST BEFORE POINT AT SPACING
		under_pt <- le_int[length(le_int)]
		
		# CREATE VECTOR BETWEEN UNDER AND NEXT POINT
		v <- x[under_pt+1, ] - x[under_pt, ]
		
		# SCALE TO UNIT VECTOR		
		v <- v*(1 / (cd[under_pt+1] - cd[under_pt]))

		# SCALE VECTOR TO FINAL PORTION OF INTERVAL
		v <- v*(s*i - cd[under_pt])
		
		# ADD VECTOR TO UNDER POINT POSITION TO GET NEXT POINT
		r[i+1, ] <- x[under_pt, ] + v
	}
	
	# ASSIGN OLD ROWNAMES TO NEW MATRIX IF NUMBER OF NEW POINTS IS LESS THAN ORIGINAL NUMBER OF POINTS
	if(n < nrow(x)) rownames(r) <- row_names[1:nrow(r)]

	r
}