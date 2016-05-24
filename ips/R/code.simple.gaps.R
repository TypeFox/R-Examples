code.simple.gaps <- function(x, append = TRUE){
	
	ncha <- dim(x)[2]
	ntax <- dim(x)[1]
	hh <- vector()
	for (i in 1:ncha){
		if ("-" %in% as.character(x[,i]))
		hh <- c(hh, i)
	}
	
	# list of gaps
	L <- list(length = 100)
	i <- 1
	while (length(hh) >= 1){
	
	gap <- hh[1]
	hh <- hh[-1]
	while (hh[1] - tail(gap, n = 1) == 1 && length(hh) > 0){
		gap <- c(gap, hh[1])
		hh <- hh[-1]
	}
	L[[i]] <- gap
	i <- i + 1
	}
	
	simple.gaps <- vector()
	gap.content <- matrix(nrow = ntax)
	for (j in seq(along = L)){
		gap <- L[[j]]
		
		if (length(gap) == 1){
			simple.gaps <- c(simple.gaps, j)
			ff <- as.character(x[, gap])
			ff[ff != "-"] <- "a"
			ff[ff == "-"] <- "g"
			gap.content <- cbind(gap.content, ff)
		}
		else {
			ff <- as.character(x[, gap])
			ff[ff != "-"] <- "a"
			ff[ff == "-"] <- "g"
			ff <- t(ff)
			colnames(ff) <- NULL
			ff <- unique(ff)
			if (dim(ff)[1] == 1){
				simple.gaps <- c(simple.gaps, j)
				gap.content <- cbind(gap.content, t(ff))
			}
		}	
	}
	gap.content <- gap.content[, -1]
	gap.content <- as.alignment(gap.content)
	gap.content$nam <- rownames(x)
	gap.content <- as.DNAbin(gap.content)
	if (append){
		L <-L[simple.gaps]
	nsg <- length(L)
	names(L) <- paste("Gap", 1:nsg, sep = "_")
	simple.gaps <- unlist(L)
	x <- x[, -simple.gaps]
	x <- cbind(x, gap.content)
	cat("\n", nsg, " simple gaps have been coded.\n", 		sep = "")
	cat("\nGap positions:")
	cat("\n--------------\n")
	print(L)
	}												else {
		binary <- as.character(gap.content)
		binary[binary == "a"] <- 0
		binary[binary == "g"] <- 1
		x <- binary
	}
	x
}