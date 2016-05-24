# make.cm() interactively solicits variable names and linakges between
# variables to produce a community matrix. It takes:
# N: the number of variables in the community matrix

make.cm <- function(n=NA) {

	N <- n
	
	while (!is.numeric(N) || !identical(floor(N),N) || N < 2) {
		N <- readline("\nPlease specify how many parameters are in the system represented by the \ncommunity matrix by typing an integer greater than or equal to two (or type \n\"q\" to quit): ")
		if (N == "q" || N == "Q") {
			return(cat("\n"))
			}
		N <- as.numeric(N)
		if (N > 12) {
			warning("\nmake.cem will take a very long time to compute community matrices larger than 12 parameters.")
			}
		}

 	CommunityMatrix <- matrix(0,N,N,dimnames = list(c(letters[1:N]), c(letters[1:N])))

	yes.name <- readline("\nWould you like to name your parameters? [y/n]\nThe default is lettered names a, b, c, etc.\n(short names of a few characters will output best) ")
	
	if ( identical(substr(yes.name,1,1),"y") || identical(substr(yes.name,1,1),"Y") ) {
		for (x in 1:N) {
			name <- readline(sprintf("\nWhat is parameter %i\'s name? ", x))
			rownames(CommunityMatrix)[x] <- name
			colnames(CommunityMatrix)[x] <- name
			}
		}
	cat("\n")

	# input qualitative causal relationships between each parameter
	for (i in 1:N) {
		for (j in 1:N) {
			i <- as.integer(i)
			j <- as.integer(j)
			while (N) {
				Nij <- readline(sprintf("Relation from %s to %s (\"q\" to quit): ", rownames(CommunityMatrix)[j], rownames(CommunityMatrix)[i]))
				if (identical(Nij,"1") || identical(Nij,"+")) {
					CommunityMatrix[i,j] <- 1
					break
					}
				if (identical(Nij,"-1") || identical(Nij,"-")) {
					CommunityMatrix[i,j] <- -1
					break
					}
				if (identical(Nij,"0") || identical(Nij," ") || identical(Nij,"")) {
						CommunityMatrix[i,j] <- 0
					break
					}
				if (Nij == "q" || Nij == "Q") {
					return(cat("\n"))
					}
				# alert if input isn't in accepted form, and rerequest
				if ( !(Nij %in% c("1","-1","0","+","-",""," ") ) ) {
					cat(sprintf("Relation from %s to %s: must be:\n increasing:  1, or +\n no relation: 0, <space> or <return>\n decreasing:  -1, or 1", letters[i], letters[as.integer(j)]))
					}
				}
			}
		}

	# alert if matrix is fully specified.
	if (!(0 %in% CommunityMatrix)) {
		stop("\nA fully connected Community Matrix is not analyzable via loop analysis.")
		}
		
	# alert if a paramter has no parents or children, and contract system.
	cmcollapse <- function(CM) {
		for (x in 1:N) {
			if (
			   (!(-1 %in% CM[1:N,x]) & 
			    !( 1 %in% CM[1:N,x]) ) || 
			   (!(-1 %in% CM[x,1:N]) & 
			    !( 1 %in% CM[x,1:N]) ) ) {
				CM <- CM[-x,-x]
				warning(sprintf("\nParameter %i has no parent or child parameters and has \nbeen removed from the Community Matrix.\n", x))
				N <<- N - 1
				if (identical(N,1)) {
					stop("The system is too small after removing childless/parentless parameters!")
					}
				cmcollapse(CM)
				}
			}

		return(CM)

		# end cmcollapse()
		}

	CommunityMatrix <- cmcollapse(CommunityMatrix)
	
return(CommunityMatrix)

}
