# enumerate.paths returns a list of paths (LOP) takes as it's arguments:
#   CM: a Community Matrix
#   i:  a starting parameter
#   j:  an ending parameter
enumerate.paths <- function(CM,i,j) {

	if (!is.numeric(i)) {
		if (is.na(match(i,rownames(CM)))) {
			stop("\nvalue \"", i, "\" of i is not the name of a parameter in the community matrix!\n")
			}
		match(i, rownames(CM)) -> i
		}
		
	if (!is.numeric(j)) {
		if (is.na(match(j,rownames(CM)))) {
			stop("\nvalue \"", j, "\" of j is not the name of a parameter in the community matrix!\n")
			}
		match(j, rownames(CM)) -> j
		}

	validate.i.j <- function(i,j,N) {
		if ( !(floor(i)==i) | !(floor(j)==j) | i > N | j > N | i < 1 | j < 1) {
			stop("\ni and j must be integer values between 1 and ", N, ", or the names of the \ncommunity matrix parameters!")
			}
		# end validate.i.j
		}

 
# is.path tests that a list of visited elements terminates in the jth
# parameter. It returns TRUE or FALSE and takes:
# LOVE: a list of visited elements
# j: a parameter that terminates a path from i to j

	is.path <- function(LOVE,j) {
		if (LOVE[length(LOVE)] == j) {
			return(TRUE)
			} else {
		 	return(FALSE)
		 	}
		 	
		# end is.path()
		} 


# make.ENVY() returns a list of Elements Not Visted Yet (ENVY) from:
# CM: a valid Community Matrix
# Term: a valid Termination matrix
# LOVE: a list of visited elements within CM 
	make.ENVY <- function() {

		ENVY <- NULL
		for (x in 1:N) {
			if (!(identical((CM[x,LOVE[length(LOVE)]]),0))) {
				if (identical(Term[x,LOVE[length(LOVE)]],0)) {
					ENVY <- c(ENVY,x)
					}
				}
			}
		for (x in LOVE) {
			ENVY <- ENVY[ENVY != x]
			}

		return(ENVY)
				
		# end make.ENVY()	
		}
		
	# update.Term() clears rows in the termination matrix Term of any
	# parameters not in the List Of Visited Elements (LOVE). It
	# modifies global Term:
	update.Term <- function() {
		for (x in 1:N) {
			if (!(x %in% LOVE)) {
				Term[c(1:N),x] <<- 0
				}
			}
		# end update.Term()
		}
		
	# add.LOVE.to.LOP adds the List Of Visited Elements (LOVE) to the
	# List Of Paths (LOP) and returns LOP. It takes:
	# LOVE: a list of visited elements
	# LOP: a list of paths
	add.LOVE.to.LOP <- function() {
		if (is.null(LOP)) {
			LOP <<- list(LOVE)
			} else {
		 			LOP <<- c(LOP,list(LOVE))
		 			}
		 # end add.LOVE.to.LOP()
		 }

	# search.step is the main logic of the breadth-wide path search. It
	# returns a list of paths (LOP), and takes:
	# CM: a valid Community Matrix
	# Term: a valid Termination matrix
	# LOVE: a list of visited elements within CM
	# LOP: a list of paths within CM 
	# j: a parameter within CM that terminates a path from i to j
	search.step <- function() {
		ENVY <- make.ENVY()
		# when there are no unvisited elements & LOVE has more than 1 element, 
		# terminate the last element of LOVE for the second to last element of 
		# LOVE.
		if (length(ENVY) == 0) {
			Term[LOVE[length(LOVE)],LOVE[length(LOVE)-1]] <<- 1
			LOVE <<- LOVE[-length(LOVE)]
		   update.Term()
			# exit SearchStep if the last element of LOVE is i & i is terminated
			# or if LOVE is empty and return List Of Paths (LOP)
			if (length(LOVE) == 0 ) {
				incomplete <<- FALSE
				}
			if (length(LOVE) == 1) {
				test <- FALSE
				for (val in c(1:N)[-i]) {
					if (Term[val,i]==TRUE) {
						test <- TRUE
						}
					}
				if (test) {
					return()
					}
				Term[1:N,LOVE[length(LOVE)]] <<- 1
				LOVE <<- LOVE[-length(LOVE)]
			   update.Term()
				incomplete <<- FALSE
				}
			} else {
					# append the first element (breadthwise search) of ENVY to LOVE
		 			LOVE <<- c(LOVE,ENVY[1])
					# test whether new LOVE is a path from i to j and respond accordingly
					if (is.path(LOVE,j) == TRUE) {
						add.LOVE.to.LOP()
						Term[LOVE[length(LOVE)],LOVE[length(LOVE)-1]] <<- 1
						LOVE <<- LOVE[-length(LOVE)]
						}
						}
	
		# end search.step()
		}
	
	# Set N = rowsize of CM
	N <- nrow(CM)

	validate.i.j(i,j,N)

	 # take care of the simple case of i = j
	 if (identical(i,j)) {
	 		LOP <- list(c(i,j))
	 		return(LOP)
	   	}

	# initialize search termination matrix
	Term <- matrix(c(0),N,N)

	 # initialize list of paths (LOP), list of visited elements (LOVE)
	 LOP <- NULL
	 LOVE <- c(i)
	 
	incomplete <- TRUE
	while (incomplete) {
		search.step()
		}
 	
 	return(LOP)
 	
 	# end enumerate.paths()
 	}
