# validate.cm validates a community matrix by testing that it:
#   is a square matrix, 
#   has only elements of -1, 0 or 1, 
#   has at least two parameters,
#   with at least one direct or indirect path from each variable to all other variables, and
#   isnot fully specified. 
# It produces an error if the community matrix is 
# invalid by these tests, and returns nothing otherwise. It takes:
# CM: a potential community matrix
validate.cm <- function(CM) {

	is.isolated <- function(CM) {
	
		# one.path returns the first path from a depth-first search. It takes:
		#   CM: a Community Matrix
		#   i:  a starting parameter
		#   j:  an ending parameter
		one.path <- function(CM,i,j) {

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
					if (!(identical((CM[LOVE[length(LOVE)],x]),0))) {
						if (identical(Term[LOVE[length(LOVE)],x],0)) {
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
						Term[x,c(1:N)] <<- 0
						}
					}
				# end update.Term()
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
					Term[LOVE[length(LOVE)-1],LOVE[length(LOVE)]] <<- 1
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
							if (Term[i,val]==TRUE) {
								test <- TRUE
								}
							}
						if (test) {
							return()
							}
						Term[LOVE[length(LOVE)],1:N] <<- 1
						LOVE <<- LOVE[-length(LOVE)]
					   update.Term()
						incomplete <<- FALSE
						}
					} else {
							# append the first element (breadthwise search) of ENVY to LOVE
				 			LOVE <<- c(LOVE,ENVY[1])
							# test whether new LOVE is a path from i to j and respond accordingly
								if (is.path(LOVE,j) == TRUE) {
								LOP <<- list(LOVE)
								incomplete <<- FALSE
								}
							}
	
			# end search.step()
			}
	
		# Set N = rowsize of CM
		N <- nrow(CM)

		# initialize search termination matrix
		Term <- matrix(c(0),N,N)

		 # initialize list of paths (LOP), list of visited elements (LOVE)
		 LOP <- NULL
		 LOVE <- c(i)
	 
		incomplete <- TRUE
		while (incomplete) {
			search.step()
			}
 	
 		return(LOP[[1]])
 	
 		# end one.path()
	 	}
	
		N <- nrow(CM)
	
		#list path from 1 to elements 2 to N
		reached.elements <- NULL
		reaching.elements <- NULL
		for (x in 2:N) {
			if (!(x %in% reached.elements)) {
				path <- one.path(CM,1,x)
				reached.elements <- c(reached.elements, path[-1])
				}
			if (!(x %in% reaching.elements)) {
				path <- one.path(CM,x,1)
				reaching.elements <- c(reaching.elements, path[-length(path)])
				}
			}

		return(!(length(unique(reached.elements)) == N-1) || !(length(unique(reaching.elements)) == N-1) ) 	
		# end is.isolated()
		}


	# Is CM a matrix? a square matrix?
 	if (!(is.matrix(CM)) | !(identical( nrow(CM), ncol(CM) ) ) ) {
	 	stop("\nA Community Matrix must be a square matrix with elements \nof values of only 1, 0 and -1.")
	 	}

	# Is CM big enough?
	if (nrow(CM) == 1) {
	 	stop("\nA Community Matrix must have two or more parameters.")
	 	}

	# Does CM contain only values = 1, 0 or -1?
	for (i in 1:nrow(CM)) {
		for (j in 1:ncol(CM)) {
	 		if ( !( (CM[i,j] == 1) | (CM[i,j] == 0) | (CM[i,j] == -1) ) ) {
	 			stop("\nA Community Matrix must be a square matrix with elements \nof values of only 1, 0 and -1.")
	 			}
	 		}
	 	}

	# Is CM matrix is fully specified?
	if (!(0 %in% CM)) {
		stop("\nA fully connected Community Matrix is not analyzable via loop analysis.")
		}
		
	N <- nrow(CM)

	# Is a subsystem of the Community Matrix isolated?
	if (is.isolated(CM)) {
		stop("\nTwo or more variables in the system are isolated from one another.")
		}

	# end validate.cm
	}
