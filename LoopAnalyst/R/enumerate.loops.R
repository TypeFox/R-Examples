# enumerate.loops returns a list of simple loops (cycles) in CM. It 
# takes as arguments:
# CM: a (-1, 0, 1) matrix
enumerate.loops <- function(CM) {

# is.loop tests that a list of visited elements begins at and 
# terminates in the ith parameter. It returns TRUE or FALSE
	is.loop <- function() {
		if ( (LOVE[length(LOVE)] == LOVE[1]) & (length(LOVE) > 1) ) {
			return(TRUE)
			} else {
		 	return(FALSE)
		 	}
	 		
		# end is.loop()
		} 


# make.ENVY() returns a list of Elements Not Visted Yet (ENVY), excluding
# the first element from:
	make.ENVY <- function() {
		ENVY <- NULL
		for (z in 1:N) {
			if (!(identical((CM[z,LOVE[length(LOVE)]]),0))) {
				if (identical(Term[z,LOVE[length(LOVE)]],0)) {
					ENVY <- c(ENVY,z)
					}
				}
			}
		for (z in LOVE[-1]) {

			ENVY <- ENVY[ENVY != z]
			}

		return(ENVY)
				
		# end make.ENVY()	
		}

	# update.Term() clears rows in the termination matrix Term of any
	# parameters not in the List Of Visited Elements (LOVE).
	update.Term <- function() {
		for (y in 1:N) {
			if (!(y %in% LOVE) ) {
				Term[c(1:N),y] <<- 0
				}
			}

		Term[1:N,1:x-1] <<- 1
		Term[1:x-1,1:N] <<- 1
			
		# end update.Term()
		}

	# add.LOVE.to.LOL adds the List Of Visited Elements (LOVE) to the
	# List Of Loops (LOL)
	add.LOVE.to.LOL <- function() {
		if (is.null(LOL)) {
			LOL <<- list(LOVE)
			} else {
		 	LOL <<- c(LOL,list(LOVE))
		 	}
		 	
		 # end add.LOVE.to.LOL()
		 }

	# search.step is the main logic of the breadth-wide path search.
	search.step <- function() {
		
		ENVY <- make.ENVY()
		if (length(ENVY) == 0) {
			Term[LOVE[length(LOVE)],LOVE[length(LOVE)-1]] <<- 1
			LOVE <<- LOVE[-length(LOVE)]
		   update.Term()

			# exit SearchStep if the last element of LOVE is i or if LOVE is empty 
			# and return List Of Loops (LOL)
			if (length(LOVE) == 0 ) {
				incomplete <<- FALSE
				return()
				}
			} else {

			# append the first element (breadthwise search) of ENVY to LOVE
 			LOVE <<- c(LOVE,ENVY[1])	

			# test whether new LOVE is a loop and respond accordingly
			if (is.loop() == TRUE) {
				add.LOVE.to.LOL()
				Term[LOVE[length(LOVE)],LOVE[length(LOVE)-1]] <<- 1
				LOVE <<- LOVE[-length(LOVE)]
				}
 			}
	
		# end search.step()
		}

#	validate.cm(CM)

	# Set N = rowsize of CM
	N <- nrow(CM)
	# initialize search termination matrix
	Term <- matrix(c(0),N,N)
   # initialize list of loop (LOL), and list of visited elements (LOVE)
   LOL <- NULL
  
	for (x in 1:(N-1)) {
		LOVE <- c(x)
		Term[1:x-1,1:N] <- 1
		Term[1:N,1:x-1] <- 1

		# find all loops starting at x
		incomplete <- TRUE
		while (incomplete) {
			search.step()
			}
		}

	if ( !identical(CM[N,N],0) ) {
		LOL <- c( LOL,list( c(N,N) ) )
		}


	return(LOL)
	
	# end enumerate.loops()
	}
