# feedback returns the adjusted sum of products of SOSLs. It takes:
# C: the complimentary subsystem of a path through a CM

feedback <- function(C) {

# validate.c validates a matrix by testing that it is
# square, and has only elements of -1, 0 or 1 It produces an error 
# if the community matrix is invalid by these tests, and returns 
# nothing otherwise. It takes:
# C: a matrix

	validate.c <- function(C) {

		# Is C a matrix? a square matrix?
 		if (!(is.matrix(C)) | !(identical( nrow(C), ncol(C) ) ) ) {
		 	stop("\nThe system must be specified by a square matrix with elements \nof values of only 1, 0 and -1.")
	 		}

		# Does C contain only values = 1, 0 or -1?
		rcC <- as.vector(C)
		indexend <- length(rcC)
		for (i in 1:indexend) {
			if ( !( (rcC[i] == 1) | (rcC[i] == 0) | (rcC[i] == -1) ) ) {
				stop("\nThe system must be a square matrix with elements \nof values of only 1, 0 and -1.")
				}
			}
		# end validate.c
		}

# sosl.prod returns the sign product of a set of spanning loops
# It takes:
# SOSL: a single set of spanning loop(s)

	sosl.prod <- function(C,SOSL) {

		loop.prod <- function(C,loop) {
			lprod <- 1
			if (length(loop) >1 ) {
				for (edge in 1:(length(loop))) {
					if (edge < length(loop)) {
						lprod <- lprod * C[loop[edge],loop[edge+1]]
						} else {
							lprod <- lprod * C[loop[edge],loop[1]]
						 	}
					 }
				}
			if (length(loop) == 1) {
				lprod <- lprod*C[loop[1],loop[1]]
				}
			return(lprod)
				
			# end loop.prod()
			}
	
		sprod <- 1
		for (loop in SOSL) {
			sprod <- sprod*loop.prod(C,loop)
			}
		return(sprod)
	
		# end sosl.prod()
		}


# make.MOSL returns a matrix of searchable loops (MOSL) a ragged 3-D data
# structure that can be understood as an upper reverse diagonal N by N matrix, 
# where the first dimension is the starting parameter of a loop, the second is
# the length of the loop, and the third is the set of all the passed list of 
# loops. For example:
#
#   make.MOSL(LOL,3)[[2]][[3]][[1]]
# 
# would produce the first loop of length 3 beginning with parameter 2.
# moreover:
#
#   length(make.MOSL(LOL,3)[[2]][[3]])
#
# would count the number of loops of length 3 beginning with parameter 2. 
# make.MOSL() takes:
# LOL: a list of loops
# N: a scalar number of parameters in the system

	make.MOSL <- function(LOL,N) {

		MOSL <- rep(list(NA),N)

		for (j in 1:N) {
			MOSL[[j]]<-rep(list(NA),(N-(j-1)))
			}

		for (a in 1:length(LOL)) {
			if (identical((MOSL[[ LOL[[a]][1] ]] [[ length(unique(LOL[[a]])) ]]),NA)) {
				MOSL[[ LOL[[a]][1] ]] [[ length(unique(LOL[[a]])) ]] <- list(as.double(unique(LOL[[a]])))
				} else {
				MOSL[[ LOL[[a]][1] ]] [[ length(unique(LOL[[a]])) ]] <- c(MOSL[[ LOL[[a]][1] ]] [[ length(as.double(unique(LOL[[a]]))) ]], list(unique(LOL[[a]])))
				}
			}

		return(MOSL)

		# end make.MOSL()
		}


# enumerate.SOSL returns all sets of spanning loops (SOSL). It takes:
# MOSL: a 3-D ragged matrix of searchable of loops
# N: the number of parameters (some potential rows in MOSL may be empty)

enumerate.SOSL <- function(MOSL,N) {

# set.size returns the number of parameters in PLOS. It takes:
# PLOS: potential list of sets

	set.size <- function(PLOS) {
  	size <- 0
  	if (length(PLOS)==0) {
		return(size)
		}

	for (q in 1:length(PLOS)) {
		size <- size+length(PLOS[[q]])
		}

	return(size)
		
	# end set.size()
	}


# make.loopENVY returns a acceptable starting paramenters (search.row) of the 
# next acceptable loop in MOSL based on which parameters remain to be searched. 
# It takes:
# PLOS: a list of loops in N
# N: a scalar indicating the number of parameters in the system.

		make.loopENVY <- function(PLOS, MOSL) {

			if (length(PLOS) == 0) {
				return(seq(1:length(MOSL)))
			  	}
	
			max.search.space <- seq(1:length(MOSL))

			search.row <- max.search.space
  			for (x in 1:length(PLOS)) {
				search.row <- setdiff(search.row,PLOS[[x]])
				}

			return(search.row)

			# end make.loopENVY
			}


# initialize.term returns a termination data structure for MOSL.
		initialize.term <- function(MOSL) {
	
			N1 <- length(MOSL)
	
			Term <- MOSL
			for (i in 1:length(MOSL)) {
				for (j in 1:length(MOSL[[i]])) {
					Term[[i]][[j]] <- c(rep(0,length(MOSL[[i]][[j]])))
					}
				}
			
			return(Term)
	
			# end initialize.term()
			}

	   N.mosl <- length(MOSL)
		Term <- initialize.term(MOSL)
		PLOS <- NULL
		SOSL <- NULL
		k.last <- NULL

		search.row <- function(PLOS) {
 
			if (length(PLOS) == 0) {    				return(1)
				}
 
			row <- seq(1:N.mosl)
			for (x in 1:length(PLOS)) {
				row <- setdiff(row,PLOS[[x]])
				}

			return(row[[1]])
     
			# end search.row
			}  
 
		search.over <- function(row) {
			if (row == 1) {
				for (j in 1:length(MOSL[[1]])) {
					for (k in 1:length(MOSL[[1]][[j]])) {
						if (Term[[1]][[j]][[k]] == 0) {
							return(FALSE)
							}
						}
					}
				return(TRUE)
				}
 
			return(FALSE)
     
			# end search.over()
			}
 
		next.loop <- function(row) {

			for (j in 1:length(MOSL[[row]])) {
				for (k in 1:length(MOSL[[row]][[j]])) {
					if ((Term[[row]][[j]][[k]] == 0) & (!is.na(MOSL[[row]][[j]][[k]][[1]])) & (length(MOSL[[row]][[j]][[k]])<=N.mosl-set.size(PLOS))) {
						return(list(MOSL[[row]][[j]][[k]],k))
						}
					}
				}
     
			return(list(NULL,NULL))
     
			# end next.loop()
			}
      
		while (TRUE) {

			row <- search.row(PLOS)
     
			# is the search over?
			if (search.over(row)) {
				return(SOSL)
				}
     
			# if there's no valid search row...
			if (length(row) == 0) {
				return(SOSL)    
				}
   
			loop <- next.loop(row)

			# if there's no valid loop in the row and it is the first...
			if (is.null(loop[[1]]) & row == 1) {
				return(SOSL)
				} 
   
			# if there's no valid loop in the row and it's not the first...
			if (is.null(loop[[1]]) & row != 1) {
 
				# clear remaining loops
				for (i in make.loopENVY(PLOS, MOSL)) {
					for (j in 1:length(MOSL[[i]])) {
						for (k in 1:length(MOSL[[i]][[j]])) {
							Term[[i]][[j]][[k]] <- 0
							}
						}
					}
        
				# clear the row in Term
				i <- row
				j <- length(MOSL[[row]])
				for (k in 1:length(MOSL[[i]][[j]])) {
					Term[[i]][[j]][[k]] <- 0
					}
     
				# terminate the last loop in PLOS
				i <- PLOS[[length(PLOS)]][[1]]
				j <- length(PLOS[[length(PLOS)]])
				k <- loop[[2]]
				for (jj in 1:j) {
					for (kk in 1:length(Term[[i]][[jj]])) {
						if (jj == j & kk <= k.last[[length(k.last)]]) {
							Term[[i]][[jj]][[kk]] <- 1
							}
						}
					}
 
				# remove the last loop from PLOS
				PLOS <- PLOS[1:length(PLOS)-1]
				k.last <- k.last[1:length(k.last)-1]
				next
				}
 
			# add next.loop to PLOS
			if (is.null(PLOS)) {
				PLOS <- list(loop[[1]])
				k.last <- c(k.last, list(loop[[2]]))
				} else {
						PLOS <- c(PLOS,list(loop[[1]]))
						k.last <- c(k.last, list(loop[[2]]))
						}
	     
			# test if PLOS spans N and add to SOSL if it does
			if (set.size(PLOS) == N) {
 
				# add PLOS to SOSL
				if (is.null(SOSL)) {
					SOSL <- list(PLOS)
					} else {
							SOSL <- c(SOSL, list(PLOS))
							}

				# terminate the last loop in PLOS
				i <- PLOS[[length(PLOS)]][[1]]
				j <- length(PLOS[[length(PLOS)]])
				k <- loop[[2]]
				Term[[i]][[j]][[k.last[[length(k.last)]]]] <- 1
 
				# remove the last loop from PLOS
				PLOS <- PLOS[1:length(PLOS)-1]
				k.last <- k.last[1:length(k.last)-1]
				next
				}
			# end while TRUE    
			}

		# end enumerate.SOSL()
		}

	validate.c(C)

	N <- nrow(C)
	if (is.null(nrow(C))) {
		return(C)
		}

   LOL <- (enumerate.loops(C))
   if (is.null(LOL)) {
   		return(0)
   		}
	Sum <- 0
	for (SOSL in enumerate.SOSL(make.MOSL(LOL,N),N) ) {
		adjust <- (-1)^(length(SOSL)+1)
		sprod <- sosl.prod(C,SOSL)*adjust
		if (Sum == 0) {
			Sum <- sprod
			} else {
			if (Sum == -1*sprod) {
				return(NA)
				}
			}
		}

	return(Sum)
	
	# end feedback()
	}
	
	
