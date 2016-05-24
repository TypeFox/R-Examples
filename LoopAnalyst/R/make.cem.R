# make.cem returns a community effect matrix. It takes:
# CM: a valid community effect matrix

make.cem <- function(CM, status=FALSE, out=FALSE) {

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


# path.compliment takes as its input a path vector (Path) from i to j and returns the 
# complimentary subsystem within the Path (C.ij) assuming there exists an N of the CM
# in its calling environment.

	path.compliment <- function(Path,N) {
	
		C.ij <- setdiff(seq(1:N),Path)
		return(C.ij)
	
		# end path.compliment()
		}


# looplist.elements returns a list of all variables appearing in
# loops in LISTSET. It takes:
# LOOPLIST: a list of loops, for example PLOS.

	looplist.elements <- function(LOOPLIST) {
		if (length(LOOPLIST) == 0) {
			return(NULL)
			}
		elements <- NULL
		for (q in 1:length(LOOPLIST)) {
			elements <- c(elements,LOOPLIST[[q]])
			}
		return(unique(elements))
	
		# end looplist.elements
		}


# sum.path.x.C returns a sign sum of paths times the feedback of their 
# complimentary subsystems. It takes:
# CM: a community matrix
# i: the starting parameter
# j: the ending parameter
# N: the size of CM 

	sum.path.x.C <- function(CM,i,j,N) {


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


# enumerate.SOSL returns all sets of spanning loops (SOSL). It takes:
# MOSL: a 3-D ragged matrix of searchable of loops
# N: the number of parameters (some potential rows in MOSL may be empty)

enumerate.SOSL <- function(MOSL,N) {


# initialize.term returns a termination data structure for MOSL.
		initialize.term <- function(MOSL) {
	
			N1 <- length(MOSL)
	
			Term <- MOSL
			for (i in 1:N1) {
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
 
		# next.loop searches each list k of each column of a row in MOSL. If the list is:
		# 1. not terminated, 
		# 2. not NA, 
		# 3. less than or equal to the length of N - PLOS, and
		# 4. candidate loop does not contain elements in PLOS, then
		# next.loop returns that list.
		next.loop <- function(row) {

			for (j in 1:length(MOSL[[row]])) {
				for (k in 1:length(MOSL[[row]][[j]])) {
					if ((Term[[row]][[j]][[k]] == 0) & (!is.na(MOSL[[row]][[j]][[k]][[1]])) & (length(MOSL[[row]][[j]][[k]])<=N.mosl-set.size(PLOS)) & !is.element(TRUE,is.element(unique(MOSL[[row]][[j]][[k]]),looplist.elements(PLOS)))) {
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


# sosl.prod returns the sign product of a set of spanning loops
# It takes:
# SOSL: a single set of spanning loop(s)

	sosl.prod <- function(C,SOSL) {

		loop.prod <- function(C,loop) {
			lprod <- 1
			if (length(loop) >1 ) {
				for (edge in 1:(length(loop))) {
					if (edge < length(loop)) {
						lprod <- lprod * C[loop[edge+1],loop[edge]]
						} else {
							lprod <- lprod * C[loop[1],loop[edge]]
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


# feedback returns the adjusted sum of products of SOSLs. It takes:
# C: the complimentary subsystem of a path through a CM

	feedback.C <- function(C) {
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


# path.prod() returns a sign product of a path. It takes:
# CM: a community matrix
# path: a path vactor through the community matrix
		path.prod <- function(CM,path) {
			pprod <- CM[path[2],path[1]]
			if (length(path) == 2) {
				if (path[1] == path[2]) {
					return(1)
					}
				return(pprod)
				}
			for (edge in 2:(length(path)-1) ) {
				pprod <- pprod * CM[path[edge+1],path[edge]]
				}	

			return(pprod)
		
			# end path.prod()
			}
	
		Sum <- 0
		for (path in enumerate.paths(CM,i,j)) {

			# if the complimentary subsystem C has zero elements, assign F.C -1
			if (length(CM[path.compliment(path,N)]) == 0) {
				F.C <- -1
			# otherwise, produce feedback F.C from the complimentary subsystem C
			 } else {
					F.C <- feedback.C(CM[path.compliment(path,N),path.compliment(path,N)])
					if (is.na(F.C)) {
						return(NA)
						}
					}

			P.F.C <- path.prod(CM,path)*F.C
			# if this is the first time through the loop
			if (Sum == 0) {
				Sum <- P.F.C
				} else {
						if (Sum == -1*P.F.C){
						return(NA)
						}
				}
			}

		return(Sum)
	
		# end sum.path.x.C()
		}

	N <- nrow(CM)
	namerows <- rownames(CM)
	F.N <- -1
	df.dc <- 1
	
	CEM <- matrix(c(NA),N,N,dimnames=list(namerows,namerows))

	if (!status) {
		for (i in 1:N) {
			for (j in 1:N) {
					Sum <- sum.path.x.C(CM,i,j,N)
					CEM[i,j] <- (df.dc*Sum)/F.N
				}
			}
		}

	if (status) {
		cat(" ",namerows,"\n")
		for (i in 1:N) {
			cat(namerows[i])
			for (j in 1:N) {
					Sum <- sum.path.x.C(CM,i,j,N)
					CEM[i,j] <- (df.dc*Sum)/F.N
					cat(" .")
				}
			cat("\n")
			}
		}
	
	if (out) {
		out.cm(CEM)
		}
	return(CEM)
	
	# end make.cem()
	}
