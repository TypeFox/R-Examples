#############################################################################
#These R functions were taken from the old package "kinship" by Atkinson and Therneauthat.
#This package is no longer in the CRAN repository:   http://cran.r-project.org/src/contrib/Archive/kinship/
# We have put it here instead of installing the archieve package to make it easier for general use.
########################################################################



# $Id: kinship.s,v 1.5 2003/01/04 19:07:53 Therneau Exp $
#
# Create the kinship matrix, using the algorithm of K Lange,
#  Mathematical and Statistical Methods for Genetic Analysis,
#  Springer, 1997, p 71-72.
#
# The rows (cols) of founders are just .5 * identity matrix, no further
#    processing is needed for them.
#  Parents must be processed before their children, and then a child's
#    kinship is just a sum of the kinship's for his/her parents.
#
kinship <- function(id, father.id, mother.id) {
	n <- length(id)
	if (any(duplicated(id))) stop("All id values must be unique")
	kmat <- diag(n+1) /2
	kmat[n+1,n+1]    <- 0  #if A and B both have "unknown" dad, this ensures
	# that they won't end up 'related' in the matrix
	
	pdepth <- kindepth(id, father.id, mother.id)
	# id number "n+1" is a placeholder for missing parents
	mrow <- match(mother.id, id, nomatch=n+1) #row number of the mother
	drow <- match(father.id, id, nomatch=n+1) #row number of the dad 
	
	# Those at depth==0 don't need to be processed
	# Subjects with depth=i must be processed before those at depth i+1.
	# Any parent is guarranteed to be at a lower depth than their children
	#  The inner loop on "i" can NOT be replaced with a vectorized expression:
	# sibs' effect on each other is cumulative.
	for (depth in 1:max(pdepth)) {
		indx <- (1:n)[pdepth==depth]
		for (i in indx) {
			mom <- mrow[i]
			dad <- drow[i]
			kmat[i,]  <- kmat[,i] <- (kmat[mom,] + kmat[dad,])/2
			kmat[i,i] <- (1+ kmat[mom,dad])/2
		}
	}
	
	kmat <- kmat[1:n,1:n]
	dimnames(kmat) <- list(id, id)
	kmat
}


#  $Id: kindepth.s,v 1.3 2003/04/17 21:57:12 Therneau Exp $
#
# Mark each person as to their depth in a pedigree
#   0 = founder
#   otherwise, depth = 1 + max(father's depth, mother's depth)
#
# Algorithm: founders =0
#    children of founders =1
#    children of children of founders = 2
#    children of children of children of founders = 3
#    ...
# Max depth = n-1
#
# If align=T, go one step further and try to make both parents of each
#   child have the same depth.  (This is not always possible).  It helps
#   the drawing program by lining up pedigrees that "join in the middle"
#   via a marriage.
#
kindepth <- function(id, dad.id, mom.id, align=F) {
	n <- length(id)
	if (n==1) return (0)  # special case of a single subject 
	midx <- match(mom.id, id, nomatch=0) # has a mother in the data set
	didx <- match(dad.id, id, nomatch=0) # has a father in the data set
	parents <- (midx==0 & didx==0)  #founders
	
	depth <- rep(0,n)
	# At each iteration below, all children of the current "parents" are
	#    labeled with depth 'i', and become the parents of the next iteration
	for (i in 1:n) {
		child  <- match(mom.id, id[parents], nomatch=0) +
				match(dad.id, id[parents], nomatch=0)
		
		if (all(child==0)) break
		if (i==n) 
			stop (paste("Impossible loop in the pedegree",
							"(someone would have to be born after their own child)"))
		
		parents <- (child>0) #next generation of parents
		depth[parents] <- i
	}
	if (!align) return(depth)
	
	#
	# Try aligning the pedigree.  The process will increase the depth of
	#   some branches, never decrease it.  
	# In the case of an inbred pedigree, there may not be a "perfect"
	#   alignment.  The algorithm below is:
	#      a. Find any mother-father pairs that are mismatched in depth.
	#         We think that aligning the top of a pedigree is more important
	#         than aligning at the bottom, so choose a mismatch pair of minimal
	#         depth.
	#      b. At least one member of the pair has depth = child-1, the other
	#         has depth < child-1.  Call these the good and bad sides.
	#      c. Chase up the good side, and get a list of all subjects connected
	#         to "good", including in-laws (spouse connections) & sibs that are
	#         at this level or above.  Call this agood (ancestors of good).
	#         We do not follow any connections at a depth lower than the 
	#         marriage in question, to get the highest marriages right.
	#         For the bad side, just get ancestors.
	#      d. Avoid pedigree loops!  If the chase list contains anyone in abad,
	#         then don't try to fix the alignment, otherwise:
	#         Push abad down, then run the pushdown algorithm to
	#             repair any descendents.
	# It may be possible to do better alignment when the pedigree has loops,
	#  but it is definitely beyond this program's abilities.  One tantalizing
	#  case appears in the framingham.s file: a pair of brothers married a
	#  pair of sisters.  Pulling one brother down fixes the other at the
	#  same time; the code below says "loop! stay away!".
	chaseup <- function(x, midx, didx) {
		new <- c(midx[x], didx[x])  # mother and father
		new <- new[new>0]
		while (length(new) >1) {
			x <- unique(c(x, new))
			new <- c(midx[new], didx[new])
			new <- new[new>0]
		}
		x
	}
	
	dads <- didx[midx>0 & didx>0]   # the father side of all spouse pairs
	moms <- midx[midx>0 & didx>0]
	# Get rid of duplicate pairs
	dups <- duplicated(dads + moms*n)
	if (any(dups)) {
		dads <- dads[!dups]
		moms <- moms[!dups]
	}
	npair<- length(dads)
	done <- rep(F, npair)  #couples that are taken care of
	while (T) {
		pairs.to.fix <- (1:npair)[(depth[dads] != depth[moms]) & !done]
		if (length(pairs.to.fix) ==0) break
		temp <- pmax(depth[dads], depth[moms])[pairs.to.fix]
		who <- min(pairs.to.fix[temp==min(temp)])  # the chosen couple
		
		good <- moms[who]; bad <- dads[who]
		if (depth[dads[who]] > depth[moms[who]]) {
			good <- dads[who]; bad <- moms[who]
		}
		abad  <- chaseup(bad,  midx, didx)
		if (length(abad) ==1 && sum(c(dads,moms)==bad)==1) {
			# simple case, a solitary marry-in
			depth[bad] <- depth[good]
		}
		else {
			agood <- chaseup(good, midx, didx)  #ancestors of the "good" side
			# For spouse chasing, I need to exclude the given pair
			tdad <- dads[-who]
			tmom <- moms[-who]
			while (1) {
				# spouses of any on agood list
				spouse <- c(tmom[!is.na(match(tdad, agood))],
						tdad[!is.na(match(tmom, agood))])
				temp <- unique(c(agood, spouse))
				temp <- unique(chaseup(temp, midx, didx)) #parents
				kids <- (!is.na(match(midx, temp)) | !is.na(match(didx, temp)))
				temp <- unique(c(temp, (1:n)[kids & depth <= depth[good]]))
				if (length(temp) == length(agood)) break
				else agood <- temp
			}
			
			if (all(match(abad, agood, nomatch=0) ==0)) {
				# shift it down
				depth[abad] <- depth[abad] + (depth[good] - depth[bad])
				#
				# Siblings may have had children: make sure all kids are
				#   below their parents.  It's easiest to run through the
				#   whole tree
				for (i in 0:n) {
					parents <- (depth==i)
					child <- match(mom.id, id[parents], nomatch=0) +
							match(dad.id, id[parents], nomatch=0)
					if (all(child==0)) break
					depth[child>0] <- pmax(i+1, depth[child>0])
				}
			}
		}
		done[who] <- T
	}
	if (all(depth>0)) stop("You found a bug in kindepth's alignment code!")
	depth
}

