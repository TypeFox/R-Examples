ethnic.database <- function(ethnic, loci=NULL, afreq=NULL) {
  # Reformats allele database to include only given ethnic group and loci.
  #
  # Filters database down to a single ethnic group and the loci of interests.
  # Removes alleles which are not present in given ethnic group. Centers
  # average length across filtered database.
  # 
  # Parameters:
  #   ethnic: Name of ethnic group. Should correspond to what's in afreq.
  #   loci: Names of the loci to use. Or use all.
  #   afreq: Frequency table. If NULL, loads in frequency table provided with
  #          likeLTD package.
 
  # Load frequency database if needed. 
  if(is.null(afreq)) afreq <- load.allele.database()
  # figueres out loci
  if(is.null(loci)) loci <- t(unique(afreq['Marker']))

  # Function which recreates the input for a single locus.
  # Input of a locus consists of one row per allele type.
  # Also filters out alleles with 0 frequencies.
  filter.locus <- function(n) {
    locus <- afreq[afreq$Marker == n, ]
    result <- matrix(c(locus[[ethnic]], locus[["BP"]]), , 2)
    result[is.na(result[, 2]), 2] <- 0
    rownames(result) <- locus$Allele
    return(result[result[, 1] > 0, ])
  }
  # now apply function over all loci.
  result <- sapply(loci, filter.locus)
  # Then center around mean fragment length
  s1 = sum( sapply(result, function(n) sum(n[, 1] * n[, 2], na.rm=TRUE)) )
  s2 = sum( sapply(result, function(n) sum(n[, 1], na.rm=TRUE)) )
  for(j in 1:length(result)) result[[j]][, 2] <- result[[j]][, 2] - s1/s2
  return(result)
}

adjust.frequencies <- function(alleleDb, queriedAlleles, adj=1, fst=0.02) {
  # Adjust frequencies for current case.
  #
  # Parameters:
  #   alleleDb: Allele database for given locus, for a given ethnic group (see
  #             ethnic.database), completed for rare alleles in queried profile
  #             (see missing.alleles). The loci in both alleleDb and
  #             queriedAlleles should occur in the same order.
  #   queriedAlleles: Profile of the queried individual, with colums as locus.
  #                   This will generally be something returned by
  #                   known.alleles.
  #   adj: sampling adjustment applied to the alleles of Q to avoid very low
  #        counts for rare alleles, and to allow for the allele counts to take
  #        into account the genotype of Q.  The default value is adj = 1; there
  #        are reasons to prefer adj = 2 but the value of adj is much less
  #        important than Fst (see below) unless the database size is very
  #        small.  
  #   fst: allows for shared ancestry of Q with X. Recommended that Fst should
  #        be at least 0.02, and may need to be as high as 0.05 in some
  #        populations (e.g.  small, isolated subpopulations of the population
  #        from which the reference database has been drawn).
  adjust.per.locus <- function(alleleDbLocus, queriedLocus, adj=1, fst=0.02) {
    # Applies operations to a single locus.
    homozygote <- as.integer(queriedLocus[1] == queriedLocus[2])
    alleleDbLocus[queriedLocus, 1] <- 
      alleleDbLocus[queriedLocus, 1] + adj * (1 + homozygote)
    #  shared ancestry adjustments.
    alleleDbLocus[, 1] <-
      alleleDbLocus[, 1] / sum(alleleDbLocus[, 1]) * (1 - fst) / (1 + fst) 
    alleleDbLocus[queriedLocus, 1] <-
      alleleDbLocus[queriedLocus, 1] + fst / (1 + fst) * (1 + homozygote)
    return(alleleDbLocus)
  }
  return(mapply(adjust.per.locus, alleleDb, queriedAlleles, adj=adj,
                fst=fst))
}

# R version of computing all genotypes per locus.
all.genotypes.per.locus.R = function(nAlleles, nContrib=1) {
  # All compatible agreggate profiles for a given locus.
  #
  # Computes all compatible combined profiles from n contributors.
  # This is a permutation on combinations: it grows way too fast as
  # :math:`N=[0.5*nAlleles*(nAlleles+1)]^{nContrib}`.
  # 
  # Parameters:
  #   nAlleles: number of alleles. 
  #   nContrib: number of unprofiled contributors.
  if(nContrib == 0 || nAlleles == 0) return(matrix(nrow=0, ncol=0))
  # All compatible genotypes for a single contributor.
  singleContributor = t(combinations(nAlleles, 2, repeats.allowed=T))
  if(nContrib == 1) return(singleContributor)
  # All compatible permutations of "nContrib" contributors.
  nContribPerms = permutations(ncol(singleContributor), nContrib,
                               repeats.allowed=T)
  # The next two lines create a matrix where the first two columns is
  # singleContributor in the order given by the first column of nContribPerms.
  # The next two columns is singleContributor in the order given by
  # nContribPerms's second column. And so on and so forth.
  apply.perms <- function(n, sing, perms) sing[, perms[, n]]
  rows <- lapply(1:nContrib, apply.perms, sing=singleContributor,
                 perms=nContribPerms) 
  do.call(rbind, rows)
}
# Computes all possible genetic make-ups for the set of
# unprofiledcontributors. 
all.genotypes.per.locus <- function(nAlleles, nContrib=1) {
  if(nAlleles < 0) stop("Negative or null number of alleles.")
  if(nContrib < 0) stop("Negative number of contributors.")
  if(nContrib == 0 || nAlleles == 0) return(matrix(nrow=0, ncol=0))
  nContrib = as.integer(nContrib)
  a = t(combinations(nAlleles, 2, repeats.allowed=TRUE))
  if(is.null(nContrib) || is.null(a)) stop("something went very wrong.")
  .Call(.cpp.allGenotypesPerLocus, nContrib, a, PACKAGE="likeLTD")
}

# Profiles of unknown contributors for given locus.
compatible.genotypes = function(cspPresence, profPresence, alleleNames,
                                nUnknowns, dropin=FALSE, missingReps=NULL) {
  
  # Catch early. Shouldn't be a problem here, but will be at some point. 
  if(!is.matrix(cspPresence)) stop("Expected a matrix as input.")
  if(!is.matrix(profPresence)) stop("Expected a matrix as input.")

  # Case where there are no unknown contributors but dropin is modelled.
  # Seems it's formally equivalent to one unknown and dropin.
  if(nUnknowns == 0 && dropin == TRUE) nUnknowns = 1

  # Compute all genotype permutations. 
  genotypes = all.genotypes.per.locus(length(alleleNames), nUnknowns)
  if(dropin) return(genotypes) # Dropin case: can return immediately.

  # Case without dropin: check that there are enough unknown contributors to
  # account for the alleles in the CSP that are not in the known profile.
  # Might be able to return early here also.
  # Reduce matrices to a single logical vector
  if(is.null(missingReps)) missingReps = rep(TRUE, nrow(cspPresence))
  cspPresence  = rowSums(cspPresence[, c(!missingReps), drop=FALSE]) > 0
  profPresence = rowSums(profPresence) > 0

  # required: indices of the alleles in the crime scene which are not in the
  #           known profiles. Indices are for freqLocus rows. In the case of
  #           dropins, then there are no required alleles.
  required = which(cspPresence & !profPresence) 

  # Not enough contributors, return empty matrix.
  if(length(required) > 2*nUnknowns)  {
    stop(sprintf("Not enough unknown contributors:
       %d contributors, 
       %d required allele (%s).", 
       nUnknowns, length(required),
       paste(alleleNames[required], collapse=", ")))
  }
  if(nUnknowns == 0) return(matrix(0, nrow=1, ncol=1))

  hasRequired <- apply(genotypes, 2, function(n) all(required %in% n))
  genotypes[, hasRequired, drop=FALSE]
}

known.epg.per.locus <- function(rcont, degradation, fragmentLengths,
                                dropoutPresence) {
  # Creates "electropherogram" vector for a given locus.
  #
  # In practice, the "electropherogram" is a vector giving the dose for each
  # allele present in the known profiles.
  #
  # Parameters:
  #   rcont: relative contributions from each individual. It should be vector
  #          of the same length as the number of individuals in the profile.
  #   degradation: degradation of the DNA from each individual. 
  #   fragmentLengths: The str lengths for a given locus.
  #   dropoutPresence: Presence matrix for known profiles with dropout.
  # Returns: A vector with an element per allele. Each element is zero if is
  #          not within the known profile, or it is its dose within the
  #          candidate CSP.

  if(ncol(dropoutPresence) == 0) return(array(0.0, nrow(dropoutPresence)))
  indices = which(dropoutPresence > 0) - 1
  indivs = trunc(indices / nrow(dropoutPresence)) + 1
  alleles = indices %% nrow(dropoutPresence) + 1
  doses = mapply( function(d, r, L, het) het * r * (1.0 + d)^{-L}, 
                  degradation[indivs], 
                  rcont[indivs], 
                  fragmentLengths[alleles],
                  rep(3 - colSums(dropoutPresence), colSums(dropoutPresence)) )
  result = array(0.0, nrow(dropoutPresence))
  for(i in 1:length(alleles)) 
    result[alleles[i]] = result[alleles[i]] + doses[i]
  return(result)
}

all.epg.per.locus <- function(rcont, degradation, dropoutPresence,
                              knownFragLengths, genotypes, addProfiles) {
  # Creates "electropherogram" for each compatible set of unknown contributors.
  #
  # Parameters:
  #   rcont: relative contributions from each individual. It should be vector
  #          of the same length as the number of individuals in the profile.
  #   degradation: degradation of the DNA from each individual. 
  #   dropoutPresence: Presence matrix for known profiles with dropout.
  #   knownFragLengths: matrix with fragment lengths for each allele in each
  #                     known profile.
  #   fragLengths: matrix with fragment lengths for each allele each computed
  #                profile.
  #   genotypes: matrix with all compatible genotypes. 
  # Returns: A vector with an element per allele. Each element is zero if is
  #          not within the known profile, or it is its dose within the
  #          candidate CSP.

  # "electropherogram" from known profiles
  knownEPG <- known.epg.per.locus(rcont, degradation, knownFragLengths,
                                  dropoutPresence)

  # allEPG: for each set of unknown contributors, total EPG of known and
  #         unknowns.
  # At point of creation, contains only component from known profiles.
  allEPG = matrix(knownEPG, ncol=ncol(genotypes), nrow=length(knownEPG)) 
  # Add into allEPG the components from unknown contributors.
  if(addProfiles) {
    # construct matrices with unknown doses
    nUnknowns = nrow(genotypes)
    indices = ncol(dropoutPresence) + rep(1:(nUnknowns/2), rep(2, nUnknowns/2))
    unknownDoses = matrix(ncol=length(indices), nrow=length(knownFragLengths))
    for(i in 1:ncol(unknownDoses))
      unknownDoses[, i] = rcont[indices[i]] * 
                          (1.0 + degradation[indices[i]])^-knownFragLengths

    # The following loop is done in C. 
    #   for(j in 1:ncol(genotypes)) {
    #     for(u in 1:nUnknowns) {
    #       index = genotypes[u, j]
    #       allEPG[index, j] = allEPG[index, j]  + unknownDoses[index, u]
    #     }
    #   }
    .Call(.cpp.addProfilesToEPG, allEPG, genotypes, unknownDoses, PACKAGE="likeLTD")
  }
  return(allEPG)
}

# Defines prod.matrix from R.
prod.matrix.col <- function(x) {
  # Fast column-wise product.
  #
  # Sadly, this is faster than apply(x, 1, prod)
  y=x[1, ]
  for(i in 2:nrow(x))
  y=y*x[i, ]
  return(y)
}

selective.col.prod <- function(condition, input) {
  # Row-wise product over selected columns or elements.
  #
  # If all(condition == FALSE), returns 1. Otherwise does one of the following.
  #
  # - If condition and input matrices: 
  #   Performs column-wise multiplication of input elements, ignoring those
  #   which are *not* selected by the condition matrix.
  # - If condition is a matrix and input a vector: 
  #   Equivalent to creating a matrix whith the same dimensions as condition
  #   from the input, where row ``i`` is ``input[condition[i, ]]``, then apply
  #   product.
  # - If condition is a vector and input is a matrix: 
  #   Performs column-wise multiplication of selected input columns only.
  # - If condition and input are vectors:
  #   Sets elements in input which are FALSE to 1, and returns it.
  # 
  # Row-wise multiplication means column 1 times column 2 times... e.g. along
  # rows. 
  #
  # Parameters:
  #   condition: a logical vector or a logical matrix. See above.
  #   input: input matrix with data for which to perform calculation.
  # Return: 
  #   - if condition is always FALSE, returns the scalar 1.0
  #   - otherwise a vector with the same length as there are rows in input
  # 
  # Note: Really should be done through R's method thingie...
  #       Maybe once I've read up on it.

  if(!any(condition)) return(1)

  # condition is a matrix
  if(is.matrix(condition)) {  
    # but input is a vector. Then build a matrix from condition and input.
    if(!is.matrix(input)) {
      if(nrow(condition) != length(input)) 
        stop("condition does not have as many columns as input has elements.") 
      # Builds input matrix
      input = matrix(input, nrow=length(input), ncol=ncol(condition))
    } else if(any(dim(condition) != dim(input)))
      stop("condition and input have different dimensions.") 

    input[!condition] = 1
    return(prod.matrix.col(input))
  }

  # condition is a vector and input is a matrix.
  if(is.matrix(input)) {
    if(nrow(input) != length(condition))
      stop("condition does not have as many elements as input has columns.") 

    nbAlleles = sum(condition)
    if(nbAlleles != 1) {
      input = input[condition, , drop=FALSE] 
      input[is.na(input)] = 1 
      return(prod.matrix.col(input))
    }

    output = input[condition, ]
    output[is.na(output)] = 1
    return(output)
  }

  # condition and input are vectors
  if(length(input) != length(condition))
    stop("condition does not have as many elements as input.") 
  output = input
  output[!condition] = 1.0
  output[is.na(output)] = 1.0
  return(output)
}



# heterozygote matrix needed for linked locus match probability
hetLinked = function(pA,pB,fst)
	{
	C = 2*(1-fst)*pA*pB
	AB1 = C *  (((2*fst)+(1-fst)*(pA+pB))/(2*(1+fst)))
	AB0 = C * het(pA, pB, fst=fst)
	return(c(AB0,AB1,C))
	}

# homozygote matrix needed for linked locus match probability
homLinked = function(pA,fst)
	{
	D = pA*(fst+(1-fst)*pA)
	AB1 = D * ((2*fst)+(1-fst)*pA)/(1+fst)
	AB0 = AB1 * ((3*fst)+(1-fst)*pA)/(1+(2*fst))
	return(c(AB0,AB1,D))
	}

# Z matrix needed for linked locus match probability siblings
Zsib = function(R)
	{
	X = R^2+(1-R)^2
	Y = 2*R*(1-R)
	out = matrix(rev(c(
			(X^2)/4,
			(X*Y)/2,
			(Y^2)/4,
			(X*Y)/2,
			(X^2+Y^2)/2,
			(X*Y)/2,
			(Y^2)/4,
			(X*Y)/2,
			(X^2)/4
			)),ncol=3,nrow=3)
	return(t(out))
	}

# Z matrix needed for linked locus match probability uncle/nephew
Zuncle = function(R)
	{
	W = (0.5*(R^2+(1-R)^2)*(1-R))+0.25*R
	out = matrix(rev(c(
			0,
			0,
			0,
			0,
			W,
			0.5-W,
			0,
			0.5-W,
			W
			)),ncol=3,nrow=3)
	return(t(out))
	}

# Z matrix needed for linked locus match probability half uncle/nephew
Zhuncle = function(R)
	{
	X = R^2+(1-R)^2
	nR = 1-R
	out = matrix(rev(c(
			0,
			0,
			0,
			0,
			0.25*X*nR,
			(1-(X*nR))*0.25,
			0,
			(1-(X*nR))*0.25,
			(2+(X*nR))*0.25
			)),ncol=3,nrow=3)
	return(t(out))
	}

# Z matrix needed for linked locus match probability cousin
Zcousin = function(R)
	{
	V = (0.25*(R^2+(1-R)^2)*(1-R)^2)+(R^2)/8
	out = matrix(rev(c(
			0,
			0,
			0,
			0,
			V,
			0.25-V,
			0,
			0.25-V,
			0.5+V
			)),ncol=3,nrow=3)
	return(t(out))
	}

# Z matrix needed for linked locus match probability grandparent
Zgrandparent = function(R)
	{
	nR = 1-R
	out = matrix(rev(c(
			0,
			0,
			0,
			0,
			nR/2,
			R/2,
			0,
			R/2,
			nR/2
			)),ncol=3,nrow=3)
	return(t(out))
	}

# Z matrix needed for linked locus match probability half siblings
Zhsibs = function(R)
	{
	X = R^2+(1-R)^2
	Y = 2*R*(1-R)
	out = matrix(rev(c(
			0,
			0,
			0,
			0,
			X/2,
			Y/2,
			0,
			Y/2,
			X/2
			)),ncol=3,nrow=3)
	return(t(out))
	}

# function to find the match probability including linked loci
linkedMatchProb = function(hypothesis,linkedIndex,R)
    {
	rr = hypothesis$relatedness
	fst = hypothesis$fst
	relationship = hypothesis$relationship
	ideal.match <- c()
	for(j in 1:ncol(hypothesis$queriedProfile))
		{
		
		if(j%in%c(linkedIndex)) next
		af = hypothesis$alleleDb[j][[1]]
		kn = hypothesis$queriedProfile[,j][[1]]
		p1 = af[row.names(af)==kn[1],1]
		p2 = af[row.names(af)==kn[2],1]

        	# undo most of fst adjustment
	        p1 = (p1*(1+fst))-fst
        	p2 = (p2*(1+fst))-fst


        	if(kn[1]==kn[2])
        	    {
        	    # undo last bit of fst adjustment
        	    p1 = (p1-fst)/(1-fst)
        	    p2 = (p2-fst)/(1-fst)
        	    # get match new fst adjusted match prob
        	    ideal.match = c(ideal.match,rr[2]+(rr[1]*((2*fst+(1-fst)*p1)/(1+fst)))+((1-sum(rr))*hom(p1,fst=fst)))
        	    } else {
        	    # undo last bit of fst adjustment
        	    p1 = p1/(1-fst)
        	    p2 = p2/(1-fst)
        	    # get match new fst adjusted match prob
        	    ideal.match = c(ideal.match,rr[2]+(rr[1]*((fst+(1-fst)*((p1+p2)/2))/(1+fst)))+((1-sum(rr))*het(p1,p2,fst=fst)))
        	    }
		}
	# joint probabilities for linked loci
	for(j in 1:nrow(linkedIndex))
		{
		# allele probabilities for first locus
		af1 = hypothesis$alleleDb[linkedIndex[j,1]][[1]]
		kn1 = hypothesis$queriedProfile[,linkedIndex[j,1]][[1]]
		p1 = af1[row.names(af1)==kn1[1],1]
		p2 = af1[row.names(af1)==kn1[2],1]
		# allele probabilities for second locus	
		af2 = hypothesis$alleleDb[linkedIndex[j,2]][[1]]
		kn2 = hypothesis$queriedProfile[,linkedIndex[j,2]][[1]]
		p3 = af2[row.names(af2)==kn2[1],1]
		p4 = af2[row.names(af2)==kn2[2],1]

		# undo most of fst adjustment
	        p1 = (p1*(1+fst))-fst
        	p2 = (p2*(1+fst))-fst	
	        p3 = (p3*(1+fst))-fst
        	p4 = (p4*(1+fst))-fst

		# Z matrix
		if(relationship==2)
			{
			Ztmp = Zsib(R[j])
			} else if(relationship==3) {
			Ztmp = Zuncle(R[j])
			} else if(relationship==4) {
			Ztmp = Zhuncle(R[j])
			} else if(relationship==5) {
			Ztmp = Zcousin(R[j])
			} else if(relationship==6) {
			Ztmp = Zgrandparent(R[j])
			} else if(relationship==7) {
			Ztmp = Zhsibs(R[j])
			}

        	if(kn1[1]==kn1[2])
			{
        	    	# undo last bit of fst adjustment
        	    	p1 = (p1-fst)/(1-fst)
        	    	p2 = (p2-fst)/(1-fst)
			# A matrix
			A = matrix(homLinked(p1,fst=fst),ncol=3)
			# genotype probability
			factor1 = A[,3]
			} else {
        	    	# undo last bit of fst adjustment
        	    	p1 = p1/(1-fst)
        	    	p2 = p2/(1-fst)
			# A matrix
			A = matrix(hetLinked(p1,p2,fst=fst),ncol=3)
			# genotype probability
			factor1 = A[,3]
			}
		if(kn2[1]==kn2[2])
			{
			# undo last bit of fst adjustment
        	    	p3 = (p3-fst)/(1-fst)
        	    	p4 = (p4-fst)/(1-fst)
			# B matrix
			B = matrix(homLinked(p3,fst=fst),nrow=3)
			# genotype probability
			factor2 = B[3,]
			} else {
			# undo last bit of fst adjustment
        	    	p3 = p3/(1-fst)
        	    	p4 = p4/(1-fst)
			# B matrix
			B = matrix(hetLinked(p3,p4,fst=fst),nrow=3)
			# genotype probability
			factor2 = B[3,]
			}

		out = (A%*%Ztmp%*%B)/prod(c(factor1,factor2))
		ideal.match = c(ideal.match,out)
		}
	ideal.match = 1/prod(ideal.match)
    return(ideal.match)
    }


# Function to find the difference between the match probability with/without linkage
linkage = function(hypothesis,factor=TRUE)
	{
	# combinations of CSP markers
	combs = combinations(length(hypothesis$alleleDb),2) 
	# recombination fractions for each combination
	R = apply(combs,MARGIN=1,FUN=function(x) hypothesis$linkageInfo[which(toupper(rownames(hypothesis$linkageInfo))==toupper(names(hypothesis$alleleDb))[x[1]]),which(toupper(colnames(hypothesis$linkageInfo))==toupper(names(hypothesis$alleleDb))[x[2]])])
	# find which recombinations actually happen
	index = suppressWarnings(sapply(R,FUN=function(x) !is.na(unlist(x))&length(unlist(x)!=0)))
	# replace zero length with FALSE or NA
	zeroIndex = which(sapply(index,FUN=function(x) length(x)==0))
	for(i in 1:length(zeroIndex)) 
	    {
	    index[[zeroIndex[i]]] = FALSE
	    R[[zeroIndex[i]]] = NA
	    }
	# matrix of linked loci
	toLink = combs[unlist(index),,drop=FALSE]
	if(any(table(toLink)>1)) stop("A locus is linked to more than one other locus")
	# do not perform linkage if an unknown relationship
	if(is.null(hypothesis$relationship)) toLink = NULL
	R = unlist(R)[!is.na(unlist(R))]
	# get non-linked match probability
	nonLinked = matchProb(hypothesis,hypothesis$relatedness,hypothesis$fst)
    	# get linked match probability
	if(nrow(toLink)==0)
	    {
        linked = nonLinked
	    } else {
        linked = linkedMatchProb(hypothesis,toLink,R)
	    }
	if(!factor) return(linked)
	# return 
	out = linked/nonLinked
    return(out)
	}

# wrapper to get linked or unlinked match probability
getMatchProb = function(hypothesis)
	{
	if(hypothesis$relationship%in%c(0,1))
		{
		matchProb(hypothesis,hypothesis$relatedness,hypothesis$fst)
		} else {
		linkage(hypothesis,FALSE)
		}
	}

