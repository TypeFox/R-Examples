allExplained = function(genotype,cspAlleles,knownWithStutter,alleleNames,doDoubleStutter=FALSE,doOverStutter=FALSE)
	{
	genotypeAlleles = as.numeric(alleleNames[genotype])
	genWithStutter = c(genotypeAlleles,genotypeAlleles-1)
        if(doDoubleStutter) genWithStutter = c(genWithStutter,genotypeAlleles-2)
	if(doOverStutter) genWithStutter = c(genWithStutter,genotypeAlleles+1)
	all(cspAlleles%in%c(genWithStutter,knownWithStutter))
	}

explain.all.peaks = function(cspPresence,profPresence,knownProfs,alleleNames,nUnknowns,cspAlleles,cspHeights, doDoubleStutter=FALSE,doOverStutter=FALSE,doDropin=FALSE)
	{
  	# Catch early. Shouldn't be a problem here, but will be at some point. 
	if(!is.matrix(cspPresence)) stop("Expected a matrix as input.")
	if(!is.matrix(profPresence)) stop("Expected a matrix as input.")
	# add extra unknown contributor if dropin and no unknowns
	#if(nUnknowns == 0 && doDropin == TRUE) nUnknowns = 1
	if(ncol(profPresence)!=0)
	    {
	    # get known alleles
    	knownIndex = sapply(unlist(knownProfs),FUN=function(x) which(alleleNames==x))
    	knownAlleles = as.numeric(alleleNames[knownIndex])
	# include stuttered known alleles
	knownWithStutter = c(knownAlleles,knownAlleles-1)
	if(doDoubleStutter) knownWithStutter = c(knownWithStutter, knownAlleles-2)
	if(doOverStutter) knownWithStutter = c(knownWithStutter, knownAlleles+1)
	    } else {
        knownIndex = c()
        knownWithStutter = c()
	    }
	# get csp Alleles
	cspAlleles = alleleNames[row(cspPresence)[which(cspPresence)]]
	# if no unknowns return knowns
	if(nUnknowns==0) 
		{
		check = !all(round(as.numeric(cspAlleles),1)%in%round(as.numeric(knownWithStutter),1))
		if(doDropin) check=FALSE
		if(check) 
		    {
		    stop(paste0("Not enough contributors to explain CSP at locus ", colnames(knownProfs)))
		    } else {
		    return(matrix(knownIndex,ncol=1))
		    }
		}
	# get all genotype combinations for unknowns
	genCombs = all.genotypes.per.locus(length(alleleNames),nUnknowns)
	# if dropin add knowns and return
	if(doDropin)
		{
		# add known profiles as last contributors
		if(length(knownIndex)!=0) 
			{
			genCombs = rbind(genCombs,matrix(rep(knownIndex,times=ncol(genCombs)),
					ncol=ncol(genCombs)))
			}
		return(genCombs)
		}
	# find which combinations explain all peaks
	index = apply(genCombs,MARGIN=2,FUN=function(x) allExplained(x,cspAlleles,knownWithStutter,alleleNames,doDoubleStutter,doOverStutter))
	if(length(which(index))==0) stop(paste0("Not enough contributors to explain CSP at locus ", colnames(knownProfs)))
	genCombs = genCombs[,index,drop=FALSE]
	# add known profiles as last contributors
	if(length(knownIndex)!=0) 
		{
		genCombs = rbind(genCombs,matrix(rep(knownIndex,times=ncol(genCombs)),ncol=ncol(genCombs)))
		}
	return(genCombs)
	}

# Defines prod.matrix from R.
prod.matrix.row <- function(x) {
  # Fast column-wise product.
  #
  # Sadly, this is faster than apply(x, 1, prod)
  y=x[,1]
  for(i in 2:ncol(x))
  y=y*x[, i]
  return(y)
}

ethnic.database.lus <- function(ethnic, loci=NULL, afreq=NULL) {
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
    result <- matrix(c(locus[[ethnic]], locus[["BP"]], locus[["LUS"]]), , 3)
    result[is.na(result[, 2]), 2] <- 0
    rownames(result) <- locus$Allele
    if(any(is.na(result[,3])))
	{
    	result = fill.unknown.LUS(result)
	}
    return(result[result[, 1] > 0, ])
    return(result)
  }
  # now apply function over all loci.
  result <- sapply(loci, filter.locus)
  # Then center around mean fragment length
  s1 = sum( sapply(result, function(n) sum(n[, 1] * n[, 2], na.rm=TRUE)) )
  s2 = sum( sapply(result, function(n) sum(n[, 1], na.rm=TRUE)) )
  for(j in 1:length(result)) result[[j]][, 2] <- result[[j]][, 2] - s1/s2
  return(result)
}
