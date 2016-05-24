trans <- function(rawdata, method.transform){
	if (method.transform == "squareroot")
		return(sqrt(rawdata))
	else if (method.transform == "log")
		return(log(rawdata))
	else if (method.transform == "PA")
		return((rawdata > 0) + 0) # +0 converts T/F to 1/0
	else if (typeof(method.transform) == "double")
		return(rawdata^method.transform)
	}

# computeAverage <- function(expectedprofile.simprof, num.expected){
#   ###
#   ### To Update
#   ###
# 	return(.Call("computeAverage", expectedprofile.simprof, num.expected))
# 	}

computeAverage <- function(expectedprofile.simprof, num.expected){
  temp.sums <- rep(0,length(expectedprofile.simprof[[1]][1,]))
  # Probably not feasible to replace this loop with a native R function
  # because of the way we need to loop over the list. But maybe we could.
  for (i in 1:num.expected){
    temp.sums <- temp.sums + expectedprofile.simprof[[i]][1,]
  }
  return(temp.sums/num.expected)
}
	
# tsComparison <- function(simulatedprofile.simprof, expectedprofile.average, num.simulated, teststatistic){
# 	pi_matrix <- matrix(data = NA, nrow = 1, ncol = num.simulated)
# 	piAsExtremeAsTS <- 0
# 	for (i in 1:num.simulated){
# 		pi_matrix[i] <- computeTestStatistic(simulatedprofile.simprof[[i]], expectedprofile.average)
#     #browser()
# 		if (pi_matrix[i] >= teststatistic)
# 			piAsExtremeAsTS <- piAsExtremeAsTS + 1
# 		}
# 	pval <- piAsExtremeAsTS/num.simulated
# 	return(pval)
# 	}

tsComparison <- function(simulatedprofile.simprof, expectedprofile.average, num.simulated, teststatistic){
  pi_vec <- sapply(simulatedprofile.simprof, FUN=computeTestStatistic, expectedprofile.average,
         simplify=TRUE)
  pval <- sum(pi_vec >= teststatistic)/num.simulated
  return(pval)
}

diveDeep <- function(rawdata, num.expected, num.simulated, method.cluster, method.distance, originaldata, 
						alpha, clust.order, startrow, pMatrix, side, const, silent, increment){
							
	side.startrow <- findNextStartRow(oldstartrow=startrow, mergematrix=clust.order, side=side)
	if(side.startrow > 0){
		if (toupper(side)=="LEFT")
			side.samples <- findleftgroup(startrow=startrow, mergematrix=clust.order)
		else if (toupper(side)=="RIGHT")
			side.samples <- findrightgroup(startrow=startrow, mergematrix=clust.order)
		side.rawdata <- matrix(data = NA, nrow = length(side.samples), ncol = ncol(originaldata))
		for (i in 1:ncol(side.rawdata))
			side.rawdata[,i] <- findsamples(rawdata=originaldata, samples=side.samples, column=i)		
		  simprof.results <- simprof.body(rawdata=side.rawdata, num.expected=num.expected, 
                                      num.simulated=num.simulated, 
			                                method.cluster=method.cluster, 
                                      method.distance=method.distance,
			                                originaldata=originaldata, alpha=alpha, 
                                      clust.order=clust.order, 
			                                startrow=side.startrow, pMatrix=pMatrix, 
                                      currentsamples=side.samples, const=const, 
                                      silent=silent, increment=increment)
    # there should be an easier way to pass all of these arguments to simprof.body...
    # also, this loop could maybe be broken up into two different apply()-style functions
		
		return(simprof.results)
		}
	else{
		simprof.results <- list()
		simprof.results[["samples"]] <- c(-1*side.startrow)
		simprof.results[["pval"]] <- pMatrix
		return(simprof.results)
		}
	}

findNextStartRow <- function(oldstartrow, mergematrix, side){
	if (toupper(side) == "LEFT")
		return (mergematrix[oldstartrow,1])
	else if (toupper (side) == "RIGHT")
		return (mergematrix[oldstartrow,2])
	}

findsamples <- function(rawdata, samples, column){
	temp <- matrix(data = NA, nrow=length(samples), ncol=1)
	for (i in 1:length(samples)) # this loop may be able to be replaced
		temp[i,1] <- rawdata[samples[i],column]
	return(temp)
	}
	
# computeTestStatistic <- function(rawdata.simprof, expectedprofile.average){
#   ###
#   ### To Update
#   ###
# 	return(.Call("computeTestStatistic", rawdata.simprof, expectedprofile.average, length(rawdata.simprof[1,])))
# 	}

computeTestStatistic <- function(rawdata.simprof, expectedprofile.average){
  return(sum(abs(rawdata.simprof[1,]-expectedprofile.average)))
}

	
findComplementaryIndices <- function(n, rawdata.thisgroup, rawdata.otherindices){
	### Generate a list of all possible sample IDs, NA those that are in the other group, remove the NAs, sample from that
	availableindices<-c(1:n)
	for(i in 1:length(rawdata.otherindices)) # we could maybe replace this loop with something else
		availableindices[rawdata.otherindices[i]]<-NA
	availableindices<-na.omit(availableindices)
	return(sample(availableindices, length(rawdata.thisgroup), replace=FALSE))
	}
	
genProfile <- function(rawdata, originaldata, num.expected, method.distance, 
                       const, silent, increment, type, undef.zero){
	### "expectedprofile" is a misnomer for these variables:
	### both expected and simulated profiles are generated with genProfile()
	expectedprofile.indices <- list()
	expectedprofile.data 	<- list()
	expectedprofile.simprof <- list()

	for (i in 1:num.expected){
		if (!(silent)){
			if (i %% increment == 0){
				print(paste(type,"iteration",i))
				}
			}
		# the interior of this loop looks like it could be replaced by apply()-style functions
    # but it might be very tricky
		expectedprofile.data[[i]] <- columnPermuter(mat=rawdata)
		expectedprofile.simprof[[i]] <- genSimilarityProfile(expectedprofile.data[[i]], method.distance, const, undef.zero)
		if (!dim(expectedprofile.simprof[[i]])[2]==1)
			expectedprofile.simprof[[i]] <- expectedprofile.simprof[[i]][,order(expectedprofile.simprof[[i]][2,])] ## order the matrix by ranks and keep the distances attached
		}
	return(expectedprofile.simprof)
	}
	
genSimilarityProfile <- function(rawdata.samples, method.distance, const, undef.zero) {
	if (is.function(method.distance))
		rawdata.dist <- method.distance(rawdata.samples)
	else if (method.distance == "braycurtis")
		rawdata.dist <- braycurtis(rawdata.samples, const, undef.zero)
	else if (method.distance == "czekanowski")
	  rawdata.dist <- czekanowski(rawdata.samples, const, undef.zero)
	else if (method.distance == "actual-braycurtis")
	  rawdata.dist <- braycurtis(preBCstandardization(rawdata.samples), const, undef.zero)
	else
		rawdata.dist <- dist(rawdata.samples, method.distance)
	rawdata.distvec <- as.vector(rawdata.dist)
	
  # Ties (equal values) ARE possible in rawdata.distvec and the rank() call used the default behavior.
  # The default behavior of rank() is ties.method = "average", which can result in non-integer ranks
  # However, the original, loop-based implementation of computeAverage in C ignored these,
  # essentially making ties.method="first". Because we are re-implemented the computeAverage function,
  # we'll also change the rank() call to make this behavior explicit. If "first" were NOT used and 
  # the computeAverage() function were re-implemented so that it averaged all of the ranks matching,
  # say, 1:max, it would miss the decimal ranks and create more problems. 
  rawdata.ranks 	<- rank(rawdata.distvec, ties.method="first")
  
	rawdata.simprof <- rbind(rawdata.distvec, rawdata.ranks)
	return(rawdata.simprof)
	}
	
columnPermuter <- function(mat){
  numrows <- nrow(mat); numcols <- ncol(mat)
  newmat <- matrix(NA,nrow=numrows,ncol=numcols)
  # this loop could maybe be optimized...
  for (current.column in 1:numcols){
    newmat[,current.column] <- mat[sample(x=1:numrows, size=numrows, replace=FALSE),current.column]
    }
  return(newmat)
  ###
  ### To Update
  ###
	#return(.Call("columnPermuter", mat, nrow(mat), ncol(mat)))
	}
	
### findleftgroup, findrightgroup, and treedive are used for making sense of hclust's $merge results	
findleftgroup <- function(startrow, mergematrix){
	if (mergematrix[startrow,1] < 0) # negative means we've encountered a singleton
		samples <- mergematrix[startrow,1]
	else if (mergematrix[startrow,1] > 0) # positive means we've encountered a subgroup
		samples <- treedive(mergematrix[startrow,1], mergematrix)
	return(-1*samples)
	}
			
findrightgroup <- function(startrow, mergematrix){
	if (mergematrix[startrow,2] < 0) # negative means we've encountered a singleton
		samples <- mergematrix[startrow,1]
	else if (mergematrix[startrow,2] > 0) # positive means we've encountered a subgroup
		samples <- treedive(mergematrix[startrow,2], mergematrix)
	return(-1*samples)
	}			
			
treedive <- function(startrow, mergematrix){
	
	## left
	if (mergematrix[startrow,1] < 0)
		samples.left <- mergematrix[startrow,1]
	else 
		samples.left <- treedive(mergematrix[startrow,1], mergematrix)
		
	## right
	if (mergematrix[startrow,2] < 0)
		samples.right <- mergematrix[startrow,2]
	else
		samples.right <- treedive(mergematrix[startrow,2], mergematrix)
		
	return(append(samples.left, samples.right))
	}

# braycurtis <- function(data, const){
#   ###
#   ### To Update
#   ###
# 	return(as.dist(.Call("braycurtis", data, nrow(data), ncol(data), const)))
# 	}

# braycurtis2 <- function(data, const){
#   tmp.bc <- matrix(NA,ncol(data),ncol(data))
#   for (i in 1:(ncol(data)-1)){
#     for (j in (i+1):ncol(data)){
#       num <- 2*sum(pmin(data[,i],data[,j]))
#       denom <- sum(data[,i],data[,j])
#       denom <- denom + 2*const # zero-adjustment
#       if (denom != 0)
#         tmp.bc[j,i] <- 100*(1-num/denom) # convert from similarity to dissimilarity and to % scale
#       else
#         tmp.bc[j,i] <- 0 # double check this
#     }
#   }
#   return(as.dist(tmp.bc))
# }


## maybe rewrite using combn() and apply()?
czekanowski <- braycurtis <- function(data, const, undef.zero=TRUE){
  data <- t(data) # to retain backward compatibility with old code
  tmp.bc <- matrix(NA,ncol(data),ncol(data))
  for (i in 1:(ncol(data)-1)){
    for (j in (i+1):ncol(data)){
      num <- sum(abs(data[,i]-data[,j]))
      denom <- sum(data[,i],data[,j])
      denom <- denom + 2*const # zero-adjustment
      #browser()
      if (denom != 0){
        tmp.bc[j,i] <- 100*num/denom # convert from similarity to dissimilarity and to % scale
      }
      else if (undef.zero==TRUE){
        tmp.bc[j,i] <- 0 # double check this
        warning("Bray-Curtis resulted in undefined values (replaced by 0).",call.=FALSE)
      }
      else{
        tmp.bc[j,i] <- 100*num/denom
        warning("Bray-Curtis resulted in undefined values.",call.=FALSE)
      }
    }
  }
  return(as.dist(tmp.bc))
}

## much too slow (for BC1,2,3)
# user  system elapsed 
# 12.53    0.04   12.70 
# user  system elapsed 
# 32.12    0.22   32.58 
# user  system elapsed 
# 60.31    0.26   61.19 
#
# braycurtis3 <- function(data, const, undef.zero=TRUE){
#   data <- t(data) # to retain backward compatibility with old code
#   tmp.bc <- matrix(NA,ncol(data),ncol(data))
#   for (i in 1:(ncol(data)-1)){
#     for (j in (i+1):ncol(data)){
#       num <- 2*sum(pmin(data[,i],data[,j]))
#       num <- num + 2*const
#       denom <- sum(data[,i],data[,j])
#       denom <- denom + 2*const # zero-adjustment
#       #browser()
#       if (denom != 0){
#         tmp.bc[j,i] <- 100*(1-num/denom) # convert from similarity to dissimilarity and to % scale
#       }
#       else if (undef.zero==TRUE){
#         tmp.bc[j,i] <- 0 # double check this
#         warning("Bray-Curtis resulted in undefined values (replaced by 0).",call.=FALSE)
#       }
#       else{
#         tmp.bc[j,i] <- 100*(1-num/denom)
#         warning("Bray-Curtis resulted in undefined values.",call.=FALSE)
#       }
#     }
#   }
#   return(as.dist(tmp.bc))
# }

preBCstandardization <- function(data){
  sp.max <- apply(data, 2, max) # maximum species value across sampling sites
  prop.max.across <- sweep(data,2,sp.max,"/") # proportion of max value attained
  return(prop.max.across/apply(prop.max.across,1,sum)) # return relative scores within each sampling site
}
	
### this function will put the p-value in a column on the same row as the split it is testing.
pTracker <- function(PMmat, startrow, pval){
	PMmat[startrow,3] <- pval
	return(PMmat)
	}
	
pTrackerMerge <- function(PMmatL, PMmatR){
	for (i in 1:nrow(PMmatL))
		if (is.na(PMmatL[i,3]) && !is.na(PMmatR[i,3]))
			PMmatL[i,3] <- PMmatR[i,3]
	return(PMmatL)
	}
