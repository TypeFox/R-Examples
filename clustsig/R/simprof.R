### This function is going to act as a wrapper for the meat of the math
### Primarily to take care of things that happen only once like hclust()
simprof <- function(data, num.expected=1000, num.simulated=999, 
					          method.cluster="average", method.distance="euclidean", 
                    method.transform="identity", alpha=0.05, 
                    sample.orientation="row", const=0, 
                    silent=TRUE, increment=100, 
                    undef.zero=TRUE, warn.braycurtis=TRUE){
						
	# the basic way the data is passed is this:
	# simprof -> simprof.body -> diveDeep (splits into 'left' and 'right' subtrees -> simprof.body (same as before)
	# if you change anything (to pass it down the line), make sure you change the arguments of each of the above
  # TODO: pass arguments with ... instead of explicitly.
	
	if (!is.matrix(data)) # TODO: this could probably be re-worked to not be matrix-dependent
		data <- as.matrix(data) ### make it consistent handling of the data

	rawdata<-data # the user doesn't need to see "rawdata", but changing a lot of code isn't desirable right now
	### the code is written for rows to be the samples/sites, so if it is the other way take the transpose
	if(sample.orientation == "column")
		rawdata<-t(rawdata)
	### now we should either have a function here or something to pass to dist()	
	if (is.function(method.distance)) # allow arbitrary distance function choice
		rawdata.dist <- method.distance(rawdata)
	else if (method.distance == "braycurtis"){
		if (warn.braycurtis){
      warning("This version of the Bray-Curtis index does not use standardization.", call.=FALSE)
      warning("To use the standardized version, use \"actual-braycurtis\".", call.=FALSE)
      warning("See the help documentation for more information.", call.=FALSE)
		  }
    rawdata.dist <- braycurtis(rawdata, const, undef.zero)
		# the next bit takes care of putting row names back onto the rawdata.dist
		if (!is.null(rownames(rawdata))){
			attr(rawdata.dist, "Labels") <- rownames(rawdata)
		  }
		}
	else if (method.distance == "czekanowski"){
    rawdata.dist <- czekanowski(rawdata, const, undef.zero) # this is IDENTICAL to the braycurtis() function
    if (!is.null(rownames(rawdata))){
      attr(rawdata.dist, "Labels") <- rownames(rawdata)
    }
	}
	else if (method.distance == "actual-braycurtis"){
    rawdata.dist <- braycurtis(preBCstandardization(rawdata),const, undef.zero)
    if (!is.null(rownames(rawdata))){
      attr(rawdata.dist, "Labels") <- rownames(rawdata)
    }
	}
	else {
		rawdata.dist <- dist(rawdata, method=method.distance)
	}
	### transforming the data
	if (!method.transform == "identity")
		rawdata <- trans(rawdata, method.transform)
		
	hclust.results <- hclust(rawdata.dist, method=method.cluster) ### Generate the cluster information
	pMatrix <- cbind(hclust.results$merge, matrix(data = NA, nrow=nrow(hclust.results$merge), ncol=1))
	
	### Namely one that includes a vector indicating which significant cluster each sample belongs to
	simprof.results <- simprof.body(rawdata=rawdata, num.expected=num.expected, num.simulated=num.simulated, 
			method.cluster=method.cluster, method.distance=method.distance, 
			originaldata=rawdata, alpha=alpha, clust.order=hclust.results$merge, 
			startrow=nrow(hclust.results$merge), pMatrix=pMatrix, 
      const=const, silent=silent, increment=increment, undef.zero)
	
	results <- list()
	results[["numgroups"]] <- length(simprof.results$samples) # number of significant groups
	results[["pval"]] <- simprof.results$pval
	results[["hclust"]] <- hclust.results
	if(!is.null(rownames(rawdata))){
		for (i in 1:length(simprof.results$samples)){
			for(j in 1:length(simprof.results$samples[[i]])){
				simprof.results$samples[[i]][j] <- rownames(rawdata)[as.numeric(simprof.results$samples[[i]][j])]
				}
			}
		}
	results[["significantclusters"]] <- simprof.results$samples
	return(results)
	
	}
	
simprof.body <- function(rawdata, num.expected=1000, num.simulated=999, 
					method.cluster="average", method.distance="euclidean", originaldata=NA,
					alpha=0.05, clust.order=NA, startrow, pMatrix, currentsamples=NA, const=const, 
					silent, increment, undef.zero){

	### The below description is incorrect. Ignore and change.
	### Find all of the basic data we'll need to form a test statistic
	### To do this, first find all of the pairwise distances between sites/samples. Rank them.
	### Second, permute each column separately. Find all of the new pairwise distances. Rank them.
	### Third, repeat this num.expected (1000) times.
	### Fourth, compute the sum of the absolute values of the distances for each rank between the original data and the expected profile.
	### Fifth, repeat steps 2 num.simulated (999) times. This is the simulated data under the null for comparison.
	### Sixth, compute the test statistic between each of these profiles and the expected.
	### Seventh, compute the p-value given the original data and the num.simulated (999) simulated values.
	### Finally, if the p-value is less than alpha (0.05), run simprof on the left and right sub-clusters.
	
	### Real data simprof
 	rawdata.simprof <- list(genSimilarityProfile(rawdata, method.distance, const, undef.zero))
	if (!dim(rawdata.simprof[[1]])[2]==1) # order screws up the list if it is only a column (and if it is only a single column, there is no need to order)
		rawdata.simprof[[1]] <- rawdata.simprof[[1]][,order(rawdata.simprof[[1]][2,])]
	### maybe this is the slowdown point? 
	
	### Expected Profile
	expectedprofile.simprof <- genProfile(rawdata, originaldata, num.expected, 
                                        method.distance, const, silent=silent, 
                                        increment=increment, type="Expected",
                                        undef.zero)
	expectedprofile.average <- computeAverage(expectedprofile.simprof, num.expected)

	### Test Statistic
	teststatistic <- computeTestStatistic(rawdata.simprof[[1]], expectedprofile.average)
   
	### Simulated Profile
	simulatedprofile.simprof <- genProfile(rawdata, originaldata, num.simulated, 
                                         method.distance, const, silent=silent, 
                                         increment=increment, type="Simulated",
                                         undef.zero)
	
	### Comparison
	pval <- tsComparison(simulatedprofile.simprof, expectedprofile.average, num.simulated, teststatistic)
   
	findings <- c();
	findings[["pval"]] <- pTracker(pMatrix, startrow, pval)
	
	if (pval > alpha){
		if (!is.na(currentsamples[1])){ ### there are some warnings here... find out why
										### warnings should stop if we reference a specific index
										### none should be NA if currentsamples != NA
										# pop out and start working our way back
			findings[["samples"]] <- list(currentsamples)
			return(findings)
			}
		else{ # we failed to reject the first null: all samples are one group
			findings[["samples"]] <- list(c(1:nrow(rawdata)))
			return(findings)
			### rework the above line (eventually) to be placed in an order determined by hclust()
			### ... maybe 
			}
		}
	
	### Need to cut down rawdata to only be the next set of samples

	### See if we need to investigate subtrees
	if (pval <= alpha){ ## should this be strict inequality?
		
		simprof.left  <- diveDeep(rawdata=rawdata, num.expected=num.expected, num.simulated=num.simulated, 
				method.cluster=method.cluster, method.distance=method.distance,
				originaldata=originaldata, alpha=alpha, clust.order=clust.order, startrow=startrow, pMatrix=findings$pval,
				side="LEFT", const=const, silent=silent, increment=increment)
		simprof.right <- diveDeep(rawdata=rawdata, num.expected=num.expected, num.simulated=num.simulated, 
				method.cluster=method.cluster, method.distance=method.distance,
				originaldata=originaldata, alpha=alpha, clust.order=clust.order, startrow=startrow, pMatrix=findings$pval,
				side="RIGHT", const=const, silent=silent, increment=increment)

		findings[["samples"]] <- append(simprof.left$samples, simprof.right$samples)
		findings[["pval"]] <- pTrackerMerge(simprof.left$pval, simprof.right$pval)
		return(findings)
		}
	}
