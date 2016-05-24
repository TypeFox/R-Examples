# _________________________________________________________________________
## This function actually computes the rank of a set of algorithms using
## their performance results on a *single* problem configuration.
## It is based on the code published in the book chapter by I.G.Del Amo and D.Pelta
##
## @param data A vector with the results of the m algorithms, that is, the \code{n}
##   independent repetitions of the performance measure, for each algorithm.
##   Since we will be performing a paired test, the data must be carefully arranged
##   so that \code{data[j]} (j-th datum of the first group) is paired with all \code{data[j+Kn]}, K = 1,2,...,m
##   (these are the j-th datum of the remaining groups, where \code{n} is the number of samples of each group).
## @param group A factor vector for the data, indicating, for each corresponding
##     entry of \code{data}, the algorithm it belongs to.
## @param alpha Significance threshold for the tests. Defaults to 0.05
## @param max Boolean: should the performance measure be maximized? 
##     TRUE means that the higher the performance measure,
##     the better. FALSE means the lower, the better. Defaults to TRUE.
## @return 	A dataframe object where every row refers to an algorith, and with the following columns:
##	  First two columns are the algorithm name and its rank, named "target" and "rank"
##	  Second two columns are the mean and std.dev. of the performance of the algorithm, named "mean" and "sd".
##    Finally: as many columns as algorithms, storing the p-value of every pairwise comparison of this algorithm with the rest,
##		named as the algorithm with ".pval" at the end.
## @examples
## dat <- c( 2.5, 2.3, 1.2, 0.7, 3.5, 4.1 )
## group <- factor(c("alg1", "alg1", "alg2", "alg2", "alg3", "alg3"))
## @references I.G. del Amo and D.A. Pelta. 
##  SRCS: A Technique for Comparing Multiple Algorithms under Several Factors in Dynamic Optimization Problems
##  In: Metaheuristics for Dynamic Optimization. E. Alba A. Nakib, P. Siarry (Eds). Studies in Computational Intelligence 433, 61-77. 
##  Springer, 2013.
# _________________________________________________________________________
.rank <- function(data, group, test, fun, correction, alpha = 0.05, maximize = TRUE, paired = FALSE, ...)
{
	algorithms <- unique(group);	
	nalgorithms = length(algorithms);
  rankFrame = data.frame(matrix(NA, nrow = length(algorithms), ncol = 2 + 2 + length(algorithms))); 
  
  rankFrame[,1] = as.factor(unlist(algorithms));

  pvalColumnNames = paste0(algorithms,".pval");
  names(rankFrame) = c("target","rank","mean","sd",pvalColumnNames);
  rankFrame$rank = 0;

	alg.bool = list();
	for(algorithm in algorithms){
		v = (rankFrame$target == algorithm);
		alg.bool[[algorithm]] = which(v, arr.ind = TRUE);
	}

	# calculate the vector of means for all the algorithms' measures
	means <- tapply(data, group, mean);
	stdevs <- tapply(data, group, sd);

	# perform a Kruskal-Wallis test to assess if there are differences
	# among all the results
	dataframe <- data.frame(group, data)
	kruskal <- kruskal.test(data ~ group, data = dataframe)
	if (!is.na( kruskal$p.value) && kruskal$p.value < alpha)
	{
		# post-hoc test: perform the indicated pairwise test with the desired correction method to assess individual differences
		test.results = NULL;
		
		## ---------------------------------
		## APPLY THE SELECTED POST-HOC TEST
		## ---------------------------------
		if(test == "wilcoxon"){
			test.results = pairwise.wilcox.test(data, group, p.adjust.method = correction, exact = FALSE, paired = paired);
		}
		if(test == "t"){
			test.results = pairwise.t.test(data, group, p.adjust.method = correction, exact = FALSE, paired = paired);
		}
		if(test == "tukeyHSD"){
			group2 = as.factor(gsub("-", "@", group)); # to avoid the case that the algorithm names already contain "-"
			levels(group2) = gsub("-", "@", levels(group));
			tempgroup = group2;
			mypvalues = TukeyHSD(aov(data ~ tempgroup))$tempgroup;
			grouplevels = levels(group);
			
			## Manually compose the p.values matrix
			pvalue.matrix = matrix(NA, nrow = length(algorithms)-1, ncol = length(algorithms)-1);						
			colnames(pvalue.matrix) = grouplevels[1:(nlevels(group)-1)];
			rownames(pvalue.matrix) = grouplevels[2:nlevels(group)];
			
			separated = strsplit(rownames(mypvalues), "-");
			for(myrow in 1:nrow(mypvalues)){
				algs = gsub("@", "-", separated[[myrow]]);
				pvalue.matrix[ algs[1], algs[2] ] = mypvalues[myrow, "p adj"];
			}
			test.results = list(p.value = pvalue.matrix);
		}
		
		if(test == "custom"){
			myargs = list(c(data, group), as.list(...));
			test.results = do.call(fun, myargs);
		}
				
		for (algorithm1 in rownames(test.results$p.value))
		{
			for (algorithm2 in colnames(test.results$p.value))
			{
				if (!is.na( test.results$p.value[algorithm1,algorithm2])){
						if(test.results$p.value[algorithm1,algorithm2] < alpha){
						# there is a significant difference between algorithm1 and algorithm2;
						# we need to identify which one is the best and which one the worst,
						# we'll use the median for that purpose
						if (means[algorithm1] > means[algorithm2])
						{
							best <- algorithm1
							worst <- algorithm2
						}
						else
						{
							best <- algorithm2
							worst <- algorithm1
						}

						if (!maximize)
						{
							# swap best and worst
							tmp <- best
							best <- worst
							worst <- tmp
						}

						# update ranks
						#rankFrame[rankFrame$target == best, 2] = rankFrame[rankFrame$target == best, 2] + 1;
						#rankFrame[rankFrame$target == worst, 2] = rankFrame[rankFrame$target == worst, 2] - 1;						
						rankFrame[alg.bool[[best]], 2] = rankFrame[alg.bool[[best]], 2] + 1;
						rankFrame[alg.bool[[worst]], 2] = rankFrame[alg.bool[[worst]], 2] - 1;
					}
					# store p-values of this comparison
					#rankFrame[rankFrame$target == algorithm1, paste0(algorithm2,".pval")] = test.results$p.value[algorithm1,algorithm2];
					#rankFrame[rankFrame$target == algorithm2, paste0(algorithm1,".pval")] = test.results$p.value[algorithm1,algorithm2];
					rankFrame[alg.bool[[algorithm1]], paste0(algorithm2,".pval")] = test.results$p.value[algorithm1,algorithm2];
					rankFrame[alg.bool[[algorithm2]], paste0(algorithm1,".pval")] = test.results$p.value[algorithm1,algorithm2];
				}
			}
		}
		for(algorithm in algorithms){
			#store the mean and sd
			#rankFrame[rankFrame$target == algorithm, 3] = means[[algorithm]];
			#rankFrame[rankFrame$target == algorithm, 4] = stdevs[[algorithm]];		
			rankFrame[alg.bool[[algorithm]], 3] = means[[algorithm]];
			rankFrame[alg.bool[[algorithm]], 4] = stdevs[[algorithm]];					
		}
	}
	return(rankFrame)
}

# ___________________________________________________________________________
## Same as .rank() but for several different performance columns. It ranks every column separately
##
## Argument data is not a vector of data, but a matrix of data where each COLUMN represents
## a different performance measure. nrows(data) must equal length(group).
##
## Returns a dataframe object with 1+ncol(data). First column is named "target" and contains unique(group). 
## The rest are named "rank1", "rank2", ... "rankn", with n = ncol(data), i.e. the number of performance measures
# _________________________________________________________________________
.rankall <- function(alldata, group, test, fun, correction, alpha = 0.05, maximize = TRUE, paired = FALSE, ...){	
	algorithms <- unique(group);	
	nranks = ncol(alldata);
  rankFrame = data.frame(matrix(NA, nrow = length(algorithms), ncol = 1 + nranks)); 
  
  rankFrame[,1] = as.factor(unlist(algorithms));
  
  names(rankFrame) = c("target",paste0("rank",seq(1:nranks)));
  rankFrame[,2:(nranks + 1)] = 0;

	alg.bool = list();
	for(algorithm in algorithms){
		v = (rankFrame$target == algorithm);
		alg.bool[[algorithm]] = which(v, arr.ind = TRUE);
	}

	means = sapply(X = names(alldata), FUN = function(algname) tapply(X = alldata[[algname]], INDEX = group, FUN = mean));

	# perform a Kruskal-Wallis test to assess if there are differences
	# among all the results
	for(j in 1:nranks){
    data = alldata[,j];
    dataframe <- data.frame(group, data)
    names(dataframe) = c("group", "data");
    kruskal <- kruskal.test(data ~ group, data = dataframe)
    if (!is.na( kruskal$p.value) && kruskal$p.value < alpha)
    {      
      ## ---------------------------------
			## APPLY THE SELECTED POST-HOC TEST
			## ---------------------------------
			if(test == "wilcoxon"){
				test.results = pairwise.wilcox.test(data, group, p.adjust.method = correction, exact = FALSE, paired = paired);
			}
			if(test == "t"){
				test.results = pairwise.t.test(data, group, p.adjust.method = correction, exact = FALSE, paired = paired);
			}
			if(test == "tukeyHSD"){
				group2 = as.factor(gsub("-", "@", group)); # to avoid the case that the algorithm names already contain "-"
				levels(group2) = gsub("-", "@", levels(group));
				tempgroup = group2;
				mypvalues = TukeyHSD(aov(data ~ tempgroup))$tempgroup;
				grouplevels = levels(group);
				
				## Manually compose the p.values matrix
				pvalue.matrix = matrix(NA, nrow = length(algorithms)-1, ncol = length(algorithms)-1);						
				colnames(pvalue.matrix) = grouplevels[1:(nlevels(group)-1)];
				rownames(pvalue.matrix) = grouplevels[2:nlevels(group)];
				
				separated = strsplit(rownames(mypvalues), "-");
				for(myrow in 1:nrow(mypvalues)){
					algs = gsub("@", "-", separated[[myrow]]);			
					pvalue.matrix[ algs[1], algs[2] ] = mypvalues[myrow, "p adj"];
				}
				test.results = list(p.value = pvalue.matrix);
			}
			
			if(test == "custom"){
				myargs = list(c(data, group), as.list(...));
				test.results = do.call(fun, myargs);
			}     
			      
      for (algorithm1 in rownames(test.results$p.value))
      {
        for (algorithm2 in colnames(test.results$p.value))
        {
          if (!is.na( test.results$p.value[algorithm1,algorithm2])){
              if(test.results$p.value[algorithm1,algorithm2] < alpha){
              # there is a significant difference between algorithm1 and algorithm2;
              # we need to identify which one is the best and which one the worst,
              # we'll use the median for that purpose
              if (means[algorithm1,j] > means[algorithm2,j])
              {            
                best <- algorithm1
                worst <- algorithm2
              }
              else
              {
                best <- algorithm2
                worst <- algorithm1
              }

              if (!maximize)
              {
                # swap best and worst
                tmp <- best
                best <- worst
                worst <- tmp
              }
              
              # update ranks
              rankFrame[alg.bool[[best]],  1+j] = rankFrame[alg.bool[[best]], 1+j] + 1;
              rankFrame[alg.bool[[worst]], 1+j] = rankFrame[alg.bool[[worst]], 1+j] - 1;
            }
          }
        }
      }
    }
	};
	return(rankFrame)
}