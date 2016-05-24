#' Computes the ranks of all the algorithms from their (repeated) results measurements after
#' grouping them by several factors combined simultaneosly.
#'
#' @param data A dataframe object containing (at least) two columns for the target factor and the performance measure
#'  Additional columns are aimed at grouping the problem configuration by (at most) 3 different factors.
#' @param params A vector with the column names in \code{data} that define a problem configuration. If not already factor objects, those columns will be converted to
#'	factors inside the function (note this does not alter the ordering of the levels in case it was explicitly set before the call).
#'  Although an arbitrary number of columns can be passed, if the user intends to plot the ranks computed by this function, at most three columns should be passed.
#' @param target Name of the target column of \code{data}. For each combination of the values of \code{params}, the ranks are obtained by
#'	comparing the repeated measurements of \code{performance} associated to each level of the \code{target} column.
#' @param performance Name of the column of \code{data} containing the repeated performance measurements. If given a vector of strings, 
#'	then a separate ranking will be computed for each of the elements, and no p-values, mean or stdev columns will be returned, just the rankings together with the factors 
#'	to indicate which problem configuration corresponds to the rank.
#' @param pairing.col Name of the column which links together the paired samples, in case we have set \code{paired = TRUE}. Otherwise, this argument will be ignored.
#' @param test The statistical test to be performed to compare the performance of every level of the target variable at each problem configuration. 
#' @param fun Function performing a custom statistical test, if \code{test = "custom"}; otherwise, this argument is ignored. The function must receive exactly
#'	two vectors (the first is a vector of real numbers and the second is a factor with the level to which each real number corresponds) 
#'  and must return a \code{pairwise.htest} object with a \code{p.value} field. This must be an (N-1)x(N-1) lower-triangular matrix, with exactly the same structure 
#'	as those returned in the \code{p.value} field by a call to \code{\link{pairwise.wilcox.test}} or \code{\link{pairwise.t.test}}.
#' @param correction The p-value adjust method. Must be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" (defaults to "holm").
#'	This parameter will be ignored if \code{test = "tukeyHSD"} as Tukey HSD incorporates its own correction procedure.
#' @param alpha Significance threshold for pairwise comparisons. Defaults to 0.05.
#' @param maximize Boolean indicating whether the higher the performance measure, the better (default), or vice-versa.
#' @param ncores Number of physical CPUs available for computations. If \code{ncores} > 1, parallelization is achieved through the \code{\link{parallel}} package and
#' 	is applied to the computation of ranks for more than one problem configuration at the same time. Defaults to 1 (sequential).
#' @param paired Boolean indicating whether samples in the same problem configuration, which only differ in the target value, and in the same relative position (row) within their
#' 	respective target values are paired or not. Defaults to FALSE. This should be set to TRUE, for instance, in Machine Learning problems in which, for a fixed problem configuration,
#' 	the target variable (usually the algorithms being compared) is associated to a number of samples (results) coming from the Cross Validation process. If a K-fold CV is being done,
#' 	then we would have, for a given problem configuration, K rows for each of the algorithms being compared, all of them identical in all the columns except for the performance column. 
#' 	In that case, the performance of the i-th row (1 <= i <= K) of all of those batches (groups of K rows) for that fixed problem configuration would be related, 
#' 	hence every pairwise comparison should take into account paired samples.
#' @param ... Further arguments to be passed to the function \code{fun} that is called for every pairwise comparison.
#' @return If \code{length(performance)} equals 1, an object of classes \code{c("SRCS", "data.frame")} with the following columns: 
#' - A set of columns with the same names as the \code{params} and \code{target} arguments.
#' - Two columns called "mean" and "sd" containing the mean of the repeated peformance measurements for each problem configuration and the standard deviation.
#' - One column named "rank" with the actual rank of each level of the \code{target} variable within that problem configuration. The lower the rank, the better the algorithm.
#' - |target| additional columns containing the p-values resulting of the comparison between the algorithm and the rest for the same problem configuration,
#' 	where |target| is the number of levels of the target variable.
#'
#' If \code{length(performance)} > 1 (let \code{P = length(performance)} for the explanation that follows), an object of classes \code{c("SRCS","data.frame")}
#'	with the following columns:
#' - A set of columns with the same names as the \code{params} and \code{target} arguments.
#' - One column per element of the \code{performance} vector, named "rank1", ..., "rankP", containing, for each performance measure, 
#'	the rank of each level of the \code{target} variable within that problem configuration for that performance measure. 
#'	The higher the rank, the better the algorithm.
#' @note Although it has no effect on the results of \code{SRCSranks}, the user should preferably have set the order 
#'	of the factor levels explicitly by calling function \code{levels} before calling this function, specially if he intends to subsequently apply \code{\link{plot}} to the results,
#'	because the level order does affect the way graphics are arranged in the plot.
#' @seealso \code{\link{plot.SRCS}} for a full working example of \code{SRCSranks} and plotting facilities. Also 
#'	\code{\link{pairwise.wilcox.test}}, \code{\link{t.test}}, \code{\link{pairwise.t.test}}, \code{\link{TukeyHSD}}, \code{\link{p.adjust.methods}}
SRCSranks <- function(data, params, target, performance, pairing.col = NULL, test = c("wilcoxon", "t", "tukeyHSD", "custom"), fun = NULL, correction = p.adjust.methods, 
					alpha = 0.05, maximize = TRUE, ncores = 1, paired = FALSE, ...){

	test = match.arg(test);
	correction = match.arg(correction);
	
	colNames = names(data);
	discarded.names = colNames[!(colNames %in% c(params, target))];
	nrows = length(data);
	param.list = c(params, target, performance, pairing.col);
	for(param in param.list){
		if(!is.null(param)){
			if(sum(param %in% colNames) < length(param)){
				stop(paste("ERROR:",param[which.min(param %in% colNames)],"not found in the frame data"));
			}
		}
  }

  if(is.null(target) || is.null(performance)){
		stop("ERROR: the target column and the performance column cannot be null");
  }
  
  if(paired == TRUE && is.null(pairing.col)){
    stop("ERROR: a pairing column must always be indicated when setting paired = TRUE");
  }
    
  # ------------------------------------
  ## CHECK ALL COLUMNS ARE FACTORS
  # ------------------------------------
  for(param in c(params,target,pairing.col)){
    if(!is.null(param) && !is.factor(data[[param]])){	
        # warning("column ",param, " had to be converted into a factor"); 
        data[[param]] = as.factor(data[[param]]); 
    }
  }
  
  # ------------------------------------

  tempData = data[,(colNames %in% params)];  
  chunked.frames = split(data, tempData, drop = TRUE);

  ## Compose the arguments to either lapply or parLapply applied to chunked.frames  
  arglist = list(X = chunked.frames, FUN = function(chframe){
      if(paired == TRUE){
        myorder = order(chframe[[target]], chframe[[pairing.col]]);
        chframe = chframe[myorder,];
      }
			if(length(performance) == 1){ 
																		rankList = .rank(chframe[[performance]], group = chframe[[target]], test = test, fun = fun, correction = correction, 
																										alpha = alpha, maximize = maximize, paired = paired, ...); 
			}
			else{                         rankList = .rankall(chframe[performance], group = chframe[[target]], test = test, fun = fun, correction = correction,
																										alpha = alpha, maximize = maximize, paired = paired, ...); 
			}
			v = names(rankList);
			v[v=="target"] = target;
			names(rankList) = v;
			## Drop the performance column and delete repeated rows ...
			chframe = unique(chframe[,!(names(chframe) %in% discarded.names)]);
			## ... and join the rank column			
			return(merge(chframe, rankList));
   });
  
	if(ncores > 1){
    # ------------------------------------------------------------
    #           PARALLEL
    # ------------------------------------------------------------
    detected.cores = detectCores();
    ncores = min(ncores, detected.cores);
    cl = makeCluster(ncores);
    clusterExport(cl,list(".rank", ".rankall"), envir=environment());
    
		arglist[["cl"]] = cl;
		names(arglist)[names(arglist) == "FUN"] = "fun";	# the correct name of this argument is "fun" for parLapply
		chunked.frames = do.call(parLapply, arglist);

		stopCluster(cl);
  }
  else{  
    # ------------------------------------------------------------
    #           SEQUENTIAL
    # ------------------------------------------------------------  
		chunked.frames = do.call(lapply, arglist);
  }

	result.frame = do.call(rbind, chunked.frames);
	rownames(result.frame) = NULL; 
	
	class(result.frame) = c("SRCS","data.frame");
	
  return(result.frame);
}

# ___________________________________________________________________________

#' Compares the performance of two algorithms for a single problem configuration specified by the user.
#' @param rankdata The ranks data frame obtained by a previous call to \code{\link{SRCSranks}}.
#' @param target Name of the target column in \code{rframe} that separates the levels to be compared, probably "Algorithm" or similar.
#' @param alpha Significance threshold to consider two set of measurements coming from two algorithms as statistically significant
#' @param pvalues Boolean. TRUE indicates that the pairwise comparison table should contain p-values. 
#' 	FALSE means only ">","<" or "=" (the latter for non-significant difference) will be displayed in the table. Defaults to FALSE.
#' @param ... The rest of the columns in \code{rframe} and the values to fully specify a single problem configuration for which algorithms will be compared.
#'  Must be indicated as named arguments, like in "severity" = 4. 
#' @return A square matrix of the same dimension as algorithms found in the data. An entry i,j contains either the p-value of the Wilcoxon test between 
#'	algorithms i and j (if \code{pvalues} was set to TRUE), or the qualitative result (">", "<" or "=") of the statistical comparison (if \code{pvalues} was set to FALSE).
#' @seealso \code{\link{SRCSranks}, \link{plot.SRCS}} for a full working example of SRCScomparison.
SRCScomparison <- function(rankdata, target, alpha = 0.05, pvalues = FALSE, ...){
	params = list(...);
	nms = names(params);
	if(length(target) > 1){
		stop("the target column must be a single name");
	}
	if(!(target %in% names(rankdata))){
		stop("target column not found on this object");
	}
	if(sum(!(nms %in% names(rankdata))) > 0){
		bad = params[[which.min(nms %in% names(rankdata))]];
		stop(paste0("column ", bad, " was not found on this object"));
	} 
	vbool = rep(TRUE, nrow(rankdata));
	
	# Create boolean vector of rows that fulfill all conditions simultaneously
	for(i in 1:length(params)){
		temp = (rankdata[[ nms[i] ]] == params[[i]]);
		vbool = vbool & temp;
	};
	
	if(sum(vbool) == 0){
		stop("the data subset meeting all the specified column values is empty");
	}
	
	rsubset = rankdata[vbool, ];
	rownames(rsubset) = rsubset[[target]];
	targetval = unique(rankdata[,target]);
	if(length(targetval) != nrow(rsubset)){
		stop("too few values provided: the <column = value> pairs specified do not restrict the target sufficiently.\nTry passing more <column = value> pairs as arguments");
	}
	
	comparisons = matrix("", nrow = nrow(rsubset), ncol = nrow(rsubset));
	pval.mat = matrix(NA, nrow = nrow(rsubset), ncol = nrow(rsubset));
	rownames(comparisons) = targetval; 
	colnames(comparisons) = targetval;
	rownames(pval.mat) = targetval; 
	colnames(pval.mat) = targetval;
	
	for(alg1 in targetval){
    for(alg2 in targetval){
      if(!is.na(rsubset[alg1, paste0(alg2,".pval")])){
        if(rsubset[alg1, paste0(alg2,".pval")] < alpha){        
          ## Note we are comparing the performance measure, disregarding whether it is to be maximized or minimized
          ## It is the user who has to interpret what a "<" and a ">" mean
          if(rsubset[alg1, "mean"] < rsubset[alg2, "mean"]){ comparisons[alg1, alg2] = "<"; }
          else{ 																						 comparisons[alg1, alg2] = ">"; }
        }
        else{
          comparisons[alg1, alg2] = "=";
        }
      }
      else{
        comparisons[alg1, alg2] = NA;
      }
      pval.mat[alg1, alg2] = rsubset[alg1, paste0(alg2,".pval")];
    }
	}
	if(pvalues){ return(pval.mat);    }
	else       { return(comparisons);	}
}

# ___________________________________________________________________________
