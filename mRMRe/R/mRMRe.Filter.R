## Definition

setClass("mRMRe.Filter", representation(filters = "list", scores = "list", mi_matrix = "matrix", causality_list = "list",
                sample_names = "character", feature_names = "character", target_indices = "integer", levels = "integer"))

## Wrappers

`mRMR.ensemble` <- function(solution_count, feature_count, ...)
{
    return(new("mRMRe.Filter", levels = c(solution_count, rep(1, feature_count - 1)), ...))
}

`mRMR.classic` <- function(feature_count, ...)
{
    return(new("mRMRe.Filter", levels = rep(1, feature_count), ...))
}

## initialize

setMethod("initialize", signature("mRMRe.Filter"),
        function(.Object, data, prior_weight, target_indices, levels, 
          method = c("exhaustive", "bootstrap"), 
          continuous_estimator = c("pearson", "spearman", "kendall", "frequency"), 
          outX = TRUE,
          bootstrap_count = 0)
{
    method <- match.arg(method)
    continuous_estimator <- match.arg(continuous_estimator)
    
    if (class(data) != "mRMRe.Data")
        stop("data must be of type mRMRe.Data")
    
    ## Prior Processing
    
    if (length(priors(data)) != 0)
    {
        if (missing(prior_weight))
            stop("prior weight must be provided if there are priors")
        else if  (prior_weight < 0 || prior_weight > 1)
            stop("prior weight must be a value ranging from 0 to 1")
    }
    else
        prior_weight <- 0
    
    ## Target Processing

    if (sum(sapply(target_indices, function(index) index < 1 || index > featureCount(data))) > 1)
        stop("target_indices must only contain values ranging from 1 to the number of features in data")
            
    ## Level Processing
    
    if (missing(levels))
        stop("levels must be provided")
    else if ((prod(levels) - 1) > choose(featureCount(data) - 1, length(levels)))
        stop("user cannot request for more solutions than is possible given the data set")
    
    .Object@target_indices <- as.integer(c(target_indices))
    .Object@levels <- as.integer(c(levels))
    
    target_indices <- as.integer(.expandFeatureIndices(data, target_indices)) - 1
    
    ## Filter; Mutual Information and Causality Matrix

    mi_matrix <- as.numeric(matrix(NA, ncol = ncol(data@data), nrow = ncol(data@data)))
    
	if(method == "exhaustive")
    	result <- .Call(.C_export_filters, as.integer(.Object@levels), as.numeric(data@data),
        	    as.numeric(data@priors), as.numeric(prior_weight), as.integer(data@strata), as.numeric(data@weights),
            	as.integer(data@feature_types), as.integer(nrow(data@data)), as.integer(ncol(data@data)),
            	as.integer(length(unique(data@strata))), as.integer(target_indices),
            	as.integer(.map.continuous.estimator(continuous_estimator)), as.integer(outX),
            	as.integer(bootstrap_count), mi_matrix)
	else if(method == "bootstrap")
		result <- .Call(.C_export_filters_bootstrap, as.integer(.Object@levels[1]), as.integer(length(.Object@levels)),
				as.numeric(data@data), as.numeric(data@priors), as.numeric(prior_weight), as.integer(data@strata),
				as.numeric(data@weights), as.integer(data@feature_types), as.integer(nrow(data@data)),
				as.integer(ncol(data@data)), as.integer(length(unique(data@strata))), as.integer(target_indices),
				as.integer(.map.continuous.estimator(continuous_estimator)), as.integer(outX),
				as.integer(bootstrap_count), mi_matrix)
	else
		stop("Unrecognized method: use exhaustive or bootstrap")
    
    .Object@filters <- lapply(result[[1]], function(solutions) matrix(.compressFeatureIndices(data, solutions + 1),
                        nrow = length(levels), ncol = prod(levels)))
	
    names(.Object@filters) <- .Object@target_indices
    .Object@causality_list <- result[[2]]
	.Object@scores <- lapply(result[[3]], function(scores) matrix(scores,	nrow = length(levels), ncol = prod(levels)))
	names(.Object@scores) <- .Object@target_indices
    
	cols_to_drop <- duplicated(.compressFeatureIndices(data, seq(ncol(data@data))))
    
    .Object@causality_list <- lapply(result[[2]], function(causality_array) causality_array[!cols_to_drop])
    names(.Object@causality_list) <- .Object@target_indices
    
    .Object@mi_matrix <- .compressFeatureMatrix(data, matrix(mi_matrix, ncol = ncol(data@data), nrow = ncol(data@data)))
    .Object@feature_names <- featureNames(data)
    .Object@sample_names <- sampleNames(data)

    return(.Object)
})

## show

setMethod("show", signature("mRMRe.Filter"), function(object)
{
    str(object)
})

## sampleCount

setMethod("sampleCount", signature("mRMRe.Filter"), function(object)
{
    return(length(object@sample_names))
})

## sampleNames

setMethod("sampleNames", signature("mRMRe.Filter"), function(object)
{
    return(object@sample_names)
})


## featureCount

setMethod("featureCount", signature("mRMRe.Filter"), function(object)
{
    return(length(object@feature_names))
})

## featureNames

setMethod("featureNames", signature("mRMRe.Filter"), function(object)
{
    return(object@feature_names)
})

## solutions

setMethod("solutions", signature("mRMRe.Filter"), function(object, mi_threshold = -Inf, causality_threshold = Inf)
{
    # filters[[target]][solution, ] is a vector of selected features
    # in a solution for a target; missing values denote removed features
            
    filters <- lapply(object@target_indices, function(target_index)
    {
        result_matrix <- object@filters[[as.character(target_index)]]
        causality_dropped <- which(object@causality_list[[as.character(target_index)]] > causality_threshold &
                        !is.na(object@causality_list[[as.character(target_index)]]))
        mi_dropped <- which(-.5 * log(1 - object@mi_matrix[, target_index, drop = TRUE]^2) < mi_threshold)
        result_matrix[result_matrix %in% c(causality_dropped, mi_dropped)] <- NA

        pre_return_matrix <- apply(as.matrix(result_matrix), 2, rev)

        return(pre_return_matrix)
    })
    names(filters) <- object@target_indices

    return(filters)
})

## scores
setMethod("scores", signature("mRMRe.Filter"), function(object)
{
	mi_matrix <- mim(object)
	targets <- as.character(target(object))
	scores <- lapply(targets, function(target) {
				apply(solutions(object)[[target]], 2, function(solution) {
							sapply(1:length(solution), function(i) {
										feature_i <- solution[i] 
										if(i == 1)
											return(mi_matrix[as.numeric(target), feature_i])
										
										ancestry_score <- mean(sapply((i-1):1, function(j) mi_matrix[feature_i, solution[j]]))
										return(mi_matrix[as.numeric(target), feature_i] - ancestry_score)
									})
							
						})
			})
	names(scores) <- targets
	return(scores)
})

## mim

setMethod("mim", signature("mRMRe.Filter"), function(object, method=c("mi", "cor"))
{
    method <- match.arg(method)
    # mi_matrix[i, j] contains the biased correlation between
    # features i and j (i -> j directionality)
    return(object@mi_matrix)
})

## causality

setMethod("causality", signature("mRMRe.Filter"), function(object)
{
    # causality_matrix[[target]][feature] contains the causality coefficient
    # between feature and target (feature -> target directionality)
    return(object@causality_list)
})
    
## target

setMethod("target", signature("mRMRe.Filter"), function(object)
{
    return(object@target_indices)
})
