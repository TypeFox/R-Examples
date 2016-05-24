## Definition

setClass("mRMRe.Data", representation(sample_names = "character", feature_names = "character", feature_types = "numeric", data = "matrix", strata = "numeric", weights = "numeric", priors = "matrix"))

## Wrapper

`mRMR.data` <- function(...)
{
    return(new("mRMRe.Data", ...))
}

## initialize

setMethod("initialize", signature("mRMRe.Data"), function(.Object, data, strata, weights, priors)
{
    ## Data Processing
    
    if (!is.data.frame(data))
        stop("data must be of type data frame")
        
    if(ncol(data) > (sqrt((2^31) - 1)))
        stop("Too many features, the number of features should be <= 46340!")
    
    feature_types <- sapply(data, function(feature) paste(class(feature), collapse = "_"))
    
    if (any(!is.element(feature_types, c("numeric", "ordered_factor", "Surv"))))
        stop("data columns must be either of numeric, ordered factor or Surv type")

    .Object@sample_names <- rownames(data)
    .Object@feature_names <- colnames(data)
    .Object@feature_types <- unlist(lapply(feature_types, switch, "Surv" = c(2, 3), "ordered_factor" = 1, 0))
    names(.Object@feature_types) <- NULL
    # Optimize the case when all features are continuous
    if(sum(.Object@feature_types) == 0)
	.Object@data <- as.matrix(data) 
    else
        .Object@data <- do.call(cbind, lapply(seq(feature_types), function(i) switch(feature_types[[i]],
                                "Surv" = cbind(event = data[, i][, "status"], time = data[, i][, "time"]),
                                "ordered_factor" = as.numeric(as.integer(data[, i]) - 1),
                                as.numeric(data[, i]))))
    
    rownames(.Object@data) <- rownames(data)
    colnames(.Object@data)[!.Object@feature_types %in% c(2, 3)] <- colnames(data)[feature_types != "Surv"]
    colnames(.Object@data)[.Object@feature_types %in% c(2, 3)] <- paste(rep(colnames(data)[feature_types == "Surv"],
                    each = 2), rep(c("event", "time"), sum(feature_types == "Surv")), sep = "_")
    
    ## Sample Stratum Processing
    
    if (missing(strata)) 
        .Object@strata <- rep.int(0, nrow(data))
    else
        sampleStrata(.Object) <- strata
    
    ## Sample Weight Processing
    
    if (missing(weights)) 
        .Object@weights <- rep(1, nrow(data))
    else
        sampleWeights(.Object) <- weights

    ## Prior Feature Matrix Processing
    
    if (!missing(priors) && !is.null(priors))
        priors(.Object) <- priors
    
    return(.Object)
})

## show

setMethod("show", signature("mRMRe.Data"), function(object)
{
    str(object)
})

## featureData

setMethod("featureData", signature("mRMRe.Data"), function(object)
{
    data <- lapply(seq(object@feature_types), function(i) switch(as.character(object@feature_types[[i]]),
                        "3" = Surv(time = object@data[, i], event = object@data[, i - 1]),
                        "2" = NULL,
                        "1" = object@data[, i] + 1,
                        "0" = object@data[, i],
                        NULL))
    data <- data.frame(data[!sapply(data, is.null)])
    colnames(data) <- object@feature_names
    
    return(data)
})

## subsetData

setMethod("subsetData", signature("mRMRe.Data"), function(object, row_indices, column_indices)
{
    if(missing(row_indices) && missing(column_indices))
        return(object)
        
    if(missing(row_indices))
        row_indices <- 1:sampleCount(object)
    if(missing(column_indices))
        column_indices <- 1:featureCount(object)
    
    data <- featureData(object)[row_indices, column_indices, drop=FALSE]
    strata <- factor(sampleStrata(object)[row_indices])
    weights <- sampleWeights(object)[row_indices]
    priors <- priors(object)
    if(length(priors) > 0)
        priors <- priors[column_indices, column_indices, drop=FALSE]
    else
        priors <- NULL
    
    return(new("mRMRe.Data", data = data, strata = strata, weights = weights, priors = priors))
})

## sampleCount

setMethod("sampleCount", signature("mRMRe.Data"), function(object)
{
    return(nrow(object@data))
})

## sampleNames

setMethod("sampleNames", signature("mRMRe.Data"), function(object)
{
    return(object@sample_names)
})

## featureCount

setMethod("featureCount", signature("mRMRe.Data"), function(object)
{
    return(length(object@feature_names))
})

## featureNames

setMethod("featureNames", signature("mRMRe.Data"), function(object)
{
    return(object@feature_names)
})

## sampleStrata

setMethod("sampleStrata", signature("mRMRe.Data"), function(object)
{
    strata <- object@strata
    names(strata) <- rownames(object@data)
    
    return(strata)
})

## sampleStrata<-

setReplaceMethod("sampleStrata", signature("mRMRe.Data"), function(object, value)
{
    if (length(value) != nrow(object@data))
        stop("data and strata must contain the same number of samples")
    else if (!is.factor(value))
        stop("strata must be provided as factors")
    else if (sum(is.na(value)) > 0)
        stop("cannot have missing values in strata")
    else
        object@strata <- as.integer(value) - 1
    
    return(object)
})

## sampleWeights

setMethod("sampleWeights", signature("mRMRe.Data"), function(object)
{
    weights <- object@weights
    names(weights) <- rownames(object@data)
    
    return(weights)
})

## sampleWeights<-

setReplaceMethod("sampleWeights", signature("mRMRe.Data"), function(object, value)
{
    if (length(value) != nrow(object@data))
        stop("data and weight must contain the same number of samples")
    else if (sum(is.na(value)) > 0)
        stop("cannot have missing values in weights")
    else
        object@weights <- as.numeric(value)
    
    return(object)
})

## priors

setMethod("priors", signature("mRMRe.Data"), function(object)
{
    if (length(object@priors) == 0)
        return(object@priors)
    else
        return(.compressFeatureMatrix(object, object@priors))
})

## priors<-

setReplaceMethod("priors", signature("mRMRe.Data"), function(object, value)
{
    if (ncol(value) != ncol(object@data) || nrow(value) != ncol(object@data))
        stop("priors matrix must be a symmetric matrix containing as many features as data")
    else
        object@priors <- .expandFeatureMatrix(object, value)
    
    return(object)
})

## mim

setMethod("mim", signature("mRMRe.Data"),
        function(object, prior_weight = 0, continuous_estimator = c("pearson", "spearman", "kendall", "frequency"), outX = TRUE, bootstrap_count = 0)
{
    continuous_estimator <- match.arg(continuous_estimator)
    if (length(object@priors) != 0)
    {
        if (missing(prior_weight))
            stop("prior weight must be provided if there are priors")
        else if  (prior_weight < 0 || prior_weight > 1)
            stop("prior weight must be a value ranging from 0 to 1")
    }
    else
        prior_weight <- 0
    
    mi_matrix <- as.numeric(matrix(NA, ncol = ncol(object@data), nrow = ncol(object@data)))
    
    .Call(.C_export_mim, as.numeric(object@data), as.numeric(object@priors),
            as.numeric(prior_weight), as.integer(object@strata), as.numeric(object@weights),
            as.integer(object@feature_types), as.integer(nrow(object@data)), as.integer(ncol(object@data)),
            as.integer(length(unique(object@strata))),
            as.integer(.map.continuous.estimator(continuous_estimator)),
            as.integer(outX), as.integer(bootstrap_count), mi_matrix)
    
    mi_matrix <- matrix(mi_matrix, ncol = ncol(object@data), nrow = ncol(object@data))
    
    mi_matrix <- .compressFeatureMatrix(object, mi_matrix)

    # mi_matrix[i, j] contains the biased correlation between
    # features i and j (i -> j directionality)
    
    return(mi_matrix)
})

## expandFeatureMatrix

setMethod(".expandFeatureMatrix", signature("mRMRe.Data"), function(object, matrix)
{
    adaptor <- which(object@feature_types != 3)
    matrix <- do.call(cbind, lapply(seq(adaptor), function(i)
    {
        column <- do.call(rbind, lapply(seq(adaptor), function(j)
        {
            item <- matrix[j, i]
            
            if (object@feature_types[[adaptor[[j]]]] == 2)
                return(rbind(item, item, deparse.level = 0))
            else
                return(item)
        }))
        
        if (object@feature_types[[adaptor[[i]]]] == 2)
            return(cbind(column, column, deparse.level = 0))
        else
            return(column)
    }))

    return(matrix)
})

## compressFeatureMatrix

setMethod(".compressFeatureMatrix", signature("mRMRe.Data"), function(object, matrix)
{
    adaptor <- which(object@feature_types != 3)
    matrix <- matrix[adaptor, adaptor]
    colnames(matrix) <- object@feature_names
    rownames(matrix) <- object@feature_names
    
    return(matrix)
})

## expandFeatureIndices

setMethod(".expandFeatureIndices", signature("mRMRe.Data"), function(object, indices)
{
    adaptor <- which(object@feature_types == 3)
    if (length(adaptor) > 0 && any(indices >= adaptor))
        indices <- sapply(indices, function(i) i + sum(sapply(1:length(adaptor), function(j) i >= (adaptor[[j]] - j + 1))))

    return(as.integer(indices))
})

## compressFeatureIndices

setMethod(".compressFeatureIndices", signature("mRMRe.Data"), function(object, indices)
{
    adaptor <- which(object@feature_types == 3)
    
    if (length(adaptor) > 0)
        indices <- sapply(indices, function(i) i - sum(i >= adaptor))
    
    return(as.integer(indices))
})


setMethod("scores", signature("mRMRe.Data"), function(object, solutions)
{
    mi_matrix <- mim(object)
    targets <- names(solutions)
    scores <- lapply(targets, function(target) {
      apply(solutions[[target]], 2, function(solution) {
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
