setGeneric("featureData", function(object) standardGeneric("featureData"))

setGeneric("subsetData", function(object, ...) standardGeneric("subsetData"))

setGeneric("sampleNames", function(object) standardGeneric("sampleNames"))

setGeneric("sampleCount", function(object) standardGeneric("sampleCount"))

setGeneric("featureCount", function(object) standardGeneric("featureCount"))

setGeneric("featureNames", function(object) standardGeneric("featureNames"))

setGeneric("sampleStrata", function(object) standardGeneric("sampleStrata"))

setGeneric("sampleStrata<-", function(object, value) standardGeneric("sampleStrata<-"))

setGeneric("sampleWeights", function(object) standardGeneric("sampleWeights"))

setGeneric("sampleWeights<-", function(object, value) standardGeneric("sampleWeights<-"))

setGeneric("priors", function(object) standardGeneric("priors"))

setGeneric("priors<-", function(object, value) standardGeneric("priors<-"))

setGeneric("mim", function(object, method = c("mi", "cor"), ...)
{
    method <- match.arg(method)
    matrix <- standardGeneric("mim")
    
    if (method == "mi")
        matrix <- -.5 * log(1 - (matrix^2))
    
    return(matrix)
})

setGeneric(".expandFeatureMatrix", function(object, ...) standardGeneric(".expandFeatureMatrix"))

setGeneric(".compressFeatureMatrix", function(object, ...) standardGeneric(".compressFeatureMatrix"))

setGeneric(".expandFeatureIndices", function(object, ...) standardGeneric(".expandFeatureIndices"))

setGeneric(".compressFeatureIndices", function(object, ...) standardGeneric(".compressFeatureIndices"))

setGeneric("solutions", function(object, ...) standardGeneric("solutions"))

setGeneric("scores", function(object, ...) standardGeneric("scores"))

setGeneric("causality", function(object, ...) standardGeneric("causality"))

setGeneric("target", function(object) standardGeneric("target"))

setGeneric("adjacencyMatrix", function(object) standardGeneric("adjacencyMatrix"))

setGeneric("adjacencyMatrixSum", function(object) standardGeneric("adjacencyMatrixSum"))

setGeneric("visualize", function(object) standardGeneric("visualize"))

`.map.continuous.estimator` <- function(continuous_estimator)
{
    value <- switch(continuous_estimator, "pearson" = 0L, "spearman" = 1L, "kendall" = 2L, "frequency" = 3L, -1L)
    
    if (value < 0L || value > 4L || !is.character(continuous_estimator))
        stop("estimator must be of the following: pearson, spearman, kendall, frequency")
    
    return(value)
}

`correlate` <- function(X, Y, method = c("pearson", "spearman", "kendall", "frequency", "cramersv", "cindex"), strata, weights, outX = TRUE, bootstrap_count = 0, alpha = 0.05, alternative=c("two.sided", "less", "greater"))
{
    method <- match.arg(method)
    alternative <- match.arg(alternative)
    
    if((is.Surv(X) || is.Surv(Y)) && method != "cindex") { stop("method should be cindex when dealing with survival data") }
    
    if (method == "pearson" || method == "spearman" || method == "kendall" || method == "frequency")
    {
        X <- as.numeric(X)
        Y <- as.numeric(Y)
    }
    else if (method == "cramersv")
    {
        X <- as.factor(X)
        Y <- as.factor(Y)
    }
    else if (method != "cindex")
        stop("estimator must be of the following: pearson, spearman, kendall, frequency, cramersv, cindex")
    
    if(is.Surv(X)) { ll <- nrow(X) } else { ll <- length(X) }

    if (missing(strata)) {
      strata <- factor(rep(0, ll))
      names(strata) <- names(X)
    }
    
    if (missing(weights)) {
      weights <- rep(1, ll)
      names(weights) <- names(X)
    } 
    
    data <- mRMR.data(data = data.frame(X, Y), strata = strata, weights = weights)

    if (method == "cindex")
    {
        empty <- vector(mode = "numeric", length = 0)
        
        if (length(data@feature_types) == 2)
            input <- list(data@data[, 1], data@data[, 2], empty, empty)
        else if (length(data@feature_types) == 3)
        {
            if (data@feature_types[[1]] == 2)
                input <- list(data@data[, 1], data@data[, 3], data@data[, 2], empty)
            else if (data@feature_types[[2]] == 2)
                input <- list(data@data[, 2], data@data[, 1], data@data[, 3], empty)
        }
        else if (length(data@feature_types) == 4)
            input <- list(data@data[, 1], data@data[, 3], data@data[, 2], data@data[, 4])
        
				ratio <- vector(mode = "numeric", length = 1)
				ch <- vector(mode = "numeric", length = length(input[[1]]))
				dh <- vector(mode = "numeric", length = length(input[[1]]))
				uh <- vector(mode = "numeric", length = length(input[[1]]))
				rh <- vector(mode = "numeric", length = length(input[[1]]))

        .Call(.C_export_concordance_index, as.numeric(input[[1]]), as.numeric(input[[2]]),
                as.numeric(input[[3]]), as.numeric(input[[4]]), as.integer(data@strata), as.numeric(data@weights),
                as.integer(length(unique(data@strata))), outX, ratio, ch, dh, uh, rh)  
        
              cindex <- ratio
              myx <- complete.cases(featureData(data), sampleStrata(data), sampleWeights(data))
              N <- sum(weights[myx])
              
              cscount <- sum(ch + dh) ## comparable pairs
              if (sum(ch) == 0 || sum (dh) ==0 || sum(ch * (ch - 1)) == 0 || sum(dh * (dh - 1)) == 0 || sum(ch * dh) == 0 || cscount < 10)
                return(list("cindex"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p"=NA, "n"=N))
              
              ## FIXME: N and subsequent calculation should be done withing strata
              pc <- (1 / (N * (N - 1))) * sum(ch)
              pd  <- (1 / (N * (N - 1))) * sum(dh)
              pcc <- (1 / (N * (N - 1) * (N - 2))) * sum(ch * (ch - 1))
              pdd <- (1 / (N * (N - 1) * (N - 2))) * sum(dh * (dh - 1))
              pcd <- (1 / (N * (N - 1) * (N - 2))) * sum(ch * dh)
              varp <- (4 / (pc + pd)^4) * (pd^2 * pcc - 2 * pc * pd * pcd + pc^2 * pdd)
              if((varp / N) > 0) {
                se <- sqrt(varp / N)
                ci <- qnorm(p=alpha / 2, lower.tail=FALSE) * se
                lower <- cindex - ci
                upper <- cindex + ci
                switch(alternative, 
                "two.sided"={ p <- pnorm((cindex - 0.5) / se, lower.tail=cindex < 0.5) * 2 }, 
                "less"={ p <- pnorm((cindex - 0.5) / se, lower.tail=TRUE) }, 
                "greater"={  p <- pnorm((cindex - 0.5) / se, lower.tail=FALSE) }
                )
              } else { se <- lower <- upper <- p <- NA } 
        
        return(list("estimate"=cindex, "se"=se, "lower"=lower, "upper"=upper, "p"=p, "n"=N))
    }
    else if (method == "cramersv")
        return(list(statistic = mim(data, method = "cor", outX = outX, bootstrap_count = bootstrap_count)[1, 2]))
    else
        return(list(statistic = mim(data, method = "cor", continuous_estimator = method, outX = outX,
                            bootstrap_count = bootstrap_count)[1, 2]))
}

`get.thread.count` <- function()
{
    thread_count <- vector(mode = "integer", length = 1)
    
    .Call(.C_get_thread_count, thread_count)
    
    return(thread_count)
}

`set.thread.count` <- function(thread_count)
{
    thread_count <- as.integer(thread_count)
    
    .Call(.C_set_thread_count, thread_count)
    
    return(thread_count)
}
