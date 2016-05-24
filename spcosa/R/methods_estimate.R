setMethod(
    f = "estimate",
    signature = signature(
        statistic = "character",
        stratification = "CompactStratificationEqualArea",
        samplingPattern = "SamplingPatternRandomComposite",
        data = "data.frame"
    ),
    definition = function(statistic, stratification, samplingPattern, data, ...) {
        
        # check statistic
        statistic <- match.arg(
            arg = statistic,
            choices = c("spatial mean", "sampling variance", "standard error"),
            several.ok = FALSE
        )
        
        # delegation
        switch(
            statistic,
            "spatial mean"      = estimate(new("SpatialMean"), stratification, samplingPattern, data, ...),
            "sampling variance" = estimate(new("SamplingVariance"), stratification, samplingPattern, data, ...),
            "standard error"    = estimate(new("StandardError"), stratification, samplingPattern, data, ...)
        )
    }
)



setMethod(
    f = "estimate",
    signature = signature(
        statistic = "character",
        stratification = "CompactStratification",
        samplingPattern = "SamplingPatternRandomSamplingUnits",
        data = "data.frame"
    ),
    definition = function(statistic, stratification, samplingPattern, data, ...) {
        
        # check statistic
        statistic <- match.arg(
            arg = statistic,
            choices = c("spatial mean", "spatial variance", "sampling variance", "standard error", "scdf"),
            several.ok = FALSE
        )
        
        # delegation
        switch(
            statistic,
            "spatial mean"      = estimate(new("SpatialMean"), stratification, samplingPattern, data, ...),
            "sampling variance" = estimate(new("SamplingVariance"), stratification, samplingPattern, data, ...),
            "standard error"    = estimate(new("StandardError"), stratification, samplingPattern, data, ...),
            "spatial variance"  = estimate(new("SpatialVariance"), stratification, samplingPattern, data, ...),
            "scdf"              = estimate(new("SpatialCumulativeDistributionFunction"), stratification,
                                           samplingPattern, data, ...)
        )
    }
)



setMethod(
    f = "estimate",
    signature = signature(
        statistic = "SamplingVariance",
        stratification = "CompactStratificationEqualArea",
        samplingPattern = "SamplingPatternRandomComposite",
        data = "data.frame"
    ),
    definition = function(statistic, stratification, samplingPattern, data, ...) {
        
        # check if data is available for each sampling location
        sampleSize <- getSampleSize(samplingPattern)
        if (nrow(data) != sampleSize) {
            stop("number of data should be equal to the number of composites,", call. = FALSE)
        }
        
        # estimate the sampling variance (eq. 7.5, De Gruijter et. al, 2006)
        apply(X = data, MARGIN = 2, FUN = var) / sampleSize
    }
)



setMethod(
    f = "estimate",
    signature = signature(
        statistic = "SamplingVariance",
        stratification = "CompactStratification",
        samplingPattern = "SamplingPatternRandomSamplingUnits",
        data = "data.frame"
    ),
    definition = function(statistic, stratification, samplingPattern, data, ...) {
        
        # check if data is available for each sampling location
        sampleSize <- getSampleSize(samplingPattern)
        if (nrow(data) != sampleSize) {
            stop("number of data should be equal to the number of sampling locations,", call. = FALSE)
        }
        
        # get relative area 'a_h' of each stratum 'h'
        a_h <- getRelativeArea(stratification)
        
        # retrieve stratum id 'h' for each sampling point
        H <- getNumberOfStrata(stratification)
        n <- getSampleSize(samplingPattern)
        n_h <- n / H
        h <- rep(x = seq_len(H), each = n_h)
        
        # compute the sample variance for each stratum h
        tmp <- lapply(X = data, FUN = function(x) {
            tapply(X = x, INDEX = h, FUN = var)})
        var_z_h <- as.matrix(as.data.frame(tmp)) / n_h
        
        # compute the sampling variance var_z
        drop(crossprod(a_h * a_h, var_z_h))
    }
)



setMethod(
    f = "estimate",
    signature = signature(
        statistic = "SpatialCumulativeDistributionFunction",
        stratification = "CompactStratification",
        samplingPattern = "SamplingPatternRandomSamplingUnits",
        data = "data.frame"
    ),
    definition = function(statistic, stratification, samplingPattern, data, ...) {
        lapply(
            X = data,
            FUN = function(y) {
                ys <- sort(unique(y))
                cbind(
                    value = ys,
                    cumFreq = as.vector(
                        sapply(
                            X = ys,
                            FUN = function(threshold) {
                                estimate(
                                    statistic = new("SpatialMean"),
                                    stratification = stratification,
                                    samplingPattern = samplingPattern,
                                    data = data.frame(
                                        i = ifelse(y < threshold, 1, 0)
                                    )
                                )
                            }
                        )
                    )
                )
            }
        )
    }
)



setMethod(
    f = "estimate",
    signature = signature(
        statistic = "SpatialMean",
        stratification = "CompactStratificationEqualArea",
        samplingPattern = "SamplingPatternRandomComposite",
        data = "data.frame"
    ),
    definition = function(statistic, stratification, samplingPattern, data, ...) {
        
        # check if data is available for each sampling location
        sampleSize <- getSampleSize(samplingPattern)
        if (nrow(data) != sampleSize) {
            stop("number of data should be equal to the number of composites,", call. = FALSE)
        }
        
        # estimate the spatial mean (eq. 7.4, De Gruijter et. al, 2006)
        colMeans(data)
    }
)



setMethod(
    f = "estimate",
    signature = signature(
        statistic = "SpatialMean",
        stratification = "CompactStratification",
        samplingPattern = "SamplingPatternRandomSamplingUnits",
        data = "data.frame"
    ),
    definition = function(statistic, stratification, samplingPattern, data, ...) {
        
        # check if data is available for each sampling location
        sampleSize <- getSampleSize(samplingPattern)
        if (nrow(data) != sampleSize) {
            stop("number of data should be equal to the number of sampling locations,", call. = FALSE)
        }
        
        # get relative area 'a_h' of each stratum 'h'
        a_h <- getRelativeArea(stratification)
        
        # retrieve stratum id 'h' for each sampling point
        H <- getNumberOfStrata(stratification)
        n <- getSampleSize(samplingPattern)
        n_h <- n / H
        h <- rep(x = seq_len(H), each = n_h)
        
        # compute the sample mean for each stratum h
        tmp <- lapply(X = data, FUN = function(x) {
            tapply(X = x, INDEX = h, FUN = mean)})
        mean_z_h <- as.matrix(as.data.frame(tmp))
        
        # compute the spatial mean mean_z
        drop(crossprod(a_h, mean_z_h))
    }
)



setMethod(
    f = "estimate",
    signature = signature(
        statistic = "SpatialVariance",
        stratification = "CompactStratification",
        samplingPattern = "SamplingPatternRandomSamplingUnits",
        data = "data.frame"
    ),
    definition = function(statistic, stratification, samplingPattern, data, ...) {
        mean_z <- estimate("spatial mean", stratification, samplingPattern, data, ...)
        mean_z2 <- estimate("spatial mean", stratification, samplingPattern, data * data, ...) 
        var_z <- estimate("sampling variance", stratification, samplingPattern, data, ...) 
        mean_z2 - (mean_z)^2 + var_z
    }
)



setMethod(
    f = "estimate",
    signature = signature(
        statistic = "StandardError",
        stratification = "CompactStratification",
        samplingPattern = "SamplingPatternRandomSamplingUnits",
        data = "data.frame"
    ),
    definition = function(statistic, stratification, samplingPattern, data, ...) {
        sqrt(callNextMethod())
    }
)
