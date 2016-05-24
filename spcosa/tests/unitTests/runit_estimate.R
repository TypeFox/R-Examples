# Unit test: "estimate"
# Tests the estimate method against independent
# (but less efficient) alternatives
test_estimate <-
function() {

    # load sp
    if (suppressWarnings(!require(sp))) {
        stop("This unit test requires package 'sp'.\nThis package is currently not available. Please install 'sp' first.", call. = FALSE)
    }    
    
    # construct map
    map <- expand.grid(s1 = 1:15, s2 = 1:15)
    coordinates(map) <- ~ s1 * s2
    gridded(map) <- TRUE

    # initialize pseudo random number generator
    set.seed(700124)

    # loop over number of strata
    for (k in c(1:10, 1:5 * 20)) {

        # loop over number of samples per stratum
        for (n in 1:10) {

            # sample map
            myStratification <- stratify(map, nStrata = k)
            mySamplingPattern <- spsample(myStratification, n = n)
    
            # construct data.frame
            sampleSize <- getSampleSize(mySamplingPattern)
            myData <- data.frame(obs = runif(n = sampleSize))
    
            # Estimate statistics by ignoring implemented convenience methods
            # (getters like getRelativeArea, etc). In this way, the estimated
            # statistics are tested (almost) independently from the actual
            # implementation.
            #
            # notation:
            # h: stratum
            # nh: number of strata
            # ah: relative area of stratum h
            # zh: mean of stratum h
            # varzh: variance of zh
            # z: spatial mean
            # varz: sampling variance
            # S: spatial variance
            cellId <- as(mySamplingPattern, "SpatialPoints") %over% map       
            h <- myStratification@stratumId[cellId]
            nh <- nrow(myData) / length(unique(h))
            ah <- table(myStratification@stratumId) / length(myStratification@stratumId)
            zh <- tapply(X = myData$obs, INDEX = h, FUN = mean)
            z2h <- tapply(X = myData$obs^2, INDEX = h, FUN = mean)
            z <- sum(ah * zh)
            z2 <- sum(ah * z2h)
            varzh <- tapply(X = myData$obs, INDEX = h, FUN = var) / nh
            varz <- sum(ah^2 * varzh)
            S <- z2 - z^2 + varz
    
            # check spatial mean (and partly, the SCDF)
            checkEqualsNumeric(
                target  = z,
                current = estimate("spatial mean", myStratification, mySamplingPattern, myData)
            )
    
            # check standard error
            checkEqualsNumeric(
                target  = sqrt(varz),
                current = estimate("standard error", myStratification, mySamplingPattern, myData)
            )
        
            # check sampling variance
            checkEqualsNumeric(
                target  = varz,
                current = estimate("sampling variance", myStratification, mySamplingPattern, myData)
            )
        
            # check spatial variance
            checkEqualsNumeric(
                target  = S,
                current = estimate("spatial variance", myStratification, mySamplingPattern, myData)
            )
        }
    }
}
