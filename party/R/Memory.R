
# $Id: Memory.R 532 2014-02-01 08:56:55Z thothorn $

ctree_memory <- function(object, MPinv = FALSE) {

    inputs <- object@inputs
    responses <- object@responses
    q <- ncol(responses@test_trafo)
    nobs <- inputs@nobs
    ninputs <- inputs@ninputs
    RET <- new("TreeFitMemory")
    expcovinf <- new("ExpectCovarInfluence", q)
    expcovinfss <- new("ExpectCovarInfluence", as.integer(1))
    linexpcov2sample  <- new("LinStatExpectCovar", as.integer(1), q) 
    splitstatistics <- as.double(rep(0, nobs))
    varmemory <- vector(mode = "list", length = ninputs)
    dontuse <- rep(FALSE, ninputs)
    dontusetmp <- rep(FALSE, ninputs)

    for (j in 1:inputs@ninputs) {

        p <- ncol(inputs@transformations[[j]])

        if (MPinv) {
            varmemory[[j]] <- new("LinStatExpectCovarMPinv", p, q)
        } else {
            varmemory[[j]] <- new("LinStatExpectCovar", p, q)
        }
    }

    RET@expcovinf <- expcovinf
    RET@expcovinfss <- expcovinfss
    RET@linexpcov2sample  <- linexpcov2sample
    RET@weights <- as.double(rep(0, nobs))
    RET@varmemory <- varmemory
    RET@dontuse <- as.logical(dontuse)
    RET@dontusetmp <- as.logical(dontusetmp)
    RET@splitstatistics <- splitstatistics
    RET
}
