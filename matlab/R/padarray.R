###
### $Id: padarray.R 55 2014-02-06 16:41:28Z plroebuck $
###
### Pad a vector.
###


##-----------------------------------------------------------------------------
setGeneric("padarray",
           function(A, padsize, padval = 0, direction = c("both",
                                                          "pre",
                                                          "post")) {
               #cat("generic", match.call()[[1]], "\n")
               standardGeneric("padarray")
           })

setMethod("padarray",
          signature(A         = "array",
                    padsize   = "numeric",
                    padval    = "missing",
                    direction = "missing"),
          function(A, padsize, padval, direction) {
              #cat(match.call()[[1]], "(array, numeric, missing, missing)", "\n")
              callGeneric(A, padsize, padval, direction)
          })

setMethod("padarray",
          signature(A         = "array",
                    padsize   = "numeric",
                    padval    = "ANY",
                    direction = "character"),
          function(A, padsize, padval, direction) {
              #cat(match.call()[[1]], "(array, numeric, ANY, character)", "\n")
              method <- "constant"
              padarray0(A, method, padsize, padval, match.arg(direction))
          })

setMethod("padarray",
          signature(A         = "array",
                    padsize   = "numeric",
                    padval    = "character",
                    direction = "character"),
          function(A, padsize, padval, direction) {
              #cat(match.call()[[1]], "(array, numeric, character, character)", "\n")

              ## third parameter is overloaded
              if (padval %in% c("circular",
                                "replicate",
                                "symmetric")) {
                  method <- padval
                  padval <- NA
              } else {
                  method <- "constant"
              }
              padarray0(A, method, padsize, padval, match.arg(direction))
          })

setMethod("padarray",
          signature(A         = "vector",
                    padsize   = "numeric",
                    padval    = "ANY",
                    direction = "ANY"),
          function(A, padsize, padval, direction) {
              #cat(match.call()[[1]], "(vector, numeric, ANY, ANY)", "\n")
              callGeneric(matrix(A, nrow = 1), padsize, padval, direction)
          })


##-----------------------------------------------------------------------------
padarray0 <- function(a,
                      method  = c("constant",
                                  "circular",
                                  "replicate",
                                  "symmetric"),
                      padsize,
                      padval,
                      direction = c("both",
                                    "pre",
                                    "post")) {
    #cat(match.call()[[1]], "(array, character, numeric, ANY, character)", "\n")
    #cat("a         =", a, "\n")
    #cat("method    =", method, "\n")
    #cat("padsize   =", padsize, "\n")
    #cat("padval    =", padval, "\t", "(", data.class(padval), ")", "\n")
    #cat("direction =", direction, "\n")

    method <- match.arg(method)
    direction <- match.arg(direction)

    if (length(padsize) < matlab::ndims(a)) {
        padsize[matlab::ndims(a)] <- 0
    }

    if (!(length(padval) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("padval")))
    }

    if (method == "constant" &&
        !(is.numeric(a) || is.logical(a))) {
        stop(sprintf("argument %s must be numeric or logical for constant padding",
                     sQuote("a")))
    }

    b <- if (matlab::isempty(a)) {
             sizeB <- if (direction == "both") {
                          matlab::size(a) + 2*padsize
                      } else {
                          matlab::size(a) + padsize
                      }
             mkconstarray(data.class(a), padval, sizeB)
         } else {
             switch(EXPR = method,
                    constant  = constantPad(a, padsize, padval, direction),
                    circular  = circularPad(a, padsize, direction),
                    symmetric = symmetricPad(a, padsize, direction),
                    replicate = replicatePad(a, padsize, direction))
         }

     if (is.logical(a)) {
         mode(b) <- "logical"
     }

     return(b)
}


##-----------------------------------------------------------------------------
constantPad <- function(a, padsize, padval, direction) {
    numDims <- matlab::numel(padsize)

    ## Form index vectors to subassign input array into output array.
    ## Also compute the size of the output array.
    idx <- matlab::cell(1, numDims)
    sizeB <- matlab::zeros(1, numDims)
    for (k in seq(1, numDims)) {
        M <- matlab::size(a, k)
        switch(EXPR = direction,
               pre  = {
                          idx[[k]] <- (1:M) + padsize[k]
                          sizeB[k] <- M + padsize[k]
                      },
               post = {
                          idx[[k]] <- 1:M
                          sizeB[k] <- M + padsize[k]
                      },
               both = {
                          idx[[k]] <- (1:M) + padsize[k]
                          sizeB[k] <- M + 2*padsize[k]
                      })
    }

    ## Initialize output array with padding value.
    ## Make sure output array is same type as the input.
    b <- mkconstarray(mode(a), padval, sizeB)

    return(do.call("[<-", c(list(b), idx, list(a))))
}


##-----------------------------------------------------------------------------
circularPad <- function(a, padsize, direction) {
    numDims <- matlab::numel(padsize)

    ## Form index vectors to subassign input array into output array.
    ## Also compute the size of the output array.
    idx <- matlab::cell(1, numDims)
    for (k in seq(1, numDims)) {
        M <- matlab::size(a, k)
        dimNums <- 1:M
        p <- padsize[k]
        switch(EXPR = direction,
               pre  = {
                          idx[[k]] <- dimNums[matlab::mod(-p:(M-1), M) + 1]
                      },
               post = {
                          idx[[k]] <- dimNums[matlab::mod(0:(M+p-1), M) + 1]
                      },
               both = {
                          idx[[k]] <- dimNums[matlab::mod(-p:(M+p-1), M) + 1]
                      })
    }

    return(do.call("[", c(list(a), idx)))
}


##-----------------------------------------------------------------------------
symmetricPad <- function(a, padsize, direction) {
    numDims <- matlab::numel(padsize)

    ## Form index vectors to subassign input array into output array.
    ## Also compute the size of the output array.
    idx <- matlab::cell(1, numDims)
    for (k in seq(1, numDims)) {
        M <- matlab::size(a, k)
        dimNums <- c(1:M, seq(from = M, to = 1, by = -1))
        p <- padsize[k]
        switch(EXPR = direction,
               pre  = {
                          idx[[k]] <- dimNums[matlab::mod(-p:(M-1), 2*M) + 1]
                      },
               post = {
                          idx[[k]] <- dimNums[matlab::mod(0:(M+p-1), 2*M) + 1]
                      },
               both = {
                          idx[[k]] <- dimNums[matlab::mod(-p:(M+p-1), 2*M) + 1]
                      })
    }

    return(do.call("[", c(list(a), idx)))
}


##-----------------------------------------------------------------------------
replicatePad <- function(a, padsize, direction) {
    numDims <- matlab::numel(padsize)

    ## Form index vectors to subassign input array into output array.
    ## Also compute the size of the output array.
    idx <- matlab::cell(1, numDims)
    for (k in seq(1, numDims)) {
        M <- matlab::size(a, k)
        p <- padsize[k]
        onesVector <- if (p > 0) {
                          matlab::ones(1, p)
                      } else {
                          NULL
                      }
        switch(EXPR = direction,
               pre  = {
                          idx[[k]] <- c(onesVector, 1:M)
                      },
               post = {
                          idx[[k]] <- c(1:M, M*onesVector)
                      },
               both = {
                          idx[[k]] <- c(onesVector, 1:M, M*onesVector)
                      })
    }

    return(do.call("[", c(list(a), idx)))
}

