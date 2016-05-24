##--------------------------------------------------------------------##
##                                                                    ##
## Split sample a test set along gradient of interest                 ##
##                                                                    ##
## split gradient into chunkk sections and sample randomly nc         ##
## observations from within each section to act as the test sample    ##
##                                                                    ##
## env    - gradient to sample along                                  ##
## chunk  - number of chunks to split env into                        ##
## n      - number of observations to sample in total                 ##
## nchunk - number of observations per group                          ##
##                                                                    ##
##--------------------------------------------------------------------##
splitSample <- function(env, chunk = 10, take, nchunk,
                        fill = c("head","tail","random"),
                        maxit = 1000) {
    sampFun <- function(ind, x, nchunk) {
        sample(x[[ind]], min(length(x[[ind]]), nchunk[[ind]]))
    }
    if(take < chunk) {
        stop("Number of samples to 'take' < number of 'chunk's")
    }
    fill <- match.arg(fill)
    env <- as.numeric(env)
    N <- length(env)
    cuts <- cut(env, chunk)
    splt <- split(seq_along(env), cuts)
    lens <- sapply(splt, length)
    if(missing(nchunk)) {
        ## base number of samples per chunk
        each <- floor(take / chunk)
        ## expand to a vector
        nchunk <- rep(each, chunk)
        ## fill in the remainder samples according to fill type
        tooSmall <- lens < nchunk
        nchunk[tooSmall] <- lens[tooSmall]
        if(sum(nchunk) < take) {
            i <- iter <- 1
            ## vector of chunks that *aren't* too small expanded to length 100
            vec <- if(isTRUE(all.equal(fill, "head"))) {
                rep(sort(which(!tooSmall)), 100)
            } else if(isTRUE(all.equal(fill, "tail"))) {
                rep(sort(which(!tooSmall), decreasing = TRUE), 100)
            } else {
                sample(rep(which(!tooSmall), 100), 100)
            }
            ## fill in chunks
            while(sum(nchunk) < take) {
                want <- vec[i]
                i <- i + 1
                iter <- iter + 1
                if(iter == maxit) {
                    warning("Failed to allocate all 'take' samples in 'maxit' iterations.")
                    break
                }
                if(lens[want] <= nchunk[want])
                    next
                nchunk[want] <- nchunk[want] + 1
                if(i > 100) ## if used all vec, start again
                    i <- 1
            }
        }
    } else {
        if(!is.numeric(env))
            stop("'env' must be a numeric vector.")
        if(!isTRUE(all.equal(sum(nchunk), length(env))))
            stop("'length(env)' and 'sum(nchunk)' do not match.")
    }
    ## indicator to loop over
    ind <- seq_along(splt)
    ## loop over gradient chunks
    samp <- lapply(ind, FUN = sampFun, x = splt, nchunk = nchunk)
    ## grab the number of samples in each chunk
    lengths <- sapply(samp, length)
    ## turn sample list into a vector
    samp <- unlist(samp, use.names = FALSE)
    ## assign lengths as attribute
    attr(samp, "lengths") <- lengths
    samp
}
