#' @include distancematrix.R
NULL

#'Nonbipartite Matching
#'
#'The nonbinmatch function creates the set of pairwise matches that minimizes
#'the sum of distances between the pairs.
#'
#'The nonbinmatch function calls the Fortran code (Derigs) and set of pairwise
#'matches that minimizes the sum of distances between the pairs.
#'
#'@aliases nonbimatch nonbimatch,distancematrix-method nonbimatch-class
#'@param mdm A distancematrix object.  See the distancematrix function.
#'@param threshold An numeric value, indicating the distance needed to create
#'chameleon matches.
#'@param precision The largest value in the matrix will have at most this many
#'digits.  The default value is six.
#'@param \dots Additional arguments, these are not used.
#'@return nonbimatch S4 object with several elements
#'
#'  \item{matches}{data.frame containing matches}
#'
#'  \item{halves}{data.frame containing each match}
#'
#'  \item{total}{sum of the distances across all pairs}
#'
#'  \item{mean}{mean distance for each pair}
#'@exportClass nonbimatch
#'@exportMethod nonbimatch
#'@author Cole Beck
#'@seealso \code{\link{distancematrix}}
#'@examples
#'
#'plainmatrix<-as.matrix(dist(sample(1:25, 8, replace=TRUE)))
#'diag(plainmatrix) <- 99999  # setting diagonal to an infinite distance for
#'                            # pedagogical reasons (the diagonal may be left
#'                            # as zero)
#'mdm<-distancematrix(plainmatrix)
#'res<-nonbimatch(mdm)
#'

setGeneric("nonbimatch", function(mdm, threshold=NA, precision=6, ...) standardGeneric("nonbimatch"))
setMethod("nonbimatch", "distancematrix", function(mdm, threshold=NA, precision, ...) {
    if(any(is.na(as.numeric(mdm)))) {
        stop("Elements of a distance matrix must be numeric")
    }
    n <- nrow(mdm)
    if(n == 0) {
        stop("distance matrix has no rows")
    }
    if(n%%2 == 1) {
        stop("There must be an even number of elements")
    }
    if(is.null(rownames(mdm))) {
        rownames(mdm) <- seq_len(n)
    }
    ids <- rownames(mdm)
    if(!is.na(threshold) && threshold >= 0) {
        # extend the distance matrix for chameleons
        newvals <- rep(threshold, n)
        # node stack overflow may occur if distancematrix
        suppressWarnings(class(mdm) <- 'matrix')
        # add columns
        mdm <- do.call("cbind", c(list(mdm), newvals))
        # add rows
        mdm <- do.call("rbind", c(list(mdm), newvals))
        n <- n*2
    } else threshold <- NA

    wt <- c(mdm)
    nmatch <- seq_len(n)
    if(!is.numeric(precision) || precision < 1) {
        precision <- 6
        warning("Precision value is too small.  Setting precision to six.")
    } else if(precision >= 10) {
        precision <- 6
        warning("Precision value is too large.  Setting precision to six.")
    }
    myinf <- is.infinite(wt)
    numdigits <- max(0, floor(log10(max(wt[!myinf])))) + 1
    shift <- 10^(precision-numdigits)
    if(shift != 1) wt <- wt*shift
    # replace Inf with new max
    if(any(myinf)) {
        wt[myinf] <- 2*10^precision
        # update precision based on new max value
        precision <- precision + 1
    }
    # the real test should be:
    # max(wt) > .Machine$integer.max
    # the largest number will have at most [precision] digits (defaulting to six)
    # warning: if the vector is large and there are too many digits, the Fortran call will crash
    if(shift < 1) {
        print(sprintf("Note: Distances scaled by %s to ensure all data can be handled", shift))
    }
    match <- .Fortran(mwrap, n=as.integer(n), wt=as.integer(wt), nmatch=as.vector(nmatch), prcn=as.integer(precision))$nmatch

    # remove chameleon to chameleon matches
    if(!is.na(threshold)) {
        c2c <- which(match[seq(from=(n/2)+1, to=n)] > n/2)
        # number of chameleons that match elements
        ncham <- n/2 - length(c2c)
        ids <- append(ids, sprintf('chameleon%s', seq_len(ncham)))
        if(length(c2c) > 0L) {
            match <- match[-(c2c + n/2)]
            index.to.replace <- which(match > n/2)
            # re-label chameleon matches from 1 to ncham
            if(length(index.to.replace) > 0L) {
                match[index.to.replace] <- seq(from=(n/2+1), to=(n/2+ncham))
                match[seq(from=(n/2+1), to=(n/2+ncham))] <- index.to.replace
            }
            n <- length(match)
        }
    }

    i <- seq_len(n)
    matches <- data.frame("Group1.ID"=ids, "Group1.Row"=i, "Group2.ID"=ids[match], "Group2.Row"=match, "Distance"=mdm[cbind(i, match)])
    halves <- matches[matches[,2] < matches[,4],]
    distance <- sum(matches[,5])
    total <- distance/2
    mean <- distance/n
    new("nonbimatch", list(matches=matches, halves=halves, total=total, mean=mean))
})

setClass("nonbimatch", contains = "list")

setMethod("show", signature(object="nonbimatch"), function(object) {
    print(setNames(object@.Data, names(object)))
})

# runner starts with covariate matrix, generates distances, finds matches, and reports QoM
setGeneric("runner", function(covariate, seed=101, ..., mate.random=FALSE) standardGeneric("runner"))
setMethod("runner", "data.frame", function(covariate, seed=101, ..., mate.random=FALSE) {
    step1 <- gendistance(covariate, ...)
    step2 <- distancematrix(step1, ...)
    if(mate.random) {
        if(exists(".Random.seed", envir = .GlobalEnv)) {
            save.seed <- get(".Random.seed", envir= .GlobalEnv)
            on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
        } else {
            on.exit(rm(".Random.seed", envir = .GlobalEnv))
        }
        if(!is.numeric(seed)) seed <- 101
        set.seed(seed)
        n <- nrow(step2)
        options <- seq_len(n)
        result <- numeric(n)
        for(i in seq_len(n/2)) {
            mates <- sample(options, 2)
            result[mates[1]] <- mates[2]
            result[mates[2]] <- mates[1]
            options <- options[-match(mates, options)]
        }
        i <- seq_len(n)
        ids <- rownames(step2)
        matches <- data.frame("Group1.ID"=ids, "Group1.Row"=i, "Group2.ID"=ids[result], "Group2.Row"=result, "Distance"=step2[i,result])
        halves <- matches[matches[,2] < matches[,4],]
        distance <- sum(matches[,5])
        total <- distance/2
        mean <- distance/n
        step3 <- list(matches=matches, halves=halves, total=total, mean=mean)
    } else {
        step3 <- nonbimatch(step2, ...)
    }
    step4 <- qom(step1$cov, step3$matches, ...)
    step5 <- assign.grp(step3$matches, ...)
    step6 <- cbind(step1$cov, step5[seq_len(nrow(step1$cov)),c(3,4,6)])
    return(list(setup=step1, mdm=step2, matches=step3, qom=step4, grps=step5, final=step6))
})

setGeneric("full.qom", function(covariate, iterations=NA, ...) standardGeneric("full.qom"))
setMethod("full.qom", "data.frame", function(covariate, iterations=NA, ...) {
    qom <- vector('list', ncol(covariate)+3)
    if(is.na(iterations) || !is.numeric(iterations) || iterations < 2) iterations <- 10000
    a <- runner(covariate, iterations=iterations, ...)
    ids <- rownames(a$setup$cov)
    cnames <- colnames(covariate)
    qom[[1]] <- a$qom
    qom[[2]] <- runner(covariate, iterations=iterations, mate.random=TRUE)$qom
    qom[[3]] <- runner(covariate, iterations=iterations, missingness=0)$qom
    for(i in seq_len(ncol(covariate))) {
        weights <- rep(0, ncol(covariate))
        weights[i] <- 1
        qom[[3+i]] <- tryCatch(runner(covariate, iterations=iterations, weights=weights, missingness=0)$qom, error=function(e) {})
    }
    names(qom) <- c("user.specified", "random.mates", "eq.weight", cnames)
    qom
})
