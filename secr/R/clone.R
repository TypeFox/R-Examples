############################################################################################
## package 'secr'
## clone.R
## last changed 2014-12-12
############################################################################################

clone <- function (object, type, ...) UseMethod("clone")

clone.default <- function (object,  type, ...)       {
    if (length(dim(object)) != 2)
        stop ("requires 2-D object")
    type <- tolower (type)
    n <- nrow(object)
    if (n == 0) {
        out <- object
    }
    else {
        if (type == 'constant')
            freq <- rep(..., n)
        else if (type == 'poisson')
            freq <- rpois(n, ...)
        else if (type == 'nbinom')
            freq <- rnbinom(n, ...)
        else
            stop("unrecognised type")
        index <- rep(1:n, freq)
        object[index,]
    }
}

clone.popn <- function (object, type, ...) {
    if (ms(object)) {
        out <- lapply (object, clone, type, ...)
        class (out) <- c('popn','list')
        out
    }
    else {
        type <- tolower (type)
	n <- nrow(object)
        if (n == 0) {
            out <- object
        }
        else {
            if (type == 'constant')
                freq <- rep(..., n)
            else if (type == 'poisson')
                freq <- rpois(n, ...)
            else if (type == 'nbinom')
                freq <- rnbinom(n, ...)
            else
                stop("unrecognised type")
            index <- rep(1:n, freq)
            out <- object[index,]
            if (!is.null(covariates(object)))
                covariates(out) <- covariates(object)[index,]
            attr (out, 'freq') <- freq
        }
        out
    }
}

clone.capthist <- function (object, type, ...) {
    if (ms(object)) {
        out <- lapply (object, clone, type, ...)
        class (out) <- class(object)
        out
    }
    else {
        type <- tolower (type)
	n <- nrow(object)
        if (n == 0) {
            out <- object
        }
        else {
            if (type == 'constant')
                freq <- rep(..., n)
            else if (type == 'poisson')
                freq <- rpois(n, ...)
            else if (type == 'nbinom')
                freq <- rnbinom(n, ...)
            else
                stop("unrecognised type")
            index <- rep(1:n, freq)
            if (length(dim(object))==2)
                out <- object[index,]
            else
                out <- object[index,,]

            seqn <- unlist(lapply(freq[freq>0], seq, from = 1))
            seqn <- leadingzero(seqn)
            rown <- paste(rep(rownames(object), freq), seqn, sep='.')
            rownames(out) <- rown
            traps(out) <- traps(object)
            attr(out, 'cutval') <- attr(object, 'cutval')
            if (!is.null(covariates(object))) {
                covariates(out) <- covariates(object)[index,]
                rownames(covariates(out)) <- rown
            }

            if (!is.null(xy(object)) | !is.null(signalframe(object))) {
                ## these attributes are defined for each detection
                ## so we first extract the numeric animal ID for each detection
                detectionindex <- animalID(object, names = FALSE)
                xy(out) <- xy(object)[index[detectionindex]]
                signalframe(out) <- signalframe(object)[index[detectionindex]]
            }
            attr (out, 'freq') <- freq
            class(out) <- class(object)
        }

        out
    }
}
