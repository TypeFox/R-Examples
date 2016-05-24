
partysplit <- function(varid, breaks = NULL, index = NULL, right = TRUE, 
                  prob = NULL, info = NULL) {

    ### informal class for splits
    split <- vector(mode = "list", length = 6)
    names(split) <- c("varid", "breaks", "index", "right", "prob", "info")

    ### split is an id referring to a variable
    stopifnot(is.integer(varid))
    split$varid <- varid

    if (is.null(breaks) && is.null(index))
        stop("either", " ", sQuote("breaks"), " ", "or", " ",
             sQuote("index"), " ", "must be given")

    ### vec
    if (!is.null(breaks)) {
        if (is.numeric(breaks) && (length(breaks) >= 1)) {
            ### FIXME: I think we need to make sure breaks are double in C
            split$breaks <- as.double(breaks)
        } else {
            stop(sQuote("break"), " ",
                 "should be a numeric vector containing at least one element")
        }
    }

    if (!is.null(index)) {
        if (is.integer(index)) {
            if (!(length(index) >= 2)) 
                stop(sQuote("index"), " ", "has less than two elements")
            if (!(min(index, na.rm = TRUE) == 1))
                stop("minimum of", " ", sQuote("index"), " ", "is not equal to 1")
            if (!all.equal(diff(sort(unique(index))), rep(1, max(index, na.rm = TRUE) - 1)))
                stop(sQuote("index"), " ", "is not a contiguous sequence")
            split$index <- index
        } else {
            stop(sQuote("index"), " ", "is not a class", " ", sQuote("integer"))
        }
        if (!is.null(breaks)) {
            if (length(breaks) != (length(index) - 1))
                stop("length of", " ", sQuote("breaks"), " ", 
                     "does not match length of", " ", sQuote("index"))
        }
    }

    if (is.logical(right) & !is.na(right))
        split$right <- right
    else
        stop(sQuote("right"), " ", "is not a logical")

    if (!is.null(prob)) {
        if (!is.double(prob) || 
            (any(prob < 0) | any(prob > 1) | !isTRUE(all.equal(sum(prob), 1))))
            stop(sQuote("prob"), " ", "is not a vector of probabilities")
        if (!is.null(index))
            stopifnot(max(index, na.rm = TRUE) == length(prob))
        if (!is.null(breaks) && is.null(index))
            stopifnot(length(breaks) == (length(prob) - 1))
        split$prob <- prob
    }

    if (!is.null(info))
        split$info <- info

    class(split) <- "partysplit"

    return(split)
}

varid_split <- function(split) {
    stopifnot(inherits(split, "partysplit"))
    split$varid
}

breaks_split <- function(split) {
    stopifnot(inherits(split, "partysplit"))
    split$breaks
}

index_split <- function(split) {
    stopifnot(inherits(split, "partysplit"))
    split$index
}

right_split <- function(split) {
    stopifnot(inherits(split, "partysplit"))
    split$right
}

prob_split <- function(split) {
    stopifnot(inherits(split, "partysplit"))
    prob <- split$prob
    if (!is.null(prob)) return(prob)

    ### either breaks or index must be there
    if (is.null(index <- index_split(split))) {
        if (is.null(breaks <- breaks_split(split)))
            stop("neither", " ", sQuote("prob"), " ", "nor", " ", 
                 sQuote("index"), " ", "or", sQuote("breaks"), " ", 
                 "given for", " ", sQuote("split"))
        nkids <- length(breaks) + 1
    } else {
        nkids <- max(index, na.rm = TRUE)
    }
    prob <- rep(1, nkids) / nkids
    return(prob)
}

info_split <- function(split) {
    stopifnot(inherits(split, "partysplit"))
    split$info
}

kidids_split <- function(split, data, vmatch = 1:ncol(data), obs = NULL) {

    id <- varid_split(split)
    x <- data[[vmatch[id]]]
    if (!is.null(obs)) x <- x[obs]

    if (is.null(breaks_split(split))) {
        stopifnot(storage.mode(x) == "integer")
    } else {
        x <- as.integer(cut.default(as.numeric(x), 
                 breaks = c(-Inf, breaks_split(split), Inf), 
                 right = right_split(split)))
    }
    index <- index_split(split)
    ### empty factor levels correspond to NA and return NA here
    ### and thus the corresponding observations will be treated
    ### as missing values (surrogate or random splits):
    if (!is.null(index))
        x <- index[x]
    return(x)
}

character_split <- function(split, data = NULL, digits = getOption("digits") - 2) {

    varid <- varid_split(split)

    if (!is.null(data)) {
        ## names and labels
        lev <- lapply(data, levels)[[varid]]
        mlab <- names(data)[varid]

        ## determine split type
        type <- sapply(data, function(x) class(x)[1])[varid_split(split)]
        type[!(type %in% c("factor", "ordered"))] <- "numeric"
    } else {
        ## (bad) default names and labels
	lev <- NULL
	mlab <- paste("V", varid, sep = "")
	type <- "numeric"
    }

    ## process defaults for breaks and index
    breaks <- breaks_split(split)
    index <- index_split(split)
    right <- right_split(split)

    if (is.null(breaks)) breaks <- 1:(length(index) - 1)
    if (is.null(index)) index <- 1:(length(breaks) + 1)

    ## check whether ordered are really ordered
    if (type == "ordered") {
        if (length(breaks) > 1)
            type <- "factor"
    }
    ### <FIXME> format ordered multiway splits? </FIXME>

    switch(type, 
        "factor" = {
            nindex <- index[cut(seq_along(lev), c(-Inf, breaks, Inf), right = right)]
            dlab <- as.vector(tapply(lev, nindex, paste, collapse = ", "))
        },
        "ordered" = {
            if (length(breaks) == 1) {
                if (right)
                    dlab <- paste(c("<=", ">"), lev[breaks], sep = " ")
                else
                    dlab <- paste(c("<", ">="), lev[breaks], sep = " ")
            } else {
                stop("") ### see above
            }
            dlab <- dlab[index]
        },
        "numeric" = {
            breaks <- round(breaks, digits)
            if (length(breaks) == 1) {
                if (right)
                    dlab <- paste(c("<=", ">"), breaks, sep = " ")
                else
                    dlab <- paste(c("<", ">="), breaks, sep = " ")
            } else {
                dlab <- levels(cut(0, breaks = c(-Inf, breaks, Inf), 
                                   right = right))
            }
            dlab <- as.vector(tapply(dlab, index, paste, collapse = " | "))
        }
    )

    return(list(name = mlab, levels = dlab))
}
