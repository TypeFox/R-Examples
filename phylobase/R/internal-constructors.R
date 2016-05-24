
#####################
## Labels constructor
#####################

## (formerly) recursive function to have labels of constant length
## base = a character string
## n = number of labels
.genlab <- function(base, n) {
    if(n <= 0) return("")
    s <- seq(length.out=n)
    fw <- max(nchar(as.character(s)))
    numstr <- formatC(s, flag="0", width=fw)
    paste(base, numstr, sep="")
}

.createLabels <- function(value, ntips, nnodes, use.names = TRUE,
                          type = c("all", "tip", "internal")) {

    type <- match.arg(type)

    ## set up final length of object to return
    lgthRes <- switch(type, tip=ntips, internal=nnodes, all=ntips+nnodes)

    ## create NA character vector of node labels
    res <- character(lgthRes)
    is.na(res) <- TRUE

    ## create internal names
    names(res) <- switch(type,
                         tip = 1:ntips,
                         internal = seq(from=ntips+1, length=lgthRes),
                         all = 1:(ntips+nnodes))

    ## Convert empty labels to NA
    value[!nzchar(value)] <- NA

    ## if no values are provided
    if(missing(value) || is.null(value) || all(is.na(value))) {
        ## tip labels can't be NULL
        if(!identical(type, "internal")) {
            tipLbl <- .genlab("T", ntips)
            res[1:ntips] <- tipLbl
        }
    }

    ## if labels are provided
    else {
        ## check that lengths match
        if(length(value) != lgthRes)
            stop("Number of labels does not match number of nodes.")

        ## check if vector 'value' has name, and if so match with node.label names
        if(use.names && !is.null(names(value))) {
            if(!all(names(value) %in% names(res)))
                stop("Names provided don't match internal labels names.")
            res[match(names(value), names(res))] <- value
        }
        else
            res[1:lgthRes] <- value
    }

    res
}


.createEdge <- function(value, edgeMat, type=c("lengths", "labels"),
                        use.names=TRUE) {
    type <- match.arg(type)

    lgthRes <- nrow(edgeMat)
    res <- switch(type, lengths=numeric(lgthRes), labels=character(lgthRes))
    is.na(res) <- TRUE
    names(res) <- paste(edgeMat[,1], edgeMat[,2], sep="-")

    if(!(missing(value) || is.null(value) || all(is.na(value)))) {
        if(use.names && !is.null(names(value))) {
            if(!all(names(value) %in% names(res)))
                stop("Names provided don't match internal edge labels names.")
            res[match(names(value), names(res))] <- value
        }
        else
            res[1:lgthRes] <- value
    }

    res
}
