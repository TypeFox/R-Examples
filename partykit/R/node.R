
partynode <- function(id, split = NULL, kids = NULL, surrogates = NULL, info = NULL) {

    if (!is.integer(id) || length(id) != 1) {
        id <- as.integer(id0 <- id)
        if (any(is.na(id)) || !isTRUE(all.equal(id0, id)) || length(id) != 1)
            stop(sQuote("id"), " ", "must be a single integer")
    }

    if (is.null(split) != is.null(kids)) {
        stop(sQuote("split"), " ", "and", " ", sQuote("kids"), " ", 
             "must either both be specified or unspecified")
    }

    if (!is.null(kids)) {
        if (!(is.list(kids) && all(sapply(kids, inherits, "partynode"))) 
            || length(kids) < 2)
            stop(sQuote("kids"), " ", "must be an integer vector or a list of", 
                 " ", sQuote("partynode"), " ", "objects")
    }

    if (!is.null(surrogates)) {
        if (!is.list(surrogates) || any(!sapply(surrogates, inherits, "partysplit")))
            stop(sQuote("split"), " ", "is not a list of", " ", sQuote("partysplit"), 
                 " ", "objects")
    }

    node <- list(id = id, split = split, kids = kids, surrogates = surrogates, info = info)
    class(node) <- "partynode"
    return(node)
}

is.partynode <- function(x) {
    if (!inherits(x, "partynode")) return(FALSE)
    rval <- diff(nodeids(x, terminal = FALSE))
    isTRUE(all.equal(unique(rval), 1))
}

as.partynode <- function(x, ...)
    UseMethod("as.partynode")

as.partynode.partynode <- function(x, from = NULL, recursive = TRUE, ...) {
    if(is.null(from)) from <- id_node(x)
    from <- as.integer(from)
    if (!recursive) {
        if(is.partynode(x) & 
           id_node(x) == from) return(x)
    }
    id <- from - 1L
     
    new_node <- function(x) {
        id <<- id + 1L
        if(is.terminal(x)) return(partynode(id, info = info_node(x)))
        partynode(id,
             split = split_node(x),
             kids = lapply(kids_node(x), new_node),
             surrogates = surrogates_node(x),
             info = info_node(x))
    }
    
    return(new_node(x))    
}

as.partynode.list <- function(x, ...) {

    if (!all(sapply(x, inherits, what = "list")))
        stop("'x' has to be a list of lists")

    if (!all(sapply(x, function(x) "id" %in% names(x))))
        stop("each list in 'x' has to define a node 'id'")

    ok <- sapply(x, function(x) 
              all(names(x) %in% c("id", "split", "kids", "surrogates", "info")))
    if (any(!ok))
        sapply(which(!ok), function(i) 
            warning(paste("list element", i, "defines additional elements:", 
                          paste(names(x[[i]])[!(names(x[[i]]) %in% 
                                c("id", "split", "kids", "surrogates", "info"))], 
                                collapse = ", "))))
    
    ids <- as.integer(sapply(x, function(node) node$id))
    if(any(duplicated(ids))) stop("nodeids must be unique integers")
    x <- x[order(ids)]
    ids <- ids[order(ids)]

    new_recnode <- function(i) {
        x_i <- x[[which(ids == i)]]
        if (is.null(x_i$kids))
            partynode(id = x_i$id, info = x_i$info)
        else
            partynode(id = x_i$id, split = x_i$split,
                 kids = lapply(x_i$kids, new_recnode),
		 surrogates = x_i$surrogates,
                 info = x_i$info)
    }
    
    ret <- new_recnode(ids[1L])
    ### <FIXME> duplicates recursion but makes sure
    ###    that the ids are in pre-order notation with
    ###    from defined in as.partynode.partynode 
    ### </FIXME>
    as.partynode(ret, ...)
}

as.list.partynode <- function(x, ...)
{
    ids <- nodeids(x)
    obj <- vector(mode = "list", length = length(ids))
    thisnode <- NULL
    
    nodelist <- function(node) {
        if (is.terminal(node))
            obj[[which(ids == id_node(node))]] <<- list(id = id_node(node), info = info_node(node))
        else {
            thisnode <<- list(id = id_node(node), split = split_node(node),
                 kids = sapply(kids_node(node), function(k) id_node(k)))
             if (!is.null(surrogates_node(node)))
		 thisnode$surrogates <- surrogates_node(node)
             if (!is.null(info_node(node)))
		 thisnode$info <- info_node(node)
            obj[[which(ids == id_node(node))]] <<- thisnode
            lapply(kids_node(node), nodelist)
        }
    }
    nodelist(x)
    return(obj)
}


id_node <- function(node) {
    stopifnot(inherits(node, "partynode"))
    node$id
}

kids_node <- function(node) {
    stopifnot(inherits(node, "partynode"))
    node$kids
}

info_node <- function(node) {
    stopifnot(inherits(node, "partynode"))
    node$info
}

formatinfo_node <- function(node, FUN = NULL, default = "", prefix = NULL, ...) {
    info <- info_node(node)
    
    ## FIXME: better dispatch to workhorse FUN probably needed in the future, e.g.:
    ## (1) formatinfo() generic with formatinfo.default() as below,
    ## (2) supply default FUN from party$info$formatinfo() or similar.
    if(is.null(FUN)) FUN <- function(x, ...) {
      if(is.null(x)) x <- ""
      if(!is.object(x) & is.atomic(x)) x <- as.character(x)
      if(!is.character(x)) x <- capture.output(print(x), ...)
      x
    }
    
    info <- if(is.null(info)) default else FUN(info, ...)
    if(!is.null(prefix)) {
      info <- if(length(info) > 1L) c(prefix, info) else paste(prefix, info, sep = "")
    }
    info
}

### FIXME: permutation and surrogate splits: is only the primary
### variable permuted?
kidids_node <- function(node, data, vmatch = 1:ncol(data), obs = NULL, 
                        perm = NULL) {

    primary <- split_node(node)
    surrogates <- surrogates_node(node)

    ### perform primary split
    x <- kidids_split(primary, data, vmatch, obs)

    ### surrogate / random splits if needed
    if (any(is.na(x))) {
        ### surrogate splits
        if (length(surrogates) >= 1) {
            for (surr in surrogates) {
                nax <- is.na(x)
                if (!any(nax)) break;
                x[nax] <- kidids_split(surr, data, vmatch, obs = obs)[nax]
            }
        }
        nax <- is.na(x)
        ### random splits
        if (any(nax)) {
            prob <- prob_split(primary)
            x[nax] <- sample(1:length(prob), sum(nax), prob = prob, 
                             replace = TRUE)
        }
    }

    ### permute variable `perm' _after_ dealing with surrogates etc.
    if (!is.null(perm)) {
        if (varid_split(primary) %in% perm)
            return(sample(x))
    }
    return(x)
}

fitted_node <- function(node, data, vmatch = 1:ncol(data), 
                        obs = 1:nrow(data), perm = NULL) {

    ### should be equivalent to:
    # return(.Call("R_fitted_node", node, data, vmatch, as.integer(obs), 
    #             as.integer(perm)))

    if (is.logical(obs)) obs <- which(obs)
    if (is.terminal(node))
        return(rep(id_node(node), length(obs)))
    retid <- nextid <- kidids_node(node, data, vmatch, obs, perm)
    for (i in unique(nextid)) {
        indx <- nextid == i
        retid[indx] <- fitted_node(kids_node(node)[[i]], data,
                                   vmatch, obs[indx], perm)
    }
    return(retid)
}


length.partynode <- function(x)
    length(kids_node(x))

"[.partynode" <- "[[.partynode" <- function(x, i, ...) {
    stopifnot(length(i) == 1 & is.numeric(i))
    kids_node(x)[[i]]
}

split_node <- function(node) {
    stopifnot(inherits(node, "partynode"))
    node$split
}

surrogates_node <- function(node) {
    stopifnot(inherits(node, "partynode"))
    node$surrogates
}

is.terminal <- function(x, ...)
    UseMethod("is.terminal")

is.terminal.partynode <- function(x, ...) {
    kids <- is.null(kids_node(x))
    split <- is.null(split_node(x))
    stopifnot(kids == split)
    kids
}

## ## depth generic now taken from package 'grid'
## depth <- function(x, ...)
##     UseMethod("depth")

depth.partynode <- function(x, root = FALSE, ...) {
    if (is.terminal(x)) return(as.integer(root))
    max(sapply(kids_node(x), depth, root = root)) + 1L
}

width <- function(x, ...)
    UseMethod("width")

width.partynode <- function(x, ...) {
    if (is.terminal(x)) return(1)
    sum(sapply(kids_node(x), width.partynode))
}
