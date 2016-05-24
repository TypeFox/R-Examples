## FIXME: data in party
##   - currently assumed to be a data.frame
##   - potentially empty
##   - the following are all assumed to work:
##     dim(data), names(data)
##     sapply(data, class), lapply(data, levels)
##   - potentially these need to be modified if data/terms
##     should be able to deal with data bases

party <- function(node, data, fitted = NULL, terms = NULL, names = NULL, info = NULL) {

    stopifnot(inherits(node, "partynode"))
    stopifnot(inherits(data, "data.frame"))
    ### make sure all split variables are there 
    ids <- nodeids(node)[!nodeids(node) %in% nodeids(node, terminal = TRUE)]
    varids <- unique(unlist(nodeapply(node, ids = ids, FUN = function(x) 
        varid_split(split_node(x)))))
    stopifnot(varids %in% 1:ncol(data))
    
    if(!is.null(fitted)) {
        stopifnot(inherits(fitted, "data.frame"))
        stopifnot(nrow(data) == 0L | nrow(data) == nrow(fitted))
        
	# try to provide default variable "(fitted)"
	if(nrow(data) > 0L) {
          if(!("(fitted)" %in% names(fitted))) 
            fitted[["(fitted)"]] <- fitted_node(node, data = data)
	} else {
	  stopifnot("(fitted)" == names(fitted)[1L])
	}

        nt <- nodeids(node, terminal = TRUE)
        stopifnot(all(fitted[["(fitted)"]] %in% nt))

        node <- as.partynode(node, from = 1L)
        nt2 <- nodeids(node, terminal = TRUE)
        fitted[["(fitted)"]] <- nt2[match(fitted[["(fitted)"]], nt)]
    } else {
        node <- as.partynode(node, from = 1L)
	# default "(fitted)"
	if(nrow(data) > 0L & missing(fitted)) 
          fitted <- data.frame("(fitted)" = fitted_node(node, 
            data = data), check.names = FALSE)
    }
    
    party <- list(node = node, data = data, fitted = fitted, 
                  terms = NULL, names = NULL, info = info)
    class(party) <- "party"

    if(!is.null(terms)) {
        stopifnot(inherits(terms, "terms"))
	party$terms <- terms
    }

    if (!is.null(names)) {
        n <- length(nodeids(party, terminal = FALSE))
        if (length(names) != n)
            stop("invalid", " ", sQuote("names"), " ", "argument")
        party$names <- names
    }

    party
}

length.party <- function(x)
    length(nodeids(x))

names.party <- function(x)
    .names_party(x)

"names<-.party" <- function(x, value) {
     n <- length(nodeids(x, terminal = FALSE))
     if (!is.null(value) && length(value) != n)
         stop("invalid", " ", sQuote("names"), " ", "argument")
     x$names <- value
     x
}

.names_party <- function(party) {
    names <- party$names
    if (is.null(names))
        names <- as.character(nodeids(party, terminal = FALSE))
    names
}

node_party <- function(party) {
    stopifnot(inherits(party, "party"))
    party$node
}

is.constparty <- function(party) {
    stopifnot(inherits(party, "party"))
    if (!is.null(party$fitted)) 
        return(all(c("(fitted)", "(response)") %in% names(party$fitted)))
    return(FALSE)
}

as.constparty <- function(obj, ...) {
    if(!inherits(obj, "party")) obj <- as.party(obj)
    if (!is.constparty(obj)) {
        if(is.null(obj$fitted))
	  obj$fitted <- data.frame("(fitted)" = predict(obj, type = "node"), check.names = FALSE)
	if(!("(fitted)" %in% names(obj$fitted)))
	  obj$fitted["(fitted)"] <- predict(obj, type = "node")
	if(!("(response)" %in% names(obj$fitted)))
	  obj$fitted["(response)"] <- model.response(model.frame(obj))
	if(!("(weights)" %in% names(obj$fitted))) {
	  w <- model.weights(model.frame(obj))
	  if(is.null(w) && any(w != 1L)) obj$fitted["(weights)"] <- w
	}
    }
    if (is.constparty(obj)) {
        ret <- obj
        class(ret) <- c("constparty", class(obj))
        return(ret)
    }
    stop("cannot coerce object of class", " ", sQuote(class(obj)), 
          " ", "to", " ", sQuote("constparty"))
}

"[.party" <- "[[.party" <- function(x, i, ...) {
    if (is.character(i) && !is.null(names(x)))
        i <- which(names(x) %in% i)
    stopifnot(length(i) == 1 & is.numeric(i))
    stopifnot(i <= length(x) & i >= 1)
    i <- as.integer(i)
    dat <- data_party(x, i)
    if (!is.null(x$fitted)) {
        findx <- which("(fitted)" == names(dat))[1]
        fit <- dat[,findx:ncol(dat), drop = FALSE]
        dat <- dat[,-(findx:ncol(dat)), drop = FALSE]
        if (ncol(dat) == 0)
            dat <- x$data
    } else {
        fit <- NULL
        dat <- x$data
    }
    nam <- names(x)[nodeids(x, from = i, terminal = FALSE)]

    recFun <- function(node) {
        if (id_node(node) == i) return(node)
        kid <- sapply(kids_node(node), id_node)
        return(recFun(node[[max(which(kid <= i))]]))
    }
    node <- recFun(node_party(x))

    ret <- party(node = node, data = dat, fitted = fit, 
                 terms = x$terms, names = nam, info = x$info)
    class(ret) <- class(x)
    ret
}

nodeids <- function(obj, ...)
    UseMethod("nodeids")

nodeids.partynode <- function(obj, from = NULL, terminal = FALSE, ...) {

    if(is.null(from)) from <- id_node(obj)

    id <- function(node, record = TRUE, terminal = FALSE) {
      if(!record) return(NULL)
      if(!terminal)
          return(id_node(node))
      else
          if(is.terminal(node)) return(id_node(node)) else return(NULL)
    }

    rid <- function(node, record = TRUE, terminal = FALSE) {  
        myid <- id(node, record = record, terminal = terminal)
        if(is.terminal(node)) return(myid)
        kids <- kids_node(node)    
        kids_record <- if(record)  
            rep(TRUE, length(kids))
        else
            sapply(kids, id_node) == from
        return(c(myid,
            unlist(lapply(1:length(kids), function(i)
                rid(kids[[i]], record = kids_record[i], terminal = terminal)))
        ))
    }

    return(rid(obj, from == id_node(obj), terminal))
}

nodeids.party <- function(obj, from = NULL, terminal = FALSE, ...)
    nodeids(node_party(obj), from = from, terminal = terminal, ...)

nodeapply <- function(obj, ids = 1, FUN = NULL, ...)
    UseMethod("nodeapply")

nodeapply.party <- function(obj, ids = 1, FUN = NULL, by_node = TRUE, ...) {

    stopifnot(isTRUE(all.equal(ids, round(ids))))
    ids <- as.integer(ids)

    if(is.null(FUN)) FUN <- function(x, ...) x

    if (length(ids) == 0)
        return(NULL)

    if (by_node) {
        rval <- nodeapply(node_party(obj), ids = ids, FUN = FUN, ...)
    } else {
        rval <- lapply(ids, function(i) FUN(obj[[i]], ...))
    }

    names(rval) <- names(obj)[ids]
    return(rval)
}

nodeapply.partynode <- function(obj, ids = 1, FUN = NULL, ...) {

    stopifnot(isTRUE(all.equal(ids, round(ids))))
    ids <- as.integer(ids)

    if(is.null(FUN)) FUN <- function(x, ...) x

    if (length(ids) == 0)
        return(NULL)

    rval <- vector(mode = "list", length = length(ids))
    rval_id <- rep(0, length(ids))
    i <- 1
	
    recFUN <- function(node, ...) {
        if(id_node(node) %in% ids) {
            rval_id[i] <<- id_node(node)
            rval[[i]] <<- FUN(node, ...)
            i <<- i + 1
        }
        kids <- kids_node(node)
        if(length(kids) > 0) {
            for(j in 1:length(kids)) recFUN(kids[[j]])
        }
        invisible(TRUE)
    }
    foo <- recFUN(obj)
    rval <- rval[match(rval_id, ids)]
    return(rval)
}

predict.party <- function(object, newdata = NULL, perm = NULL, ...)
{
    ### compute fitted node ids first
    fitted <- if(is.null(newdata)) {    
        object$fitted[["(fitted)"]]	
    } else {
      terminal <- nodeids(object, terminal = TRUE)
	
      if(max(terminal) == 1L) {
        rep.int(1L, NROW(newdata))
      } else {
	
        inner <- 1L:max(terminal)
        inner <- inner[-terminal]

        primary_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
            varid_split(split_node(node))
        })
        surrogate_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
            surr <- surrogates_node(node)
            if(is.null(surr)) return(NULL) else return(sapply(surr, varid_split))
        })
        vnames <- names(object$data)

        ### the splits of nodes with a primary split in perm
        ### will be permuted
        if (!is.null(perm)) {
            stopifnot(all(perm %in% vnames))
            perm <- match(perm, vnames)
        }

        ## ## FIXME: the is.na() call takes loooong on large data sets
        ## unames <- if(any(sapply(newdata, is.na))) 
        ##     vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
        ## else 
        ##     vnames[unique(unlist(primary_vars))]
	unames <- vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
	
        vclass <- structure(lapply(object$data, class), .Names = vnames)
        ndnames <- names(newdata)
        ndclass <- structure(lapply(newdata, class), .Names = ndnames)
        checkclass <- all(sapply(unames, function(x) 
          isTRUE(all.equal(vclass[[x]], ndclass[[x]]))))
        factors <- sapply(unames, function(x) inherits(object$data[[x]], "factor"))
        checkfactors <- all(sapply(unames[factors], function(x) 
          isTRUE(all.equal(levels(object$data[[x]]), levels(newdata[[x]])))))
        ## FIXME: inform about wrong classes / factor levels?
        if(all(unames %in% ndnames) && checkclass && checkfactors) {
            vmatch <- match(vnames, ndnames)
            fitted_node(node_party(object), data = newdata, 
                        vmatch = vmatch, perm = perm)
        } else {
            if (!is.null(object$terms)) {
                mf <- model.frame(delete.response(object$terms), newdata)
                fitted_node(node_party(object), data = mf, 
                            vmatch = match(vnames, names(mf)), perm = perm)
            } else
                stop("") ## FIXME: write error message
        }
      }
    }
    ### compute predictions
    predict_party(object, fitted, newdata, ...)
}

predict_party <- function(party, id, newdata = NULL, ...)
    UseMethod("predict_party")

### do nothing expect returning the fitted ids
predict_party.default <- function(party, id, newdata = NULL, FUN = NULL, ...) {

    if (length(list(...)) > 1) 
        warning("argument(s)", " ", sQuote(names(list(...))), " ", "have been ignored")

    ## get observation names: either node names or
    ## observation names from newdata
    nam <- if(is.null(newdata)) {
      if(is.null(rnam <- rownames(data_party(party)))) names(party)[id] else rnam      
    } else {
      rownames(newdata)
    }
    if(length(nam) != length(id)) nam <- NULL

    if (!is.null(FUN))
        return(.simplify_pred(nodeapply(party, 
            nodeids(party, terminal = TRUE), FUN, by_node = TRUE), id, nam))

    ## special case: fitted ids
    return(structure(id, .Names = nam))
}

predict_party.constparty <- function(party, id, newdata = NULL,
    type = c("response", "prob", "quantile", "density", "node"),
    at = if (type == "quantile") c(0.1, 0.5, 0.9),
    FUN = NULL, simplify = TRUE, ...)
{
    ## extract fitted information
    response <- party$fitted[["(response)"]]
    weights <- party$fitted[["(weights)"]]
    fitted <- party$fitted[["(fitted)"]]
    if (is.null(weights)) weights <- rep(1, NROW(response))

    ## get observation names: either node names or
    ## observation names from newdata
    nam <- if(is.null(newdata)) names(party)[id] else rownames(newdata)
    if(length(nam) != length(id)) nam <- NULL

    ## match type
    type <- match.arg(type)

    ## special case: fitted ids
    if(type == "node")
      return(structure(id, .Names = nam))

    ### multivariate response
    if (is.data.frame(response)) {
        ret <- lapply(response, function(r) {
            ret <- .predict_party_constparty(node_party(party), fitted = fitted, 
                response = r, weights, id = id, type = type, at = at, FUN = FUN, ...)
            if (simplify) .simplify_pred(ret, id, nam) else ret
        })
        if (all(sapply(ret, is.atomic)))
            ret <- as.data.frame(ret)
        names(ret) <- colnames(response)
        return(ret)
    }

    ### univariate response
    ret <- .predict_party_constparty(node_party(party), fitted = fitted, response = response, 
        weights = weights, id = id, type = type, at = at, FUN = FUN, ...)
    if (simplify) .simplify_pred(ret, id, nam) else ret[as.character(id)]
}

### functions for node prediction based on fitted / response
.pred_Surv <- function(y, w) {
    if (length(y) == 0) return(NA)
    survfit(y ~ 1, weights = w, subset = w > 0)
}

.pred_Surv_response <- function(y, w) {
    if (length(y) == 0) return(NA)
    .median_survival_time(.pred_Surv(y, w))
}
 
.pred_factor <- function(y, w) {
    lev <- levels(y)
    sumw <- tapply(w, y, sum)
    sumw[is.na(sumw)] <- 0
    prob <- sumw / sum(w)
    names(prob) <- lev
    return(prob)
}

.pred_factor_response <- function(y, w) {
    prob <- .pred_factor(y, w)
    return(factor(which.max(prob), levels = 1:nlevels(y),
                  labels = levels(y), 
                  ordered = is.ordered(y)))
    return(prob) 
}
                    
.pred_numeric_response <- function(y, w) 
    weighted.mean(y, w, na.rm = TRUE)

.pred_ecdf <- function(y, w) {
    if (length(y) == 0) return(NA)
    iw <- as.integer(round(w))
    if (max(abs(w - iw)) < sqrt(.Machine$double.eps)) {
        y <- rep(y, w)
        return(ecdf(y))
    } else {
        stop("cannot compute empirical distribution function with non-integer weights")
    }
}

.pred_quantile <- function(y, w) {
    y <- rep(y, w)
    function(p, ...) quantile(y, probs = p, ...)
}

.pred_density <- function(y, w) {
    d <- density(y, weights = w / sum(w))
    approxfun(d[1:2], rule = 2)
}

### workhorse: compute predictions based on fitted / response data
.predict_party_constparty <- function(node, fitted, response, weights,
    id = id, type = c("response", "prob", "quantile", "density"),
    at = if (type == "quantile") c(0.1, 0.5, 0.9), FUN = NULL, ...) {

    type <- match.arg(type)
    if (is.null(FUN)) {

        rtype <- class(response)[1]
        if (rtype == "ordered") rtype <- "factor"    
        if (rtype == "integer") rtype <- "numeric"

        if (type %in% c("quantile", "density") && rtype != "numeric")
            stop("quantile and density estimation currently only implemented for numeric responses")

        FUN <- switch(rtype,
            "Surv" = if (type == "response") .pred_Surv_response else .pred_Surv,
            "factor" = if (type == "response") .pred_factor_response else .pred_factor,
            "numeric" = switch(type,
                "response" = .pred_numeric_response,
                "prob" = .pred_ecdf,
                "quantile" = .pred_quantile, 
                "density" = .pred_density) 
       )
    }
      
    ## empirical distribution in each leaf
    if (all(id %in% fitted)) {
        tab <- tapply(1:NROW(response), fitted, 
                      function(i) FUN(response[i], weights[i]), simplify = FALSE)
    } else {
        ### id may also refer to inner nodes
        tab <- as.array(lapply(sort(unique(id)), function(i) {
            index <- fitted %in% nodeids(node, i, terminal = TRUE)
            ret <- FUN(response[index], weights[index])
            ### no information about i in fitted
            if (all(!index)) ret[1] <- NA
            return(ret)
        }))
        names(tab) <- as.character(sort(unique(id)))
    }
    if (inherits(tab[[1]], "function") && !is.null(at))
        tab <- lapply(tab, function(f) f(at))
    tn <- names(tab)
    dim(tab) <- NULL
    names(tab) <- tn

    tab
}


### simplify structure of predictions
.simplify_pred <- function(tab, id, nam) {

    if (all(sapply(tab, length) == 1) & all(sapply(tab, is.atomic))) {
        ret <- do.call("c", tab)
        names(ret) <- names(tab)
        ret <- if (is.factor(tab[[1]]))
            factor(ret[as.character(id)], levels = 1:length(levels(tab[[1]])),
		   labels = levels(tab[[1]]), ordered = is.ordered(tab[[1]]))
        else 
            ret[as.character(id)]
        names(ret) <- nam
    } else if (length(unique(sapply(tab, length))) == 1 & 
               all(sapply(tab, is.numeric))) {
        ret <- matrix(unlist(tab), nrow = length(tab), byrow = TRUE)
        colnames(ret) <- names(tab[[1]])
        rownames(ret) <- names(tab)
        ret <- ret[as.character(id),, drop = FALSE]
        rownames(ret) <- nam
    } else {
        ret <- tab[as.character(id)]
        names(ret) <- nam
    }
    ret
}

data_party <- function(party, id = 1L)
    UseMethod("data_party")

data_party.default <- function(party, id = 1L) {
    
    extract <- function(id) {
        if(is.null(party$fitted))
            if(nrow(party$data) == 0) return(NULL)
        else
            stop("cannot subset data without fitted ids")

        ### which terminal nodes follow node number id?
        nt <- nodeids(party, id, terminal = TRUE)
        wi <- party$fitted[["(fitted)"]] %in% nt

        ret <- if (nrow(party$data) == 0)
            subset(party$fitted, wi)
        else
            subset(cbind(party$data, party$fitted), wi)
        ret
    }
    if (length(id) > 1)
        return(lapply(id, extract))
    else 
        return(extract(id))
}

width.party <- function(x, ...) {
  width(node_party(x), ...)
}

depth.party <- function(x, root = FALSE, ...) {
  depth(node_party(x), root = root, ...)
}

getCall.party <- function(x, ...) {
  x$info$call
}

nodeprune <- function(x, ids, ...)
    UseMethod("nodeprune")

nodeprune.party <- function(x, ids, ...) {

    ### map names to nodeids
    if (!is.numeric(ids))
        ids <- match(ids, names(x))
    stopifnot(ids %in% nodeids(x))

    ### compute indices path to each node
    ### to be pruned off
    idxs <- lapply(ids, .get_path, obj = node_party(x))

    ### [[.party is NOT [[.list
    cls <- class(x)
    x <- unclass(x)
    ni <- which(names(x) == "node")

    for (i in 1:length(idxs)) {
    
        idx <- c(ni, idxs[[i]])
        ### check if we already pruned-off this node
        tmp <- try(x[[idx]], silent = TRUE)
        if (inherits(tmp, "try-error"))
            next()

        ### node ids of off-pruned daugther nodes
        idrm <- nodeids(x[[idx]])[-1]

        ### prune node by introducing a "new" terminal node
        x[[idx]] <- partynode(id = id_node(x[[idx]]),
                              info = info_node(x[[idx]]))

        ### constparty only: make sure the node ids in
        ### fitted are corrected
        if (length(idrm) > 0) {
             if(!is.null(x$fitted) && 
                 "(fitted)" %in% names(x$fitted)) {
                     j <- x$fitted[["(fitted)"]] %in% idrm
                     x$fitted[["(fitted)"]][j] <- ids[i]
             }
        }
    }

    ### reindex to 1:max(nodeid)
    class(x) <- cls
    nodeids(x) <- 1:length(nodeids(x))
    return(x)
}

nodeprune.partynode <- function(x, ids, ...) {

    stopifnot(ids %in% nodeids(x))

    ### compute indices path to each node
    ### to be pruned off
    idxs <- lapply(ids, .get_path, obj = x)

    ### [[.partynode is NOT [[.list
    cls <- class(x)
    x <- unclass(x)

    for (i in 1:length(idxs)) {
        ## path to be pruned
        idx <- idxs[[i]]
	if(!is.null(idx)) {
          ### check if we already pruned-off this node
          tmp <- try(x[[idx]], silent = TRUE)
          if(inherits(tmp, "try-error")) next()
          ### prune node by introducing a "new" terminal node
          x[[idx]] <- partynode(id = id_node(tmp), info = info_node(tmp))
	} else {
	  ## if idx path is NULL prune everything
	  x[2L:4L] <- NULL
	}
    }

    class(x) <- cls
    return(as.partynode(x, from = 1L))
}

nodeprune.default <- function(x, ids, ...)
    stop("No", sQuote("nodeprune"), "method for class", class(x), "implemented")

.list.rules.party <- function(x, i = NULL, ...) {
    if (is.null(i)) i <- nodeids(x, terminal = TRUE)
    if (length(i) > 1) {
        ret <- sapply(i, .list.rules.party, x = x)
        names(ret) <- if (is.character(i)) i else names(x)[i]
        return(ret)
    }
    if (is.character(i) && !is.null(names(x)))
        i <- which(names(x) %in% i)
    stopifnot(length(i) == 1 & is.numeric(i))
    stopifnot(i <= length(x) & i >= 1)
    i <- as.integer(i)
    dat <- data_party(x, i)  
    if (!is.null(x$fitted)) {
        findx <- which("(fitted)" == names(dat))[1]  
        fit <- dat[,findx:ncol(dat), drop = FALSE]   
        dat <- dat[,-(findx:ncol(dat)), drop = FALSE]
        if (ncol(dat) == 0)
            dat <- x$data
    } else {
        fit <- NULL  
        dat <- x$data
    }

    rule <- c()

    recFun <- function(node) {
        if (id_node(node) == i) return(NULL)   
        kid <- sapply(kids_node(node), id_node)
        whichkid <- max(which(kid <= i))
        split <- split_node(node)
        ivar <- varid_split(split)
        svar <- names(dat)[ivar]
        index <- index_split(split)
        if (is.factor(dat[, svar])) {
            slevels <- levels(dat[, svar])[index == whichkid]
            srule <- paste(svar, " %in% c(\"", 
                paste(slevels, collapse = "\", \"", sep = ""), "\")",
                sep = "")
        } else {
            if (is.null(index)) index <- 1:length(kid)
            breaks <- cbind(c(-Inf, breaks_split(split)), 
                            c(breaks_split(split), Inf))
            sbreak <- breaks[index == whichkid,]
            right <- right_split(split)
            srule <- c()
            if (is.finite(sbreak[1]))
                srule <- c(srule, 
                    paste(svar, ifelse(right, ">", ">="), sbreak[1]))
            if (is.finite(sbreak[2]))
                srule <- c(srule, 
                    paste(svar, ifelse(right, "<=", "<"), sbreak[2]))
            srule <- paste(srule, collapse = " & ")
        }
        rule <<- c(rule, srule)
        return(recFun(node[[whichkid]]))
    }
    node <- recFun(node_party(x))
    paste(rule, collapse = " & ")
}
