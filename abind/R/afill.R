"afill<-" <- function(x, ..., excess.ok=FALSE, local=TRUE, value) UseMethod("afill<-")

"afill<-.default" <- function(x, ..., excess.ok=FALSE, local=TRUE, value) {
    # The idea of afill<- is that some of the indices can be specified on the LHS
    # and the others are taken from the dimnames of the RHS, e.g., if length(dim(x))==4
    # and y is a matrix, then
    #   afill(x, , 2:3, , T) <- y
    # will use the indices list(rownames(y), 2:3, colnames(y), T) to assign
    # the values of y into x.
    x.dn <- if (length(dim(x))) dimnames(x) else list(names(x))
    x.d <- if (length(dim(x))) dim(x) else length(x)
    value.dn <- if (length(dim(value))) dimnames(value) else list(names(value))
    value.d <- if (length(dim(value))) dim(value) else length(value)
    # to find the empty anon args, must work with the unevaluated dot args
    dot.args.uneval <- match.call(expand.dots=FALSE)$...
    if (length(dot.args.uneval))
        missing.dot.args <- sapply(dot.args.uneval, function(arg) is.symbol(arg) && as.character(arg)=="")
    else
        missing.dot.args <- logical(0)
    if (length(value.d) < length(x.d)) {
        if (length(dot.args.uneval)==0) {
            stop("must supply anonymous args to show how value matches up with x")
        } else {
            if (sum(!missing.dot.args) + max(length(value.d),1) != max(length(x.d),1))
                stop("must have ", length(x.d)-max(length(value.d),1),
                     " non-missing anon args to assign a ", max(length(value.d), 1),
                     "-d value into a ", max(length(x.d), 1), "-d array")
        }
    } else if (length(value.d) == length(x.d)) {
        if (length(dot.args.uneval)==0) {
            dot.args.uneval <- vector("list", length(x.d))
            missing.dot.args <- rep(TRUE, length(x.d))
        } else if (length(dot.args.uneval) != length(x.d)) {
            stop("must have 0 or ", length(x.d), " empty arguments when 'x' and 'value' have same number of dims")
        }
    } else {
        stop("does not make sense to have more dims in value than x")
    }
    if (any(missing.dot.args) && (is.null(x.dn) || any(missing.dot.args & sapply(x.dn, length)==0 & x.d!=0)))
        stop("'x' must have names on dimensions corresponding to those in 'value'")
    if (any(missing.dot.args) && (is.null(value.dn) || any(sapply(value.dn, length)==0 & value.d!=0)))
        stop("'value' must have names on dimensions corresponding to empty args on the LHS")
    # Now we can work with evaluated dot args.
    # Can't do dot.args <- list(...) because that will
    # stop with an error for missing args.
    dot.args <- mapply(dot.args.uneval, missing.dot.args, FUN=function(arg, m) if (!m) eval(arg) else NULL)
    # construct the numeric indices
    idxs <- vector("list", length(x.d))
    strip.excess <- FALSE
    for (i in seq(len=length(x.d))) {
        if (missing.dot.args[i]) {
            j <- cumsum(missing.dot.args)[i] # dim-num in value
            idxs[[i]] <- match(value.dn[[j]], x.dn[[i]], nomatch=0)
            if (any(idxs[[i]]==0)) {
                if (!excess.ok)
                    stop("value has dimnames that are not in 'x' on dim[", i, "]: ",
                         paste("'", value.dn[[j]][which(idxs[[i]]==0)[min(3, sum(idxs[[i]]==0))]],
                               "'", sep="", collapse=", "), if (sum(idxs[[i]]==0)>3) " ...")
                strip.excess <- TRUE
            }
        } else {
            if (is.character(dot.args[[i]])) {
                if (length(x.dn[[i]]) != x.d[i])
                    stop("'x' doesn't have dimnames on dim ", i)
                idxs[[i]] <- match(dot.args[[i]], x.dn[[i]], nomatch=NA)
                if (any(is.na(idxs[[i]])))
                    stop("LHS character indicies at on dim ", i, " not matched: ", paste("'", dot.args[[i]][which(is.na(idxs[[i]]))[seq(len=min(3, sum(is.na(idxs[[i]]))))]], "'", collapse=", "))
            } else if (is.logical(dot.args[[i]])) {
                idxs[[i]] <- seq(len=x.d[i])[dot.args[[i]]]
                if (any(is.na(idxs[[i]])))
                    stop("LHS logical indicies at on dim ", i, " have NA value")
            } else if (is.numeric(dot.args[[i]]) & all(dot.args[[i]] >= 0)) {
                if (any(ii <- dot.args[[i]] == 0))
                    idxs[[i]] <- dot.args[[i]][!ii]
                else
                    idxs[[i]] <- dot.args[[i]]
                if (any(idxs[[i]] > x.d[i]))
                    stop("LHS numeric indicies at on dim ", i, " values too large")
            } else if (is.numeric(dot.args[[i]]) & all(dot.args[[i]] <= 0)) {
                idxs[[i]] <- seq(len=x.d[i])[dot.args[[i]]]
            } else {
                stop("LHS args for indices at dim ", i, " must be character, logical, numeric>0 or numeric<=0")
            }
        }
    }
    if (strip.excess) {
        value <- eval(as.call(c(list(as.name("["), as.name("value")), lapply(idxs[missing.dot.args], function(i) which(i!=0)))))
    }
    # replicate value appropriately if needed
    # look in the examples/tests for afill for the 4-d case for an example
    # that explains the logical here.
    if (prod(value.d)>1 && length(unique(value))>1) {
        j <- 0
        need.rep <- 1
        for (i in seq(along=missing.dot.args)) {
            if (missing.dot.args[i]) {
                j <- j+1
                if (need.rep > 1)
                    value <- asub(value, rep(seq(len=value.d[j]), each=need.rep), dims=j)
                need.rep <- 1
            } else {
                need.rep <- need.rep * length(idxs[[i]])
            }
        }
    }
    if (length(value)) {
        # Construct a skeleton call that we can pull an empty arg out of (xic[[3]])
        xic <- Quote(x[,drop=drop])
        # Find the name of x in the caller's frame
        x.caller <- substitute(x)
        if (local || !is.name(x.caller)) {
            # Evaluate the assignment in the frame of the function.  This
            # will create a duplicate of 'x', but trying to evaluate in
            # the frame of the caller is tricky...
            subcall <- call("<-", as.call(c(list(as.name("["), as.name("x")), idxs)), as.name("value"))
            if (length(i <- which(sapply(idxs, is.null))+2))
                subcall[[2]][i] <- xic[[3]]
            eval(subcall)
            return(x)
        } else {
            # Attempt to evaluate in the frame of the caller
            subcall <- call("<-", as.call(c(list(as.name("["), x.caller), idxs)), value)
            if (length(i <- which(sapply(idxs, is.null))+2))
                subcall[[2]][i] <- xic[[3]]
            eval(subcall, sys.parent(1))
            return(eval(x.caller, sys.parent(1)))
        }
    }
    x
}
