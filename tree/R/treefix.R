#
# file tree/R/treefix.R copyright (C) 1994-2006 B. D. Ripley
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
prune.tree <-
    function(tree, k = NULL, best = NULL, newdata, nwts,
             method = c("deviance", "misclass"),
             loss = 1 - diag(nc), eps = 1e-3)
{
    if(inherits(tree, "singlenode")) stop("can not prune singlenode tree")
    if(!inherits(tree, "tree")) stop("not legitimate tree")
    method <- match.arg(method)
    nc <- length(attr(tree, "ylevels"))
    if(method == "misclass" && !nc)
        stop("misclass only for classification trees")
    frame <- tree$frame
    node <- row.names(frame)
    nodes <- as.numeric(node)
    nnode <- length(node)
    ndim <- ceiling(nnode/2)

    if(is.null(y <- tree$y))
        y <- model.extract(model.frame(tree), "response")
    if(is.null(w <- tree$weights))
        w <- model.extract(model.frame(tree), "weights")
    if(is.null(w)) w <- rep(1, length(y))
    if(method == "misclass") {
        Z <- .C(VR_dev1,
                as.integer(nnode),
                as.integer(nodes),
                integer(nnode),
                dev = double(nnode),
                sdev = double(nnode),
                as.integer(y),
                as.integer(length(y)),
                as.integer(frame$yval),
                as.integer(tree$where),
                as.double(w),
                as.integer(nc), as.double(loss))
        dev <- Z$dev; sdev <- Z$sdev
    } else {
        dev <- tree$frame$dev
        sdev <- if(!nc) {
            .C(VR_dev3,
               as.integer(nnode),
               as.integer(nodes),
               integer(nnode),
               dev = double(nnode),
               sdev = double(nnode),
               as.double(y),
               as.integer(length(y)),
               as.double(frame$yval),
               as.integer(tree$where),
               as.double(w))$sdev
        } else  {
            -2 * .C(VR_dev2,
                    as.integer(nnode),
                    as.integer(nodes),
                    integer(nnode),
                    dev = double(nnode),
                    sdev = double(nnode),
                    as.integer(y),
                    as.integer(length(y)),
                    as.double(frame$yprob),
                    as.integer(tree$where),
                    as.double(w))$sdev
        }
    }
    if(missing(newdata) || is.null(newdata)) {
        ndev <- dev
        nsdev <- sdev
    } else {
        if(is.null(attr(newdata, "terms")))
            nd <- model.frame(tree$terms, newdata, na.action=na.pass,
                              xlev = tree$xlevels)
        else nd <- newdata
        if (!is.null(cl <- attr(tree$terms, "dataClasses")))
            .checkMFClasses(cl, nd)
        y <- model.extract(nd, "response")
        if(missing(nwts)) nwts <- rep(1, length(y))
        where <- pred1.tree(tree, tree.matrix(nd))
        if(method == "misclass") {
            Z <- .C(VR_dev1,
                    as.integer(nnode),
                    as.integer(nodes),
                    integer(nnode),
                    dev = double(nnode),
                    sdev = double(nnode),
                    as.integer(y),
                    as.integer(length(y)),
                    as.integer(frame$yval),
                    as.integer(where),
                    as.double(nwts),
                    as.integer(nc), as.double(loss))
            ndev <- Z$dev; nsdev <- Z$sdev
        } else {
            if(!nc) {
                Z <- .C(VR_dev3,
                        as.integer(nnode),
                        as.integer(nodes),
                        integer(nnode),
                        dev = double(nnode),
                        sdev = double(nnode),
                        as.double(y),
                        as.integer(length(y)),
                        as.double(frame$yval),
                        as.integer(where),
                        as.double(nwts))
                ndev <- Z$dev; nsdev <- Z$sdev
            } else {
                yp <- frame$yprob
                yp[yp==0] <- max(0,eps)
                Z <- .C(VR_dev2,
                        as.integer(nnode),
                        as.integer(nodes),
                        integer(nnode),
                        dev = double(nnode),
                        sdev = double(nnode),
                        as.integer(y),
                        as.integer(length(y)),
                        as.double(yp),
                        as.integer(where),
                        as.double(nwts))
                ndev <- -2 * Z$dev; nsdev <- -2 *Z$sdev
            }
        }
    }
    zp <- .C(VR_prune2,
             n=as.integer(nnode),
             as.integer(nodes),
             as.integer(frame$var == "<leaf>"),
             as.double(dev), as.double(sdev),
             as.double(ndev), as.double(nsdev),
             keep=integer(nnode),
             as.integer(order(nodes)),
             double(nnode),
             integer(nnode),
             double(nnode),
             alpha=double(ndim),
             inode=integer(ndim),
             size=integer(ndim),
             deviance=double(ndim),
             newdev=double(ndim))
    n <- zp$n
    alpha <- zp$alpha[1L:n]
    size <- zp$size[1L:n]
    index <- 0L
    if(missing(k) || is.null(k)) {
        ind <- drop(outer(unique(alpha), alpha, ">=") %*% rep(1L, length(alpha)))
        k <- alpha[ind]
        k[1L] <- -Inf
        deviance <- zp$newdev[ind]
        size <- size[ind]
        if(!missing(best) && !is.null(best)) {
            index <- ind[sum(best <= size)]
            if(length(index) == 0L) {
                warning("best is bigger than tree size")
                index <- 1L
            }
        }
    } else {
        if(length(k) == 1) index <- sum(k >= alpha)
        else {
            k <- pmax(k, -1e+100)
            ind <- drop(outer(k, alpha, ">=") %*% rep(1, length(alpha)))
            deviance <- zp$newdev[ind]
            size <- size[ind]
        }
    }
    if(index == 1L) return(tree)
    if(index > 1L) {
        pnodes <- zp$inode[-1L]
        tree <- snip.tree(tree, pnodes[seq(index-1)])
        tree$call$tree <- match.call()$tree
        return(tree)
    }
    obj <- list(size = size, dev = deviance, k = k, method = method)
    class(obj) <- c("prune", "tree.sequence")
    obj
}

predict.tree <-
    function(object, newdata = list(),
             type = c("vector", "tree", "class", "where"),
             split = FALSE, nwts, eps = 1e-3, ...)
{
    which.is.max <- function(x)
    {
        y <- seq(length(x))[x == max(x)]
        if(length(y) > 1L) sample(y, 1)
        else y
    }

    pred2.tree  <- function(tree, x)
    {
        frame <- tree$frame
        if(!length(frame$yprob)) stop("only for classification trees")
        dimx <- dim(x)
        ypred <- .C(VR_pred2,
                    as.double(x),
                    as.integer(unclass(frame$var) - 1),#0 denotes leaf node
                    as.character(frame$splits[, "cutleft"]),
                    as.character(frame$splits[, "cutright"]),
                    as.integer(sapply(attr(tree, "xlevels"), length)),
                    as.integer(row.names(frame)),
                    as.integer(frame$n),
                    as.integer(nf <- dim(frame)[1L]),
                    as.integer(dimx[1L]),
                    where = double(nf*dimx[1L]),
                    NAOK = TRUE)
        ypred <- matrix(ypred$where, nf)
        dimnames(ypred) <- list(row.names(frame),dimnames(x)[[1L]])
        ypred
    }

    if(!inherits(object, "tree") && !inherits(object, "singlenode"))
        stop("not legitimate tree")
    type <- match.arg(type)
    if(type == "class" && is.null(attr(object, "ylevels")))
        stop("type \"class\" only for classification trees")
    if((missing(newdata) || is.null(newdata)) && type == "tree")
        return(object)                  #idiot proofing
    if(missing(newdata) || is.null(newdata)) {
        where <- object$where
        newdata <- model.frame(object)
        if(!is.null(object$call$weights))
            nwts <- model.extract(model.frame(object), "weights")
    } else {
        if(is.null(attr(newdata, "terms"))) {
            # newdata is not a model frame.
            Terms <- object$terms
            if(type == "tree") {
                # test if response can be extracted from newdata
                response.vars <- all.vars(formula(Terms)[[2L]])
                response.exists <-
                    sapply(response.vars, function(nm, newdata)
                           eval(substitute(exists(nm), list(nm=nm)),
                                envir = newdata),
                           newdata)
                if(!all(response.exists)) Terms <- delete.response(Terms)
            } else Terms <- delete.response(Terms)
            newdata <- model.frame(Terms, newdata, na.action = na.pass,
                                   xlev = object$xlevels)
            if (!is.null(cl <- attr(Terms, "dataClasses")))
                .checkMFClasses(cl, newdata)
        }
        where <- pred1.tree(object, tree.matrix(newdata))
    }
    if(type == "where") return(where)
    frame <- object$frame
    node <- row.names(frame)
    nodes <- as.numeric(node)
    nnode <- length(node)
    if(type != "tree")
        if(is.null(lev <- attr(object, "ylevels"))) {
            if(!split) {
                frame <- frame$yval[where]
                names(frame) <- names(where)
                return(frame)
            } else {
                where <- pred2.tree(object, tree.matrix(newdata))
                leaf <- frame$var=="<leaf>"
                frame <- t(where[leaf, , drop = FALSE]) %*% frame$y[leaf]
                names(frame) <- names(where)
                return(frame)
            }
        } else {
            if(!split) {
                pr <- frame$yprob[where,  , drop = FALSE]
                dimnames(pr)[[1L]] <- names(where)
            } else {
                where <- pred2.tree(object, tree.matrix(newdata))
                leaf <- frame$var=="<leaf>"
                pr <- t(where[leaf,,drop = FALSE]) %*% frame$yprob[leaf,,drop=FALSE]
                dimnames(pr) <- list(names(where), lev)
            }
            if(type=="class") {
                cl <- apply(pr, 1L, which.is.max)
                return(factor(lev[cl], levels=lev))
            } else return(pr)
        }
    # now must be type = "tree"
    which <- descendants(as.numeric(row.names(frame)))[, where, drop = FALSE]
    if(!all(response.exists)) dev <- rep(NA, nrow(frame))
    else {
        y <- model.extract(newdata, "response")
        if(missing(nwts)) nwts <- rep(1, length(y))
        if(!length(attr(object, "ylevels"))) {
#
#  handle NAs in y separately.
#
            drp <- is.na(y); nwts[drp] <- 0; y[drp] <- 0
            dev <- .C(VR_dev3,
                      as.integer(nnode),
                      as.integer(nodes),
                      integer(nnode),
                      dev = double(nnode),
                      sdev = double(nnode),
                      as.double(y),
                      as.integer(length(y)),
                      as.double(frame$yval),
                      as.integer(where),
                      as.double(nwts))$dev
            dev[which %*% drp > 0] <- NA
        } else {
            yp <- frame$yprob
            yp[yp==0] <- max(0,eps)
            drp <- is.na(y); nwts[drp] <- 0; y[drp] <- levels(y)[1L]
            dev <- -2 * .C(VR_dev2,
                           as.integer(nnode),
                           as.integer(nodes),
                           integer(nnode),
                           dev = double(nnode),
                           sdev = double(nnode),
                           as.integer(y),
                           as.integer(length(y)),
                           as.double(yp),
                           as.integer(where),
                           as.double(nwts))$dev
            dev[which %*% drp > 0] <- NA
        }
    }
    object$frame$dev <- as.vector(dev)
    object$frame$n <- as.vector(which %*% rep(1, length(where)))
    object$where <- where
    object$call <- match.call()
    object$y <- object$x <- NULL
    object
}

pred1.tree <- function(tree, x)
{
    frame <- tree$frame
    dimx <- dim(x)
    ypred <- .C(VR_pred1,
                as.double(x),
                as.integer(unclass(frame$var) - 1),#0 denotes leaf node
                as.character(frame$splits[, "cutleft"]),
                as.character(frame$splits[, "cutright"]),
                as.integer(sapply(attr(tree, "xlevels"), length)),
                as.integer(row.names(frame)),
                as.integer(frame$n),
                as.integer(dim(frame)[1L]),
                as.integer(dimx[1L]),
                as.integer(dimx[2L]),
                where = integer(dimx[1L]),
                NAOK = TRUE)
    ypred <- ypred$where
    names(ypred) <- dimnames(x)[[1L]]
    ypred
}


na.tree.replace <- function(frame)
{
    if(!is.null(j <- attr(attr(frame, "terms"), "response"))) {
        x <- frame[[j]]
        pos <- is.na(x)
        if(any(pos)) {
            frame <- frame[!pos,  , drop = FALSE]
            warning(sum(pos),
                    " observations omitted due to missing values in the response")
        }
    }
    if(!is.na(j <- match("(weights)", names(frame)))) {
        x <- frame[[j]]
        pos <- is.na(x)
        if(any(pos)) {
            frame <- frame[!pos,  , drop = FALSE]
            warning(sum(pos),
                    " observations omitted due to missing values in the supplied weights")
        }
    }
    vars <- names(frame)
    names(vars) <- vars
    for(j in names(frame)) {
        x <- frame[[j]]
        pos <- is.na(x)
        if(any(pos)) {
            if(!length(levels(x)))
                stop("continuous variable ", j, " contained NAs")
            else {
                cl <- class(x)
                class(x) <- NULL
                lev <- c(attr(x, "levels"), "NA")
                x[pos] <- length(lev)
                levels(x) <- lev
                class(x) <- cl
            }
            frame[[j]] <- x
        }
    }
    frame
}

prune.misclass <- function(tree, k = NULL, best = NULL, newdata,
                           nwts, loss, eps = 1e-3)
{
    oc <- match.call()
    oc$method <- "misclass"
    oc[[1L]] <- as.name("prune.tree")
    eval.parent(oc)
}
