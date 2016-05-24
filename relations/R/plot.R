### * plot.relation

plot.relation <-
function(x,
         attrs = list(graph = list(rankdir = "BT"),
                      edge = list(arrowsize = NULL),
                      node = list(shape = "rectangle", fixedsize = FALSE)),
         limit = 6L, labels = NULL, main = NULL,
         type = c("simplified", "raw"),
         ...)
{
    type <- match.arg(type)
    plot(as.relation_ensemble(x),
         attrs = list(attrs), type = type, limit = limit, labels = list(labels), ..., main = main)
}


.make_unique_labels <-
function(x)
{
    I <- which(duplicated(x))
    x[I] <- paste(x[I], ".", seq_along(I), sep="")
    x
}


### * plot.relation_ensemble

plot.relation_ensemble <-
function(x, attrs = list(list(graph = list(rankdir = "BT"),
                              edge = list(arrowsize = NULL),
                              node = list(shape = "rectangle",
                                          fixedsize = FALSE))),
         type = "simplified",
         limit = 6L, labels = NULL,
         ..., layout = NULL, main = NULL)
{
    if(!is.relation_ensemble(x))
        stop("Wrong class.")
    if(!all(sapply(x,
                   function(e)
                   (isTRUE(relation_is_crisp(e)) &&
                    relation_is_endorelation(e)))))
        stop("Plotting only available for ensembles of crisp endorelations without missings.")

    ## Make things a bit more efficient.
    x <- unclass(x)
    if(system.file(package = "Rgraphviz") == "")
        stop("Plotting requires package 'Rgraphviz'.")

    ## Number of elements.
    n <- length(x)

    ## expand and check type argument.
    type <- sapply(rep(type, length.out = n),
                   match.arg, c("simplified", "raw"))

    ## expand limit
    limit <- rep(limit, length.out = n)

    ## expand main
    if(is.null(main))
        main <- rep.int(NA_character_, n)

    ## Layout.
    if (n > 1L) {
        byrow <- TRUE
        if(is.null(layout)) {
            nc <- ceiling(sqrt(n))
            nr <- ceiling(n / nc)
        }
        else {
            layout <- c(as.list(layout), byrow)[seq_len(3)]
            if(is.null(names(layout)))
                names(layout) <- c("nr", "nc", "byrow")
            nr <- layout[["nr"]]
            nc <- layout[["nc"]]
            byrow <- layout[["byrow"]]
        }
        op <- if(byrow)
            par(mfrow = c(nr, nc))
        else
            par(mfcol = c(nr, nc))
        on.exit(par(op))
    }

    attrs <- rep(attrs, length.out = length(x))
    for(i in seq_along(x)) {
        if(relation_is_acyclic(x[[i]]) && type[i] == "simplified") {
            ## possibly, transform relation to obtain a poset,
            ## and change main title if not specified.
            if (relation_is_linear_order(x[[i]])) {
                if (is.na(main[i]))
                    main[i] <- "Linear Order"
            } else if (relation_is_partial_order(x[[i]])) {
                if (is.na(main[i]))
                    main[i] <- "Partial Order"
            } else if(relation_is_weak_order(x[[i]])) {
                ## If x is a preference, use dual instead.
                x[[i]] <- dual(x[[i]])
                if (is.na(main[i]))
                    main[i] <- "Weak Order"
            } else { ## extract asymmetric part
                x[[i]] <- x[[i]] & dual(x[[i]])
                if (is.na(main[i]))
                    main[i] <- "Strict Preference Part"
            }
            if(is.null(attrs[[i]]$edge$arrowsize))
                attrs[[i]]$edge$arrowsize <- "0"
        } else {
            if(is.null(attrs[[i]]$edge$arrowsize))
                attrs[[i]]$edge$arrowsize <- "1"
        }

        if (type[i] == "simplified")
            x[[i]] <- transitive_reduction(x[[i]])

        ## Compute transitive reduction to avoid cluttered graph and
        ## extract incidence.
        I <- unclass(relation_incidence(x[[i]], limit = limit[i]))

        ## Perform reflexive reduction.
        diag(I) <- 0

        l <- if (is.null(labels) || is.null(labels[[i]]))
            labels(I)
        else
            labels[[i]]

        ## Transform to graphViz-compatible incidence.
        dimnames(I) <- lapply(l, .make_unique_labels)

        if(is.na(main[i]))
            main[i] <- ""

        Rgraphviz::plot(methods::as(I, "graphNEL"),
                        attrs = attrs[[i]], main = main[i], ...)
    }
}


####### plot incidence matrices

plot.relation_incidence <-
function(x, oma = c(3, 3, 1, 1), ...)
{
    o <- rev(order(colSums(x, na.rm = TRUE)))
    parhold <- par(no.readonly = TRUE)
    on.exit(par(parhold))
    par(las = 2L, oma = oma)
    x <- x[o, o]
    image.default(1:dim(x)[2L], 1:dim(x)[1L], t(x)[, dim(x)[1L]:1L],
                  axes = FALSE, col = gray(c(1, 0.5)),
                  xlab = "", ylab = "", ...)
    axis(1L, at = 1:dim(x)[1L], labels = rev(labels(x)[[1L]]))
    axis(2L, at = 1:dim(x)[2L], labels = rev(labels(x)[[2L]]))
}

### plot rankings

plot.ranking <-
function(x, ...)
    plot(as.relation(x), ...)


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***

