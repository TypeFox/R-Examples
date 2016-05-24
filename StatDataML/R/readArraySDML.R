readArraySDML <- function(x)
{
    if (is.null(x)) return(NULL)

    ## parse dimension
    dimension <- readDimensionSDML(x[["dimension"]])

    ## parse categories
    type <- readType(x[["type"]])

    ## parse properties
    attrib <- readProperties(x[["properties"]])

    info <- NULL

    if (!is.null(x[["data"]]) || !is.null(x[["textdata"]])) {

        if(!is.null(x[["data"]])) {
            attribs <- getAttrSDML(x[["data"]])
            vals <- if (!is.null(type$type) && type$type == "numeric" &&
                        type$mode == "complex")
                getComplexDataSDML(xmlChildren(x[["data"]]))
            else
                getDataSDML(xmlChildren(x[["data"]]))
            info <- attr(vals, "info")
        } else {
            attribs <- getAttrSDML(x[["textdata"]])
            vals <- if (is.null(xmlChildren(x[["textdata"]])))
                character()
            else
                getTextDataSDML(xmlChildren(x[["textdata"]])[[1]],
                                attribs, type$type)
        }

        if (length(vals) != prod(dimension$dim))
            stop(paste("Wrong dimension !",
                       paste(dimension$dim, collapse=", "),
                       length(vals), collapse= " "))

        if (length(dimension$dim) > 1) {
            vals <- array(vals, dim=dimension$dim, dimnames=dimension$names)
            if (!is.null(x[["data"]])) {
                if (!is.null(info))
                    info <- array(info, dim = dimension$dim,
                                  dimnames = dimension$names)
            }
        } else {
            if (!length(vals))
                vals <- vector()
            else
                if (length(dimension$names[[1]]))
                    names(vals) <- dimension$names[[1]]
        }

        ## handle types
        if (!is.null(type$type)) {

            if (!is.null(x[["textdata"]]) && type$type == "logical") {
                tr <- as.character(default(attribs, "true", "1"))
                fa <- as.character(default(attribs, "false", "0"))
                vals1 <- vals == tr
                vals2 <- vals != fa
                vals <- vals1
                vals[vals != vals2] <- NA
            }

            if (type$type == "numeric") {
                if (type$mode %in% c("integer", "complex"))
                    mode(vals) <- type$mode
                else ## default mode if none
                    mode(vals) <- "double"

                if (type$mode %in% c("integer", "real") &&
                    (min(vals, na.rm=TRUE) < type$min ||
                     max(vals, na.rm=TRUE) > type$max))
                    warning("numeric values out of specified range.", call. = FALSE)
            }

            if (type$type == "categorical") {
                vals <- structure(as.integer(vals), levels = type$labels)
                class(vals) <- if (type$mode == "ordered")
                    c("ordered", "factor")
                else
                    "factor"
            }

            if (type$type == "datetime")
                vals <- as.POSIXct(strptime(vals, format="%Y-%m-%dT%H:%M:%S"))

            if (type$type == "character")
                vals[] <- as.character(vals)
        }

        atvals <- attributes(vals)
        if (!is.null(atvals)) attrib <- c(atvals, attrib)

        ## recombine possibly splitted class attribute
        ind <- names(attrib) == "class"
        if (length(ind)) {
            newclass <- as.character(unlist(attrib[ind]))
            attrib[ind] <- NULL
            attrib[["class"]] <- newclass
        }

        if (!is.null(attrib)) attributes(vals) <- attrib
        if (!is.null(info))
            attr(vals, "info") <- info
        return(vals)
    }
}

getComplexDataSDML <- function(y) {
    y <- y[which(sapply(y[seq_along(y)], xmlName) != "text")]
    ret <- as.complex(sapply(y, function(x) {
        if (x$name=="na") return(NA)
        cs <- getDataSDML(xmlChildren(x))
        complex(1, as.double(cs[1]), as.double(cs[2]))
    }))
    i <- sapply(y, function(x) if (is.null(a <-
                                           getAttrSDML(x)[["info"]])) NA else a)
    attributes(i) <- NULL
    if (!all(is.na(i))) attr(ret, "info") <- i
    ret
}

getDataSDML <- function(y)
{
    y <- y[which(sapply(y[seq_along(y)], xmlName) != "text")]
    w <- sapply(y,
                function(x) switch (x$name,
                                    na = NA,
                                    T = TRUE,
                                    F = FALSE,
                                    min = , max = switch(x[[1]]$name,
                                           posinf = +Inf,
                                           neginf = -Inf,
                                           xmlValue(x[[1]])
                                           ),
                                    e = if (!length(x)) "" else
                                    switch(x[[1]]$name,
                                           posinf = +Inf,
                                           neginf = -Inf,
                                           nan    = NaN,
                                           xmlValue(x[[1]])
                                           ),
                                    NULL
                                    )
                )
    attributes(w) <- NULL

    i <- sapply(y, function(x) if (is.null(a <- getAttrSDML(x)[["info"]]))
                NA else a)
    attributes(i) <- NULL

    if(is.character(w)) {
        w <- gsub("&amp;", "&", w)
        w <- gsub("&lt;", "<", w)
        w <- gsub("&gt;", ">", w)
    }

    if (!all(is.na(i))) attr(w, "info") <- i
    w
}

getTextDataSDML <- function(y, attribs, type)
{
    y <- gsub("&amp;", "&", y)
    y <- gsub("&lt;", "<", y)
    y <- gsub("&gt;", ">", y)
    y[y == attribs[["na.string"]]] <- NA
    if (type == "character")
        y[y == attribs[["null.string"]]] <- ""
    if (type == "numeric") {
        y[grep(attribs[["neginf.string"]], y, fixed = TRUE)] <- "-Inf"
        y[grep(attribs[["posinf.string"]], y, fixed = TRUE)] <- "Inf"
        y[grep(attribs[["nan.string"]], y, fixed = TRUE)] <- "NaN"
    }
    if (type == "logical") {
        y[y == attribs[["true"]]] <- "1"
        y[y == attribs[["false"]]] <- "0"
    }
    y
}

default <- function (attr, name, defval)
    if (name %in% names(attr)) attr[name] else defval




