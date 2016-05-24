
# FIXME:  What should happen if a grob has BOTH group and individual hrefs?
#         Is that an error?

hyperlinkGrob <- function(x, href, show=NULL, group=TRUE) {
    if (group) 
        x$groupLinks <- href
    else
        x$links <- href

    # Determines which window the link is going to open in
    x$show <- show

    class(x) <- unique(c("linked.grob", class(x)))
    x
}

grid.hyperlink <- function(path, href, show=NULL, group=TRUE, redraw=FALSE,
                           strict=FALSE, grep=FALSE, global=FALSE) {
    grobApply(path, function(path) {
        grid.set(path, hyperlinkGrob(grid.get(path), href, show, group),
                 redraw = redraw)
    }, strict = strict, grep = grep, global = global)
    invisible()
}

link <- function(x) {
    UseMethod("link")
}

link.grob <- function(x) {
    x$name <- getID(x$name, "grob", FALSE)
    href <- x$links
    if (!is.null(href)) {
        n <- length(href)
        if (is.null(names(href)))
            names(href) <- subGrobName(x$name, 1:n)
    }
    groupHref <- x$groupLinks
    if (!is.null(groupHref))
        names(groupHref) <- x$name
    c(href, groupHref)
}

# A hopefully useful default for gTrees
link.gTree <- function(x, ...) {
    x$name <- getID(x$name, "grob", FALSE)
    href <- x$links
    if (!is.null(href)) {
        n <- length(href)
        if (is.null(names(href)))
            names(href) <- sapply((x$childrenOrder)[1:n],
                                  function(x) getID(x, "grob", FALSE))
    }
    groupHref <- x$groupLinks
    if (!is.null(groupHref))
        names(groupHref) <- x$name
    c(href, groupHref)
}

linkShow <- function(x) {
    UseMethod("linkShow")
}

linkShow.grob <- function(x, ...) {
    x$name <- getID(x$name, "grob", FALSE)
    show <- x$show
    if (is.null(show))
        return("")
    if (!is.null(x$links)) {
        n <- length(show)
        if (is.null(names(show)))
            names(show) <- subGrobName(x$name, 1:n)
    }
    if (!is.null(x$groupLinks))
        names(show) <- x$name
    show
}

linkShow.gTree <- function(x, ...) {
    x$name <- getID(x$name, "grob", FALSE)
    show <- x$show
    if (is.null(show))
        return("")
    if (!is.null(x$links)) {
        n <- length(show)
        if (is.null(names(show)))
            names(show) <- sapply((x$childrenOrder)[1:n],
                                  function(x) getID(x, "grob", FALSE))
    }
    if (!is.null(x$groupLinks))
        names(show) <- x$name
    show
}

# Set the 'links' slot in the device
# The catsvg() function in svg.R picks this up
# and matches links to element names
primToDev.linked.grob <- function(x, dev) {
    dev@links <- link(x)
    dev@show <- linkShow(x)
    NextMethod()
}

# gridToDev method for linked.grob objects
# grobToDev.linked.grob <- function(x, dev) {
#   svgStartLink(x$href, dev@dev)
#   NextMethod()
#   svgEndLink(dev@dev)
# }

# Ensure the hyperlink is retained on a forced grob
forceGrob.linked.grob <- function(x) {
    y <- NextMethod()
    if (inherits(y, "forcedgrob")) {
        y$links <- x$links
        y$groupLinks <- x$groupLinks
        y$show <- x$show
        class(y) <- unique(c("linked.grob", class(y)))
    }
    y
}
