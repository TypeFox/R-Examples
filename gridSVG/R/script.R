
scriptGrob <- function(script=NULL, filename=NULL, type="application/ecmascript",
                       inline=FALSE, name=NULL) {
    body <- ""
    href <- ""
    if (!is.null(filename)) {
        if (inline) {
            body <- paste(readLines(filename), collapse="\n")
        } else {
            href <- filename
        }
    } else if (!is.null(script)) {
        body <- script
    } else {
        stop("No script specified")
    }
    sg <- grob(type = type, href = href, body = body,
               name = name, cl="script.grob")
    sg
}

grid.script <- function(...) {
    grid.draw(scriptGrob(...))
}

grobToDev.script.grob <- function(x, dev) {
    svgScript(x$body, x$href, x$type, x$name, svgdev=dev@dev)
}

