mcmcplot <- function(mcmcout, parms=NULL, regex=NULL, random=NULL, leaf.marker="[\\[_]", dir=tempdir(), filename="MCMCoutput", extension="html", title=NULL, heading=title, col=NULL, lty=1, xlim=NULL, ylim=NULL, style=c("gray", "plain"), greek=FALSE){
    ## This must come before mcmcout is evaluated in any other expression
    if (is.null(title))
        title <- paste("MCMC Plots: ", deparse(substitute(mcmcout)), sep="")
    if (is.null(heading))
        heading <- title

    style <- match.arg(style)

    ## Turn off graphics device if interrupted in the middle of plotting
    current.devices <- dev.list()
    on.exit( sapply(dev.list(), function(dev) if(!(dev %in% current.devices)) dev.off(dev)) )

    ## Convert input mcmcout to mcmc.list object
    mcmcout <- convert.mcmc.list(mcmcout)
    nchains <- length(mcmcout)
    if (is.null(col)){
        col <- mcmcplotsPalette(nchains)
    }
    css.file <- system.file("MCMCoutput.css", package="mcmcplots")
    css.file <- paste("file:///", css.file, sep="")
    htmlfile <- .html.begin(dir, filename, extension, title=title, cssfile=css.file)

    ## Select parameters for plotting
    if (is.null(varnames(mcmcout))){
        warning("Argument 'mcmcout' did not have valid variable names, so names have been created for you.")
        varnames(mcmcout) <- varnames(mcmcout, allow.null=FALSE)
    }
    parnames <- parms2plot(varnames(mcmcout), parms, regex, random, leaf.marker, do.unlist=FALSE)
    if (length(parnames)==0)
        stop("No parameters matched arguments 'parms' or 'regex'.")
    np <- length(unlist(parnames))

    cat('\n<div id="outer">\n', file=htmlfile, append=TRUE)
    cat('<h1>', heading, '</h1>', sep="", file=htmlfile, append=TRUE)
    cat('<div id="toc">\n', file=htmlfile, append=TRUE)
    cat('\n<h2>Table of Contents</h2>', file=htmlfile, append=TRUE)
    cat('<ul id="toc_items">\n', file=htmlfile, append=TRUE)
    for (group.name in names(parnames)) {
        cat(sprintf('<li class="toc_item"><a href="#%s">%s</a></li>\n', group.name, group.name), file=htmlfile, append=TRUE)
    }
    cat('</ul></div>\n', file=htmlfile, append=TRUE)

    cat('<div class="main">\n', file=htmlfile, append=TRUE)
    htmlwidth <- 640
    htmlheight <- 480
    for (group.name in names(parnames)) {
        cat(sprintf('<h2><a name="%s">Plots for %s</a></h2>\n', group.name, group.name), file=htmlfile, append=TRUE)
        for (p in parnames[[group.name]]) {
            pctdone <- round(100*match(p, unlist(parnames))/np)
            cat("\r", rep(" ", getOption("width")), sep="")
            cat("\rPreparing plots for ", group.name, ".  ", pctdone, "% complete.", sep="")
            gname <- paste(p, ".png", sep="")
            png(file.path(dir, gname), width=htmlwidth, height=htmlheight)
            plot_err <- tryCatch({
                mcmcplot1(mcmcout[, p, drop=FALSE], col=col, lty=lty, xlim=xlim, ylim=ylim, style=style, greek=greek)
              }, error=function(e) {e})
            dev.off()
            if (inherits(plot_err, "error")) {
                cat(sprintf('<p class="plot_err">%s. %s</p>', p, plot_err),
                    file=htmlfile, append=TRUE)
            } else {
                .html.img(file=htmlfile, class="mcmcplot", src=gname,
                  width=htmlwidth, height=htmlheight)
            }
        }
    }
    cat("\r", rep(" ", getOption("width")), "\r", sep="")
    cat('\n</div>\n</div>\n', file=htmlfile, append=TRUE)
    .html.end(htmlfile)
    full.name.path <- paste("file://", htmlfile, sep="")
    browseURL(full.name.path)
    invisible(full.name.path)
}
