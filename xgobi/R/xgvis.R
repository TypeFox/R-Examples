xgvis <-
function(dmat	= NULL,
	 edges	= NULL,
	 pos	= NULL,
         rowlab     = colnames(dmat),
	 colors = NULL,
	 glyphs = NULL,
	 erase	= NULL,
	 lines	= NULL,
	 linecolors = NULL,
	 resources  = NULL,
	 display    = NULL,
         multi      = TRUE,
	 keep	= FALSE,
	 fprefix= "xgvis-")
{
    if (is.null(edges) && is.null(pos) && is.null(dmat))
        stop("One of dmat, edges, or pos must be present")

    basefile <- tempfile(fprefix)
    if(.Platform$OS.type == "windows") basefile <- gsub("\\\\", "/", basefile)

    tmpfiles <- character(0)
    if(!keep) on.exit(file.remove(tmpfiles))

    ## distance matrix ###
    if (!is.null(dmat)) {
        dmat <- eval(dmat)
        if (any(isinf <- is.infinite(dmat[!is.na(dmat)]))) {
            warning("xgvis can't handle Inf's in dmat; replaced with NA")
            dmat[isinf] <- NA
        }
        dfile <- paste(basefile, ".dist", sep="")
        write(t(dmat), file = dfile, ncolumns = ncol(dmat))
        tmpfiles <- c(tmpfiles, dfile)
    }

    ## Edges ###
    if (!is.null(edges)) {              # check data type
        if (!is.matrix(edges) || !is.numeric(edges) || dim(edges)[2] != 2)
            stop("The `edges' argument must be a numeric 2-column matrix")

        edgesfile <- paste(basefile, ".edges", sep="")
        if (nrow(edges) > 0) {
            write(t(edges), file = edgesfile, ncolumns=2)
        }
        tmpfiles <- c(tmpfiles, edgesfile)
    }

    ## position matrix ###
    if (!is.null(pos)) {
        pos <- eval(pos)
        if (any(isinf <- is.infinite(pos[!is.na(pos)]))) {
            warning("xgvis can't handle Inf's in pos; replaced with NA")
            pos[isinf] <- NA
        }
        pfile <- paste(basefile, ".pos", sep="")
        write(t(pos), file = pfile, ncolumns = ncol(pos))
        tmpfiles <- c(tmpfiles, pfile)
    }

    ## Row / Case labels ###
    if (!is.null(rowlab)) {
        if (!is.vector(rowlab) || !is.character(rowlab)) # check data type
            stop("The `rowlab' argument needs to be a character vector")
        if (!missing(rowlab) && length(rowlab) != NROW(dmat))
            stop("`rowlab' has wrong length (not matching NROW(dmat))")
        cat(rowlab, file = (rowfile <- paste(basefile, ".row", sep="")), sep="\n")
        tmpfiles <- c(tmpfiles, rowfile)
    }

    ## Colors ###
    if (!is.null(colors)) {
                                        # check data type
        if (!is.vector(colors) || !is.character(colors))
            stop("The `colors' argument needs to be a character vector")

        colorfile <- paste(basefile, ".colors", sep="")
        write(colors, file = colorfile, ncolumns=1)
        tmpfiles <- c(tmpfiles, colorfile)
    }

    ## Glyphs ###
    if (!is.null(glyphs)) {
                                        # check data type
        if (!is.vector(glyphs) || !is.numeric(glyphs))
            stop("The `glyphs' argument needs to be a numeric vector")

        glyphfile <- paste(basefile, ".glyphs", sep="")
        write(glyphs, file = glyphfile, ncolumns=1)
        tmpfiles <- c(tmpfiles, glyphfile)
    }

    ## Erase ###
    if (!is.null(erase)) {
                                        # check data type
        if (!is.vector(erase) || !is.numeric(erase))
            stop("The `erase' argument needs to be a numeric vector")

        erasefile <- paste(basefile, ".erase", sep="")
        write(erase, file = erasefile, ncolumns=1)
        tmpfiles <- c(tmpfiles, erasefile)
    }

    ## Connected lines ###
    if (!is.null(lines)) {
                                        # check data type
        if (!is.matrix(lines) || !is.numeric(lines) || dim(lines)[2] != 2)
            stop("The `lines' argument must be a numeric 2-column matrix")

        linesfile <- paste(basefile, ".lines", sep="")
        if (nrow(lines) > 0) {
            write(t(lines), file = linesfile, ncolumns=2)
            tmpfiles <- c(tmpfiles, linesfile)
        }
    }

    ## Line colors ###
    if ((!is.null(lines) || !is.null(edges)) && !is.null(linecolors)) {
                                        # check data type
        if (!is.vector(linecolors) || !is.character(linecolors))
            stop("The `linecolors' argument must be a character vector")

        linecolorfile <- paste(basefile, ".linecolors", sep="")
        write(linecolors, file = linecolorfile, ncolumns=1)
        tmpfiles <- c(tmpfiles, linecolorfile)
    }

    ## Resources ###
    if (!is.null(resources)) {
                                        # check data type
        if (!is.vector(resources) || !is.character(resources))
            stop("The `resources' argument must be a character vector")

        resourcefile <- paste(basefile, ".resources", sep="")
        write(resources, file = resourcefile, ncolumns=1)
        tmpfiles <- c(tmpfiles, resourcefile)
    }

    if(.Platform$OS.type == "windows") {
        xgpath <- system.file("scripts", "xgvis.bat", package="xgobi")
        command <- paste(xgpath, basefile)
        cat(command, "\n")
        s <- if(multi) system(command, wait=FALSE, minimized=TRUE)
        else system(command, intern=FALSE, wait=TRUE, minimized=TRUE)
    } else {
### Note to installer:
### Here you will need to specify the path to the xgvis executable
### on your system (here we assume it *is* in the user's PATH :)

        xgpath <- "xgvis"
        command <- paste(xgpath, basefile, if(multi) "&")
        cat(command, "\n")

        s <- system(command, FALSE)
    }
    ## Now wait a bit to allow gvis to start and read its files
    if(multi && !keep) Sys.sleep(3)
    invisible(s)
}
