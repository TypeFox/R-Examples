## These should really match the *brushColor[0-9]  `fallback resources' in
##  XGOBISRC/src/xgobitop.h :
## [these are ok for the "Dec. 1999" version of xgobi]:
xgobi.colors.default <-
  c("DeepPink", "OrangeRed1", "DarkOrange", "Gold", "Yellow",
    "DeepSkyBlue1", "SlateBlue1", "YellowGreen",
    "MediumSpringGreen", "MediumOrchid")

xgobi <-
function(matrx,
	 collab = colnames(matrx),
	 rowlab = rownames(matrx),
	 colors = NULL,
	 glyphs = NULL,
	 erase	= NULL,
	 lines	= NULL,
	 linecolors = NULL,
	 resources  = NULL,
	 title	= deparse(substitute(matrx)),
	 vgroups= NULL,
	 std	= "mmx",
	 nlinkable  = NULL,
	 subset = NULL,
	 display= NULL,
         multi  = TRUE,
	 keep	= FALSE,
	 fprefix= "xgobi-")
{
    x <- if(is.expression(matrx) || is.character(matrx))
	eval(matrx) else matrx
    if(is.data.frame(x)) x <- data.matrix(x)

    if (any(is.infinite(x[!is.na(x)])))
	stop("Sorry, xgobi can't handle Inf's")

    if (!is.null(title) && !is.character(title))
        stop("title must be a character string")
    dfile <- tempfile(paste(fprefix,
                            abbreviate(gsub("[^A-Za-z0-9]","",title), 5),
                            sep=""))

    if(.Platform$OS.type == "windows") dfile <- gsub("\\\\", "/", dfile)

    write.table(x, file = dfile, quote = FALSE,
		row.names = FALSE, col.names = FALSE)
    tmpfiles <- dfile
    if(!keep) on.exit(file.remove(tmpfiles))

    args <- paste("-std", std)

    ## Column / Var labels ###
    if (!is.null(collab)) {
	if (!is.vector(collab) || !is.character(collab))# check data type
	    stop("The `collab' argument needs to be a character vector")
	if (!missing(collab) && length(collab) != NCOL(x))
	    stop("`collab' has wrong length (not matching NCOL(x))")
        cat(collab, file = (colfile <- paste(dfile, ".col", sep="")), sep= "\n")
        tmpfiles <- c(tmpfiles, colfile)
    }
    ## Row / Case labels ###
    if (!is.null(rowlab)) {
	if (!is.vector(rowlab) || !is.character(rowlab))
	    stop("The `rowlab' argument needs to be a character vector")
	if (!missing(rowlab) && length(rowlab) != NROW(x))
	    stop("`rowlab' has wrong length (not matching NROW(x))")
        cat(rowlab, file = (rowfile <- paste(dfile, ".row", sep="")), sep="\n")
        tmpfiles <- c(tmpfiles, rowfile)
    }
    ## Variable groups ##
    if (!is.null(vgroups)) {
	if (!is.vector(vgroups) || !is.numeric(vgroups))
	    stop("The `vgroups' argument needs to be a numeric vector")
	cat(vgroups, file=(vgfile <- paste(dfile,".vgroups",sep="")), sep="\n")
        tmpfiles <- c(tmpfiles, vgfile)
    }
    ## Colors ##
    if (!is.null(colors)) {
	if (!is.vector(colors) || !is.character(colors))
	    stop("The `colors' argument needs to be a character vector")
	cat(colors, file = (clrfile <- paste(dfile,".colors",sep="")), sep="\n")
        tmpfiles <- c(tmpfiles, clrfile)
    }
    ## Glyphs ##
    if (!is.null(glyphs)) {
	if (!is.vector(glyphs) || !is.numeric(glyphs))
	    stop("The `glyphs' argument needs to be a numeric vector")
	glyphfile <- paste(dfile, ".glyphs", sep = "")
	cat(glyphs, file = glyphfile, sep = "\n")
        tmpfiles <- c(tmpfiles, glyphfile)
    }
    ## Erase ##
    if (!is.null(erase)) {
	if (!is.vector(erase) || !is.numeric(erase))
	    stop("The `erase' argument needs to be a numeric vector")
	erasefile <- paste(dfile, ".erase", sep = "")
	cat(erase, file = erasefile, sep = "\n")
        tmpfiles <- c(tmpfiles, erasefile)
    }
    ## Connected lines ##
    if (!is.null(lines)) {
	if (!is.matrix(lines) || !is.numeric(lines) || dim(lines)[2] != 2)
	    stop("The `lines' argument must be a numeric 2-column matrix")
	linesfile <- paste(dfile, ".lines", sep = "")
	unlink(linesfile)# in case it existed
	if (nrow(lines) > 0) {
	    for (i in 1:nrow(lines))
		cat(lines[i, ], "\n", file = linesfile, append = TRUE)
	}
        tmpfiles <- c(tmpfiles, linesfile)

	## Line colors ##
	if (!is.null(linecolors)) {
	    if (!is.vector(linecolors) || !is.character(linecolors))
		stop("The `linecolors' argument must be a character vector")
	    linecolorfile <- paste(dfile, ".linecolors", sep = "")
	    cat(linecolors, file = linecolorfile, sep = "\n")
            tmpfiles <- c(tmpfiles, linecolorfile)
	}
    }
    ## Resources ##
    if (!is.null(resources)) {
	if (!is.vector(resources) || !is.character(resources))
	    stop("The `resources' argument must be a character vector")
	resourcefile <- paste(dfile, ".resources", sep = "")
	cat(resources, file = resourcefile, sep = "\n")
        tmpfiles <- c(tmpfiles, resourcefile)
    }
    ## nlinkable ##
    if (!is.null(nlinkable)) {
	nlinkable <- as.integer(nlinkable)
	if (length(nlinkable) > 1)
	    stop("The `nlinkable' argument must be a scalar integer")
	linkablefile <- paste(dfile, ".nlinkable", sep = "")
	cat(nlinkable, "\n", file = linkablefile)
        tmpfiles <- c(tmpfiles, linkablefile)
    }
    ## subset ##
    subsetarg <- ""
    if (!is.null(subset)) {
	subset <- as.integer(subset)
	if (length(subset) > 1)
	    stop("The `subset' argument must be a scalar integer")
	if (subset == 0 || subset > nrow(x))
	    stop("The `subset' argument must be >0 and <= nrows")
	subsetarg <- paste(" -subset ", subset, sep = "")
	args <- paste(args, subsetarg, sep = " ")
    }

    if (!is.null(display)) {
	if (!is.character(display))
	    warning("display must be a character string")
	else args <- paste("-display", display, args)
    }

    if(.Platform$OS.type == "windows") {
        args <- paste("-vtitle", paste("'", title, "'", sep = ""), args)
        xgpath <- system.file("scripts", "xgobi.bat", package="xgobi")
        command <- paste(xgpath, args, dfile)
        cat(command, "\n")

        s <- if(multi) system(command, wait=FALSE, minimized=TRUE)
        else system(command, intern=FALSE, wait=TRUE, minimized=TRUE)
    } else {
        args <- paste("-title", paste("'", title, "'", sep = ""), args)
### Note to installer:
### Here you will need to specify the path to the xgobi executable
### on your system (here we assume it *is* in the user's PATH :)

        xgpath <- "xgobi"
        command <- paste(xgpath, args, dfile, if(multi) "&")
        cat(command, "\n")

        s <- system(command, FALSE)
    }
    ## Now wait a bit to allow xgobi to start and read its files
    if(multi && !keep) Sys.sleep(3)
    invisible(s)
}
