# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# phonR version 1.0-3
# Functions for phoneticians and phonologists
# AUTHOR: Daniel McCloy, drmccloy@uw.edu
# LICENSED UNDER THE GNU GENERAL PUBLIC LICENSE v3.0:
# http://www.gnu.org/licenses/gpl.html
# DEVELOPMENT OF THIS PACKAGE WAS FUNDED IN PART BY NIH-R01DC006014
#
# CHANGELOG:
# v1.0: major refactor of the entire codebase. Enhancements: more control
# over color and style. Support for most standard plot/par arguments (xlim,
# ylim, cex.axis, etc). New drawing functions for repulsive force heatmaps
# and convex hulls. Hulls, polygons, and ellipses can have color fills.
# Shortcut argument "pretty=TRUE" soothes the senses. Diphthongs with
# arrowheads.
#
# v0.4: bugfixes: poly.order now works with arbitrary labels; bug in
# s-centroid calculation fixed.  Enhancements: added user-override
# arguments for color, shape and linestyle; added support for diphthong
# plotting, argument poly.include eliminated (inferred from elements
# present in poly.order), new argument points.label allows override of
# points label when points='text'.
#
# v0.3 bugfixes: font specification on windows now works for direct-to-
# file output. Enhancements: graphics handling overhauled to use base
# graphics instead of Cairo(). Several new output formats added. Raster
# resolution and font size now specifiable. Improved error handling.
#
# v0.2 bugfixes: points.alpha and means.alpha now work for grayscale
# plots. Plots with polygons or ellipses but no shapes now get proper
# legend type (lines, not boxes). Graphical parameters now captured and
# restored when plotting to onscreen device. Vowels with no variance
# (e.g., single tokens) no longer crash ellipse function. Vowels not in
# default poly.order() no longer go unplotted when points='text'.
# Enhancements: support for custom axis titles (to accommodate pre-
# normalized values), point and mean sizes, and fonts. Custom line types
# added (11 total now).
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#' @export
plotVowels <- function(f1, f2, vowel=NULL, group=NULL,
    plot.tokens=TRUE, pch.tokens=NULL, cex.tokens=NULL, alpha.tokens=NULL,
    plot.means=FALSE, pch.means=NULL, cex.means=NULL, alpha.means=NULL,
    hull.line=FALSE, hull.fill=FALSE, hull.args=NULL,
    poly.line=FALSE, poly.fill=FALSE, poly.args=NULL, poly.order=NA,
    ellipse.line=FALSE, ellipse.fill=FALSE, ellipse.conf=0.6827, ellipse.args=NULL,
    diph.arrows=FALSE, diph.args.tokens=NULL, diph.args.means=NULL,
    diph.label.first.only=TRUE, diph.mean.timept=1, diph.smooth=FALSE,
    heatmap=FALSE, heatmap.args=NULL, heatmap.legend=FALSE, heatmap.legend.args=NULL,
    var.col.by=NULL, var.style.by=NULL, fill.opacity=0.3, label.las=NULL,
    legend.kwd=NULL, legend.args=NULL, pretty=FALSE, output='screen', ...)
{
    # # # # # # # # # # #
    # HANDLE EXTRA ARGS #
    # # # # # # # # # # #
    exargs <- list(...)
    font.specified <- 'family' %in% names(exargs)
    # two arguments get overridden no matter what
    exargs$ann <- FALSE
    exargs$type <- 'n'
    # Some graphical devices only support inches, so we convert here.
    if ("units" %in% names(exargs)) {
        if (!exargs$units %in% c("in", "cm", "mm", "px")) {
            warning("Unsupported argument value '", units, "': 'units' must be ",
                    "one of 'in', 'cm', 'mm', or 'px'. Using default ('in').")
            exargs$units <- "in"
        }
        if (output %in% c("pdf", "svg", "screen")) {
            if ("width" %in% names(exargs)) {
                if      (exargs$units == "cm") exargs$width <- exargs$width/2.54
                else if (exargs$units == "mm") exargs$width <- exargs$width/2540
                else if (exargs$units == "px") exargs$width <- exargs$width/72
            }
            if ("height" %in% names(exargs)) {
                if      (exargs$units == "cm") exargs$height <- exargs$height/2.54
                else if (exargs$units == "mm") exargs$height <- exargs$height/2540
                else if (exargs$units == "px") exargs$height <- exargs$height/72
            }
        }
    }
    # Some args are only settable by direct par() call (not via plot(), etc).
    # Not strictly true for "family" or "las", but works better this way.
    par.only <- c("ask", "fig", "fin", "las", "lheight", "mai", "mar",
                "mex", "mfcol", "mfrow", "mfg", "new", "oma", "omd", "omi",
                "pin", "plt", "ps", "pty", "usr", "xlog", "ylog", "ylbias")
    if (output == "screen") {
        par.only <- append(par.only, "family")
    } else {
        file.only <- c("filename", "width", "height", "units", "pointsize",
                    "res", "quality", "compression", "family")
        file.args <- exargs[names(exargs) %in% file.only]
        exargs <- exargs[!(names(exargs) %in% file.only)]
    }
    par.args <- exargs[names(exargs) %in% par.only]
    # allow separate "las" for axis numbers & axis labels
    if (is.null(label.las)) {
        if ("las" %in% names(par.args)) label.las <- par.args$las
        else                            label.las <- par("las")
    }
    # split out arguments for annotations. If "pretty", axes get drawn separately
    # from plot(); other annotation gets drawn in separate mtext() calls whether
    # "pretty" or not.
    if (pretty) {
        axis.only <- c("cex.axis", "col.axis", "font.axis")
        axis.args <- exargs[names(exargs) %in% axis.only]
    } else {
        axis.only <- c()
    }
    main.only <- c("cex.main", "col.main", "font.main")
    sub.only <- c("cex.sub",  "col.sub",  "font.sub")
    lab.only <- c("cex.lab",  "col.lab",  "font.lab")
    main.args <- exargs[names(exargs) %in% main.only]
    sub.args <- exargs[names(exargs) %in% sub.only]
    lab.args <- exargs[names(exargs) %in% lab.only]
    names(main.args) <- gsub(".main", "", names(main.args))
    names(sub.args) <- gsub(".sub", "", names(sub.args))
    names(lab.args) <- gsub(".lab", "", names(lab.args))
    if ("xlab" %in% names(exargs)) xlab <- exargs$xlab
    else                           xlab <- "F2"
    if ("ylab" %in% names(exargs)) ylab <- exargs$ylab
    else                           ylab <- "F1"
    if ("main" %in% names(exargs)) main <- exargs$main
    #else if (pretty)              main <- "Vowels"
    else                           main <- ""
    if ("sub" %in% names(exargs))   sub <- exargs$sub
    #else if (pretty)               sub <- "plotted with love using phonR"
    else                            sub <- ""
    exargs$xlab <- NULL
    exargs$ylab <- NULL
    exargs$main <- NULL
    exargs$sub <- NULL
    # add "las" into label args
    lab.args <- as.list(c(lab.args, las=label.las))
    # don't pass args twice!
    args.to.remove <- c(par.only, main.only, axis.only, lab.only, sub.only)
    exargs <- exargs[!names(exargs) %in% args.to.remove]

    # # # # # # # # # #
    # OUTPUT PARSING  #
    # # # # # # # # # #
    output <- tolower(output)
    if (output=="jpeg") output <- "jpg"
    if (output=="tiff") output <- "tif"
    output.types <- c("pdf", "svg", "jpg", "tif", "png", "bmp", "screen")
    output.raster <- c("jpg", "tif", "png", "bmp", "screen")
    if (!(output %in% output.types)) {
        warning("Unknown argument value '", output, "': 'output' ",
                "must be one of 'pdf', 'svg', 'png', 'tif', 'bmp', ",
                "'jpg', or 'screen'. Using default ('screen').")
        output <- "screen"
    }

    # # # # # # # # # # #
    # LEGEND KWD CHECK  #
    # # # # # # # # # # #
    legend.kwds <- c("left", "right", "top", "bottom", "center", "topleft",
                    "topright", "bottomleft", "bottomright")
    if (!is.null(legend.kwd)) {
        if (!legend.kwd %in% legend.kwds) {
            stop(paste(c("legend.kwd must be one of '",
                         paste(legend.kwds, collapse="', '"), "'."),
                 collapse=""))
    }   }

    # # # # # # # # # # # # # # # # # #
    # PRELIMINARY DIPHTHONG HANDLING  #
    # # # # # # # # # # # # # # # # # #
    if (is.vector(f1) && is.vector(f2)) {
        diphthong <- FALSE
    } else if (length(dim(f1)) == 1 && length(dim(f2)) == 1) {
        f1 <- as.vector(f1)
        f2 <- as.vector(f2)
        diphthong <- FALSE
    } else {
        if (!all(dim(f1) == dim(f2))) {
            stop("Unequal dimensions for 'f1' and 'f2'.")
        } else if (length(dim(f1)) > 2) {
            stop("Argument 'f1' has more than two dimensions.")
        } else if (length(dim(f2)) > 2) {
            stop("Argument 'f2' has more than two dimensions.")
        } else if (!is.null(vowel) && length(vowel) != dim(f1)[1]) {
            stop("First axis of 'f1' does not equal length of 'vowel'.")
        }
        diphthong <- TRUE
    }
    if (!diphthong) {
        if (length(f2) != length(f1)) stop('Unequal dimensions for "f1" and "f2".')
        else if (!is.null(vowel) && length(vowel) != length(f1)) {
            stop("Unequal dimensions for 'f1' and 'vowel'.")
        }
        l <- length(f1)
    } else {
        f1d <- f1
        f2d <- f2
        f1 <- f1d[,diph.mean.timept]
        f2 <- f2d[,diph.mean.timept]
        l <- nrow(f1d)
    }

    # # # # # # # # # # # # # #
    # DIPHTHONG ARG HANDLING  #
    # # # # # # # # # # # # # #
    if (diphthong) {
        arrow.only <- c("angle", "length")
        line.only <- c("type")
        if (pretty) {
            if (is.null(pch.means)) type <- "l"
            else                    type <- "o"
            pretty.diph.tokens <- list(length=0.1, angle=20, type=type)
            pretty.diph.means  <- list(length=0.1, angle=20, type=type,
                                       lwd=2*par("lwd"))
            # user override
            pretty.diph.tokens[names(diph.args.tokens)] <- diph.args.tokens
            pretty.diph.means[names(diph.args.means)] <- diph.args.means
            # re-unify
            diph.args.tokens <- pretty.diph.tokens
            diph.args.means  <- pretty.diph.means
        }
        if (diph.arrows) {
            diph.arrow.tokens <- diph.args.tokens[!names(diph.args.tokens) %in% line.only]
            diph.arrow.means  <- diph.args.means[!names(diph.args.means) %in% line.only]
        }
        diph.args.tokens[names(diph.args.tokens) %in% arrow.only] <- NULL
        diph.args.means[names(diph.args.means) %in% arrow.only] <- NULL
    }

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # PRELIMINARY HANDLING OF GROUPING FACTOR, COLOR, AND STYLE #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    if (is.null(vowel)) v <- rep(NA, l)
    else                v <- factor(vowel, levels=unique(vowel))
    if (is.null(group)) gf <- rep('gf', l)
    else 			    gf <- factor(group, levels=unique(group))
    # used later to set default polygon color when color varies by vowel
    if (identical(as.numeric(factor(var.col.by)), 
    			  as.numeric(factor(vowel)))) col.by.vowel <- TRUE
    else                                      col.by.vowel <- FALSE
    # var.col.by & var.style.by
    if (!is.null(var.col.by[1])) {
        if (is.na(var.col.by[1])) legend.col.lab <- NULL
        else                      legend.col.lab <- unique(as.character(var.col.by))
        var.col.by <- as.numeric(factor(var.col.by, levels=unique(var.col.by)))
    } else {
        legend.col.lab <- c()
        var.col.by <- rep(1, l)  # default to black
    }
    if (!is.null(var.style.by[1])) {
        if (is.na(var.style.by[1])) legend.style.lab <- NULL
        else                        legend.style.lab <- unique(as.character(var.style.by))
        var.style.by <- as.numeric(factor(var.style.by, levels=unique(var.style.by)))
    } else {
        legend.style.lab <- c()
        var.style.by <- rep(1, l)  # default to solid
    }
    num.col <- length(unique(var.col.by))
    num.sty <- length(unique(var.style.by))
    # misc. plotting defaults
    if (is.null(cex.tokens))    cex.tokens <- par('cex')
    if (is.null(cex.means))      cex.means <- par('cex')

    # # # # # # # # # # # # #
    # DEFAULTS FOR "PRETTY" #
    # # # # # # # # # # # # #
    if (pretty) {
        # still suppress axes if axes=FALSE
        draw.axes <- TRUE
        if ("axes" %in% names(exargs)) {
            if (!exargs$axes) draw.axes <- FALSE
        }
        # if no colors specified, use equally spaced HCL values
        # [-1] avoids duplicate hues 0 and 360
        if (num.col == 1) pretty.col <- "black"
        else {
            hue <- seq(0,  360, length.out=1+num.col)[-1]
            chr <- seq(60, 100, length.out=num.col)
            lum <- seq(60,  40, length.out=num.col)
            pretty.col <- hcl(hue, chr, lum, alpha=1)
        }
        # PCH: filled / open {circ,tri,squ,diam}, plus, x, inverted open tri
        pretty.pch <- rep(c(16,1,17,2,15,0,18,5,3,4,6), length.out=num.sty)
        # LTY: custom linetypes more readily distinguishable
        pretty.lty <- c('solid', '44', 'F4', '4313', 'F3131313', '23F3',
                        '232923', '23258385', '282823B3', '13', '82')
        pretty.args <- list(mgp=c(2,0.5,0), xaxs='i', yaxs='i', axes=FALSE,
                            fg=hcl(0,0,40), tcl=-0.25, xpd=NA,
                            pch=pretty.pch, lty=pretty.lty, col=pretty.col)
        pretty.par.args <- list(mar=c(1,1,5,5), las=1)
        # legend args
        pretty.legend.args <- list(bty="n", seg.len=1)
        # LET USER-SPECIFIED ARGS OVERRIDE "PRETTY" DEFAULTS
        pretty.args[names(exargs)] <- exargs
        pretty.par.args[names(par.args)] <- par.args
        pretty.legend.args[names(legend.args)] <- legend.args
        # RE-UNIFY TO AVOID LATER LOGIC BRANCHING
        exargs <- pretty.args
        par.args <- pretty.par.args
        legend.args <- pretty.legend.args
    }

    # # # # # # # # # #
    # OTHER DEFAULTS  #
    # # # # # # # # # #
    vary.col <- ifelse(is.na(var.col.by[1]), FALSE, TRUE)
    vary.sty <- ifelse(is.na(var.style.by[1]), FALSE, TRUE)
    # colors: use default pallete if none specified and pretty=FALSE
    if (!'col' %in% names(exargs)) exargs$col <- palette()
    if (vary.col) exargs$col <- exargs$col[var.col.by]
    # linetypes & plotting characters
    if (!'lty' %in% names(exargs)) exargs$lty <- seq_len(num.sty)
    if (!'pch' %in% names(exargs)) exargs$pch <- seq_len(num.sty)
    if (vary.sty) {
        exargs$lty <- exargs$lty[var.style.by]
        exargs$pch <- exargs$pch[var.style.by]
    }
    if (is.null(pch.tokens)) pcht <- exargs$pch
    else                     pcht <- pch.tokens
    if (is.null(pch.means))  pchm <- exargs$pch
    else                     pchm <- pch.means
    # transparency
    trans.col <- makeTransparent(exargs$col, fill.opacity)
    trans.fg <- makeTransparent(par('fg'), fill.opacity)
    # ellipse colors
    if (ellipse.line) {
        if ("border" %in% names(ellipse.args)) {
            ellipse.args$border <- rep(ellipse.args$border, length.out=num.col)
            if (vary.col) ellipse.line.col <- ellipse.args$border[var.col.by]
            else          ellipse.line.col <- ellipse.args$border
        } else            ellipse.line.col <- exargs$col
    } else                ellipse.line.col <- NA[var.col.by]
    if (ellipse.fill) {
        if ("col" %in% names(ellipse.args)) {
            if (vary.col) ellipse.fill.col <- ellipse.args$col[var.col.by]
            else          ellipse.fill.col <- ellipse.args$col
        } else            ellipse.fill.col <- trans.col
    } else                ellipse.fill.col <- NA[var.col.by]
    # hull colors
    if (hull.line) {
        if ("border" %in% names(hull.args)) {
            if (vary.col)        hull.line.col <- hull.args$border[var.col.by]
            else                 hull.line.col <- hull.args$border
        } else if (col.by.vowel) hull.line.col <- par("fg")
        else                     hull.line.col <- exargs$col
    } else                       hull.line.col <- NA[var.col.by]
    if (hull.fill) {
        if ("col" %in% names(hull.args)) {
            if (vary.col)        hull.fill.col <- hull.args$col[var.col.by]
            else                 hull.fill.col <- hull.args$col
        } else if (col.by.vowel) hull.fill.col <- trans.fg
        else                     hull.fill.col <- trans.col
    } else                       hull.fill.col <- NA[var.col.by]
    # polygon colors
    if (poly.line) {
        if ("border" %in% names(poly.args)) {
            if (vary.col)        poly.line.col <- poly.args$border[var.col.by]
            else                 poly.line.col <- poly.args$border
        } else if (col.by.vowel) poly.line.col <- par("fg")
        else                     poly.line.col <- exargs$col
    } else                       poly.line.col <- NA[var.col.by]
    if (poly.fill) {
        if ("col" %in% names(poly.args)) {
            if (vary.col)        poly.fill.col <- poly.args$col[var.col.by]
            else                 poly.fill.col <- poly.args$col
        } else if (col.by.vowel) poly.fill.col <- trans.fg
        else                     poly.fill.col <- trans.col
    } else                       poly.fill.col <- NA[var.col.by]

    # # # # # # # # # # # # # # #
    # TOKEN / MEAN TRANSPARENCY #
    # # # # # # # # # # # # # # #
    if (!is.null(alpha.tokens)) col.tokens <- makeTransparent(exargs$col, alpha.tokens)
    else                        col.tokens <- exargs$col
    if (!is.null(alpha.means))   col.means <- makeTransparent(exargs$col, alpha.means)
    else                         col.means <- exargs$col

    # # # # # # # # # # # # # # #
    # INITIALIZE OUTPUT DEVICES #
    # # # # # # # # # # # # # # #
    if      (output=='pdf') do.call(cairo_pdf, file.args)
    else if (output=='svg') do.call(svg, file.args)
    else if (output=='jpg') do.call(jpeg, file.args)
    else if (output=='tif') do.call(tiff, file.args)
    else if (output=='png') do.call(png, file.args)
    else if (output=='bmp') do.call(bmp, file.args)
    # FONT HANDLING FOR WINDOWS
    is.win <- .Platform$OS.type == 'windows'
    if (is.win && font.specified && output %in% output.raster) {
        windowsFonts(phonr=windowsFont(par.args$family))
        par.args$family <- 'phonr'
        #if (output=='screen') warning("Font specification may fail if saving ",
        #                            "as PDF from onscreen plot window menu. ",
        #                            "To ensure PDF font fidelity, run ",
        #                            "plotVowels() with output='pdf'.")
    }
    # INITIAL CALL TO PAR()
    op <- par(par.args)

    # # # # # # # # # # # # # #
    # COLLECT INTO DATAFRAMES #
    # # # # # # # # # # # # # #
    d <- data.frame(f1=f1, f2=f2, v=v, gf=factor(gf, levels=unique(gf)),
                    col.tokens=col.tokens, col.means=col.means,
                    style=var.style.by,
                    ellipse.fill.col=ellipse.fill.col,
                    ellipse.line.col=ellipse.line.col,
                    poly.fill.col=poly.fill.col,
                    poly.line.col=poly.line.col,
                    hull.fill.col=hull.fill.col,
                    hull.line.col=hull.line.col,
                    pchmeans=pchm, pchtokens=pcht,
                    stringsAsFactors=FALSE)
    if (diphthong) {
        d$f2d <- lapply(apply(f2d, 1, list), unlist, recursive=TRUE, use.names=FALSE)
        d$f1d <- lapply(apply(f1d, 1, list), unlist, recursive=TRUE, use.names=FALSE)
    }
    if (is.null(vowel) && is.null(group))  byd <- list(d=d)
    else                                   byd <- by(d, d[c('v','gf')], identity)
    # DATAFRAME OF MEANS
    m <- lapply(byd, function(i) {
        if (!is.null(i)) {
            with(i, data.frame(f2=mean(f2, na.rm=TRUE),
                               f1=mean(f1, na.rm=TRUE),
                               v=unique(v), gf=unique(gf),
                               col.means=ifelse(length(unique(col.means)) == 1,
                               					unique(col.means), par("fg")),
                               style=ifelse(length(unique(style)) == 1,
                               				unique(style), 1),
                               poly.fill.col=unique(poly.fill.col),
                               poly.line.col=unique(poly.line.col),
                               ellipse.fill.col=unique(ellipse.fill.col),
                               ellipse.line.col=unique(ellipse.line.col),
                               pchmeans=unique(pchmeans),
                               stringsAsFactors=FALSE))
    }})
    m <- do.call(rbind, m)
    m$gfn <- as.numeric(factor(m$gf, levels=unique(m$gf)))
    if (diphthong) {
        m$f2d <- do.call(rbind, lapply(byd, function(i) if (!is.null(i))
                         with(i, list(colMeans(do.call(rbind, f2d))))))
        m$f1d <- do.call(rbind, lapply(byd, function(i) if (!is.null(i))
                         with(i, list(colMeans(do.call(rbind, f1d))))))
    }
    # MEANS & COVARIANCES FOR ELLIPSE DRAWING
    if (ellipse.fill || ellipse.line) {
        mu <- lapply(byd, function(i) { if (!is.null(i)) {
            with(i, list(colMeans(cbind(f2, f1), na.rm=TRUE)))
        }})
        mu <- do.call(rbind, mu)
        sigma <- lapply(byd, function(i) { if (!(is.null(i))) {
            with(i[!(is.na(i$f2)) && !(is.na(i$f1)),],
                 list(cov(cbind(f2, f1))))
        }})
        # the covariance calculation above still may yield some NA cov. matrices,
        # due to some vowels having only 1 token. This is handled later.
        sigma <- do.call(rbind, sigma)
        m$mu <- mu
        m$sigma <- sigma
    }

    # # # # # # # # # # # # #
    # DETERMINE PLOT BOUNDS #
    # # # # # # # # # # # # #
    # TOKEN EXTREMA
    if (diphthong) plot.bounds <- cbind(range(f2d, finite=TRUE),
                                        range(f1d, finite=TRUE))
    else           plot.bounds <- apply(d[,c('f2','f1')], 2, range, finite=TRUE)
    # ELLIPSE EXTREMA
    if (ellipse.fill || ellipse.line) {
        ellipse.param <- apply(m, 1, function(i) {
            if (any(is.na(i$sigma))) {
                i$sigma <- matrix(rep(0, 4), nrow=2)
                msg <- ifelse(i$gf == "gf", as.character(i$v),
                              paste("(", i$gf, ", ", i$v, ")", sep=""))
                message("No ellipse drawn for ", msg, " because there is only one token.")
            }
            list('mu'=i$mu, 'sigma'=i$sigma, 'alpha'=1 - ellipse.conf, 'draw'=FALSE)
        })
        ellipse.points <- lapply(ellipse.param,
                                 function(i) do.call(ellipse, i))
        ellipse.bounds <- lapply(ellipse.points,
                                 function(i) apply(i, 2, range, finite=TRUE))
        ellipse.bounds <- apply(do.call(rbind, ellipse.bounds), 2,
                                range, finite=TRUE)
        plot.bounds <- apply(rbind(plot.bounds, ellipse.bounds), 2,
                             range, finite=TRUE)
    }

    # # # # # # # # # # #
    # PREPARE GARNISHES #
    # # # # # # # # # # #
    # decide plot bounds
    user.set.xlim <- "xlim" %in% names(exargs)
    user.set.ylim <- "ylim" %in% names(exargs)
    if (!user.set.xlim) exargs$xlim <- rev(plot.bounds[,1])
    if (!user.set.ylim) exargs$ylim <- rev(plot.bounds[,2])
    # calculate nice tickmark intervals
    xticks <- prettyTicks(exargs$xlim)
    yticks <- prettyTicks(exargs$ylim)
    # ensure axes begin and end at a tickmark (unless xlim/ylim overtly set)
    if (pretty && !user.set.xlim) exargs$xlim <- rev(range(xticks))
    if (pretty && !user.set.ylim) exargs$ylim <- rev(range(yticks))
    # annotation
    if (pretty) {
        x.args <- list(side=3, line=2)
        y.args <- list(side=4, line=3)
        t.args <- list(side=3, line=4)
        s.args <- list(side=3, line=3)
    } else {
        x.args <- list(side=1, line=par("mgp")[1])
        y.args <- list(side=2, line=par("mgp")[1])
        t.args <- list(side=3, line=1, outer=FALSE)
        s.args <- list(side=4, line=par("mgp")[1] + 1)
    }

    # # # # # # # # # # # # # # # #
    # PLOT THE AXES AND GARNISHES #
    # # # # # # # # # # # # # # # #
    do.call(plot, c(list(NA, NA), exargs))
    do.call(mtext, c(xlab, x.args, lab.args))
    do.call(mtext, c(ylab, y.args, lab.args))
    do.call(mtext, c(main, t.args, main.args))
    do.call(mtext, c(sub, s.args, sub.args))
    if (pretty && draw.axes) {
        if (!is.null(xticks[1]))  x.axis.args <- c(list(side=3, at=xticks), axis.args)
        else                      x.axis.args <- c(list(side=3), axis.args)
        if (!is.null(yticks[1]))  y.axis.args <- c(list(side=4, at=yticks), axis.args)
        else                      y.axis.args <- c(list(side=4), axis.args)
        do.call(axis, x.axis.args)
        do.call(axis, y.axis.args)
    }

    # # # # # # # # #
    # PLOT HEATMAP  #
    # # # # # # # # #
    if (heatmap) {
        if (pretty & !"colormap" %in% names(heatmap.args)) {
            heatmap.args$colormap <- plotrix::color.scale(x=0:100, cs1=c(0, 180),
                                                          cs2=100, cs3=c(25, 100),
                                                          color.spec='hcl')
        }
        if (!"add" %in% names(heatmap.args)) heatmap.args$add <- TRUE
        with(d, do.call(repulsiveForceHeatmap, c(list(f2, f1, type=v),
                                                 heatmap.args)))
    }
    if (heatmap.legend) {
        if (!"x" %in% names(heatmap.legend.args)) {
            heatmap.legend.args$x <- rep(exargs$xlim[1], 2)
            heatmap.legend.args$y <- exargs$ylim - c(0, diff(exargs$ylim) / 2)
        }
        if (!"colormap" %in% heatmap.legend.args) {
            heatmap.legend.args$colormap <- heatmap.args$colormap
        }
        do.call(repulsiveForceHeatmapLegend, heatmap.legend.args)
    }

    # # # # # # #
    # PLOT HULL #
    # # # # # # #
    if (hull.fill || hull.line) {
        hull.columns <- c('f2', 'f1', 'hull.fill.col', 'hull.line.col', 'style')
        hh <- by(d, d$gf, function(i) i[!is.na(i$f2) & !is.na(i$f1),])
        if (diphthong) {
            # if diphthong, hull should ignore diph.mean.timept
            hulls <- lapply(hh, function(i) {
                ipts <- with(i, data.frame(f2=do.call(c, f2d),
                                           f1=do.call(c, f1d)))
                hull <- ipts[chull(ipts),]
                hull$hull.fill.col <- unique(i$hull.fill.col)
                hull$hull.line.col <- unique(i$hull.line.col)
                hull$style <- unique(i$style)
                return (hull)
            })
        } else {
            hulls <- lapply(hh, function(i) with(i, i[chull(f2, f1),
                            hull.columns]))
        }
        lapply(hulls, function(i) {
            hull.args$border <- i$hull.line.col
            hull.args$col    <- i$hull.fill.col
            if (!"lty" %in% names(hull.args))  hull.args$lty <- i$style
            hull.args$x      <- cbind(i$f2, i$f1)
            with(i, do.call(polygon, hull.args))
            })
    }

    # # # # # # # # #
    # PLOT ELLIPSES #
    # # # # # # # # #
    if (ellipse.fill || ellipse.line) {
        invisible(lapply(seq_along(ellipse.points),
                         function(i) {
                             ellipse.args$border <- m$ellipse.line.col[i]
                             ellipse.args$col    <- m$ellipse.fill.col[i]
                             if (!"lty" %in% names(ellipse.args)) {
                                 ellipse.args$lty <- m$style[i]
                             }
                             do.call(polygon, c(list(x=ellipse.points[[i]]), ellipse.args))
                             }))
    }

    # # # # # # # # #
    # PLOT POLYGONS #
    # # # # # # # # #
    if (!is.na(poly.order[1]) && (poly.fill || poly.line)) {
        if (length(poly.order) != length(unique(poly.order))) {
            message("Duplicate entries in 'poly.order' detected; they will be ",
                    "ignored.")
        }
        poly.order <- as.character(poly.order) # as.character in case factor
        v <- unique(as.character(m$v))
        if (length(setdiff(poly.order, v)) > 0) {
            message("There are vowels in 'poly.order' that are not in ",
                    "'vowel'; they will be ignored.")
        }
        poly.order <- intersect(poly.order, v)
        pp <- m
        pp$v <- factor(pp$v, levels=poly.order, ordered=TRUE)
        pp <- pp[order(pp$v),]
        pp <- pp[!is.na(pp$v),]
        pp <- split(pp, pp$gf)
        if (poly.fill) {
            bigenough <- sapply(pp, function(i) nrow(i) > 2)
            lapply(pp[bigenough], function(i) {
                pargs <- poly.args
                pargs$x <- cbind(i$f2, i$f1)
                pargs$col <- i$poly.fill.col
                pargs$border <- NA
                with(i, do.call(polygon, pargs))
            })
        }
        if (poly.line) {
            if (plot.means) type <- 'c'
            else type <- 'l'
            bigenough <- sapply(pp, function(i) nrow(i) > 1)
            invisible(lapply(pp[bigenough], function(i) {
                pargs <- poly.args
                pargs$x <- i$f2
                pargs$y <- i$f1
                pargs$type <- type
                pargs$cex <- 1.2 * cex.means
                pargs$col <- i$poly.line.col
                if (!"lty" %in% names(pargs)) pargs$lty <- i$style
                with(i[i$v %in% poly.order,], do.call(points, pargs))
                }))
    }   }

    # # # # # # # #
    # PLOT TOKENS #
    # # # # # # # #
    if (plot.tokens) {
        if (diphthong) {
            # setup
            timepts <- length(d$f2d[[1]])
            if (diph.arrows) line.range <- 1:(timepts-1)
            else             line.range <- 1:timepts
            # no smoothing splines
            if (!diph.smooth || timepts < 4) {
                if (diph.smooth) warning("Cannot smooth diphthong traces with ",
                                         "fewer than 4 timepoints. Plotting ",
                                         "connecting segments instead.")
                # plot first point
                if (diph.label.first.only) {
                    if (!is.null(pch.tokens)) {
                        with(d, text(f2, f1, labels=pchtokens, col=col.tokens,
                                     cex=cex.tokens))
                    } else {
                        with(d, points(f2, f1, pch=pchtokens, col=col.tokens,
                                       cex=cex.tokens))
                    }
                    diph.args.tokens$type <- "l"
                }
                # plot lines
                apply(d, 1, function(i) {
                    # if diph.label.first.only, cex and pch will get ignored
                    with(i, do.call(points, c(list(t(f2d)[line.range],
                                                   t(f1d)[line.range],
                                                   pch=pchtokens,
                                                   cex=cex.tokens,
                                                   col=col.tokens),
                                              diph.args.tokens)))
                    })
                # plot arrowheads
                if (diph.arrows) {
                    apply(d, 1, function(i) {
                        with(i, do.call(arrows, c(list(x0=t(f2d)[timepts-1],
                                                       y0=t(f1d)[timepts-1],
                                                       x1=t(f2d)[timepts],
                                                       y1=t(f1d)[timepts],
                                                       col=col.tokens),
                                                  diph.arrow.tokens)
                        ))
                    })
                }
            # diphthong smoothing spline
            } else if (diph.smooth) {  # timepts > 3
                apply(d, 1, function(i) {
                    tryCatch({
                        steep <- with(i, abs(lm(f1d~f2d)$coefficients['f2d']) > 1)
                        if (steep) dat <- with(i, cbind(f1d, f2d))
                        else       dat <- with(i, cbind(f2d, f1d))
                        pc <- prcomp(dat, center=FALSE, scale.=FALSE)
                        ss <- smooth.spline(pc$x)
                        ssi <- as.matrix(as.data.frame(predict(ss))) %*% solve(pc$rotation)  #* pc$scale + pc$center
                        end <- nrow(ssi)
                        if (diph.arrows) {
                            curve.range <- 1:(end-1)
                            with(as.data.frame(ssi),
                                 do.call(arrows, c(list(x0=f2[end-1],
                                                        y0=f1[end-1],
                                                        x1=f2[end],
                                                        y1=f1[end],
                                                        col=i$col.tokens),
                                                   diph.arrow.tokens)))
                        } else {
                            curve.range <- 1:end
                        }
                        do.call(lines, c(list(ssi[curve.range]), diph.args.tokens))
                    },
                    error=function(e){
                        message("Warning: could not plot diphthong smoother. ",
                                "Plotting connecting segments instead.")
                        message(paste(e, ""))
                        if (diph.arrows) {
                            end <- nrow(i)
                            with(i, points(f2[1:end-1], f1[1:end-1], col=col.tokens,
                                           pch=pchtokens, cex=cex.tokens, type="o"))
                            with(i, do.call(arrows, c(list(x0=f2[end-1], y0=f1[end-1],
                                                           x1=f2[end], y1=f1[end],
                                                           col=col.tokens), diph.arrow.args)))
                        } else {
                            with(i, points(f2, f1, col=col.tokens, pch=pchtokens, cex=cex.tokens, type="o"))
                        }
                    },
                    warning=function(w) message(w),
                    finally={}
                    )
                })
            }
        } else {  # !diphthong
            if (is.null(pch.tokens)) {
                with(d, points(f2, f1, pch=pchtokens, cex=cex.tokens, col=col.tokens))
            } else {
                with(d, text(f2, f1, labels=pchtokens, cex=cex.tokens, col=col.tokens))
    }   }   }

    # # # # # # # #
    # PLOT MEANS  #
    # # # # # # # #
    if (plot.means) {
        if (diphthong) {
            # TODO: implement smoothing splines for means
            # setup
            timepts <- length(m$f2d[[1]])
            if (diph.arrows) line.range <- 1:(timepts-1)
            else             line.range <- 1:timepts

            # plot first point
            if (diph.label.first.only) {
                if (!is.null(pch.means)) {
                    with(m, text(f2, f1, labels=pchmeans, col=col.means,
                                 cex=cex.means))
                } else {
                    with(m, points(f2, f1, pch=pchmeans, col=col.means,
                                   cex=cex.means))
                }
                diph.args.means$type <- "l"
            }
            # plot lines
            apply(m, 1, function(i) {
                # if diph.label.first.only, cex and pch will get ignored
                with(i, do.call(points, c(list(t(f2d)[line.range],
                                               t(f1d)[line.range],
                                               pch=pchmeans,
                                               cex=cex.means,
                                               col=col.means),
                                          diph.args.means)))
            })
            # plot arrowheads
            if (diph.arrows) {
                apply(m, 1, function(i) {
                    with(i, do.call(arrows, c(list(x0=t(f2d)[timepts-1],
                                                   y0=t(f1d)[timepts-1],
                                                   x1=t(f2d)[timepts],
                                                   y1=t(f1d)[timepts],
                                                   col=col.means),
                                              diph.arrow.means)
                    ))
                })
            }
        } else {
            if (is.null(pch.means)) {
                with(m, points(f2, f1, col=col.means, pch=pchmeans, cex=cex.means))
            } else {
                with(m, text(f2, f1, labels=pchmeans, col=col.means, cex=cex.means))
            }
        }
    }

    # # # # # #
    # LEGEND  #
    # # # # # #
    if (!is.null(legend.kwd)) {
        if (is.null(legend.col.lab) & is.null(legend.style.lab)) {
            warning("Legend will not be drawn because var.col.by and var.style.by ",
                    "are both NULL or NA. You will have to use the legend() function.")
        } else {
            legend.merge <- TRUE
            # legend pch
            legend.pch <- NULL
            if (length(legend.style.lab)) {
                if (plot.means && all(grepl("[[:digit:]]", pch.means))) {
                    legend.pch <- unique(m$pchmeans)
                } else if (plot.tokens && all(grepl("[[:digit:]]", pch.tokens))) {
                    legend.pch <- unique(d$pchtokens)
            }   }
            # legend col
            legend.col <- NULL
            if (length(legend.col.lab)) {
                if (plot.means)       legend.col <- unique(m$col.means)
                else if (plot.tokens) legend.col <- unique(d$col.tokens)
            }
            # legend background fill
            legend.bgf <- NULL
            if (hull.fill || poly.fill || ellipse.fill) {
                if (!is.na(m$ellipse.fill.col[1]))   legend.bgf <- unique(m$ellipse.fill.col)
                else if (!is.na(m$poly.fill.col[1])) legend.bgf <- unique(m$poly.fill.col)
                else                                 legend.bgf <- unique(hull.fill.col)
            }
            # legend linteype & border color
            legend.lty <- NULL
            legend.brd <- NULL
            if (hull.line | poly.line | ellipse.line) {
                legend.lty <- unique(m$style)
                if (!is.na(m$ellipse.line.col[1]))   legend.brd <- unique(m$ellipse.line.col)
                else if (!is.na(m$poly.line.col[1])) legend.brd <- unique(m$poly.line.col)
                else                                 legend.brd <- unique(hull.line.col)
                if (length(legend.brd) != length(legend.bgf)) legend.brd <- NULL
            }
            # handle lty specially; needed for both style & color
            if (!is.null(legend.lty)) {
                if (!length(legend.style.lab)) {
                    legend.lty <- rep(legend.lty, length.out=length(legend.col.lab))
                } else if (length(legend.col.lab)) {
                    legend.lty <- c(legend.lty, rep(NA, length(legend.col.lab)))
                }
            }
            # reconcile
            if (identical(legend.style.lab, legend.col.lab)) {
                legend.lab <- legend.col.lab
                legend.merge <- FALSE
            } else {
                legend.lab <- c(legend.style.lab, legend.col.lab)
                legend.pch <- c(legend.pch, rep(NA, length(legend.col.lab)))
                # handle case where no lines / fills / pchs for color
                if (length(legend.col.lab) && is.null(legend.bgf) &&
                        is.null(legend.brd) && is.null(legend.lty)) {
                    legend.bgf <- c(rep(NA, length(legend.style.lab)), legend.col)
                    legend.brd <- legend.bgf
                    legend.col <- c(rep(par("fg"), length(legend.style.lab)), legend.col)
                # handle case where only hulls (only 1 fill col) but col.by.vowel
                } else if (length(legend.col) != length(legend.bgf) &&
                               !is.null(legend.bgf)) {
                    legend.bgf <- c(rep(NA, length(legend.style.lab)), legend.col)
                    legend.brd <- legend.bgf
                    legend.col <- c(rep(par("fg"), length(legend.style.lab)), legend.col)
                # handle other cases
                } else {
                    legend.bgf <- c(rep(NA, length(legend.style.lab)), legend.bgf)
                    legend.brd <- c(rep(NA, length(legend.style.lab)), legend.brd)
                    legend.col <- c(rep(par("fg"), length(legend.style.lab)), legend.col)
                }
            }
            # eliminate vacuous args
            if (identical(legend.bgf, logical(0))) legend.bgf <- NULL
            if (identical(legend.brd, logical(0))) legend.brd <- legend.bgf
            # assemble legend args
            new.legend.args <- list(legend.kwd, legend=legend.lab, pch=legend.pch,
                                    col=legend.col, lty=legend.lty)
            # user override & recombine
            new.legend.args[names(legend.args)] <- legend.args
            legend.args <- new.legend.args
            # can't always pass fill because fill=NULL triggers drawing empty boxes
            # and border=NULL draws black! :(
            if (!is.null(legend.bgf)) {
                if (!"fill" %in% names(legend.args))   legend.args$fill   <- legend.bgf
                if (!"border" %in% names(legend.args)) legend.args$border <- legend.brd
            }
            # avoid warning that merge only works when segments are drawn
            if (!is.null(legend.lty)) {
                if (!"merge" %in% names(legend.args)) legend.args$merge <- legend.merge
            }
            # draw legend
            do.call(legend, legend.args)
        }
    }

    # # # # # #
    # CLEANUP #
    # # # # # #
    # close file devices
    if (output != "screen") dev.off()
    # reset graphical parameters to defaults
    par(op)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# OMNIBUS NORMALIZATION FUNCTION (convenience function) #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#' @export
normVowels <- function(method, f0=NULL, f1=NULL, f2=NULL, f3=NULL,
                    vowel=NULL, group=NULL, ...) {
    m <- tolower(method)
    methods <- c('bark', 'mel', 'log', 'erb', 'zscore', 'lobanov',
                'logmean', 'shared', 'nearey1', 'nearey2', 'scentroid',
                'wattfabricius')
    if (!(m %in% methods)) {
        stop('Method must be one of: bark, mel, log, erb, ',
             'zscore | lobanov, logmean | nearey1, shared | nearey2, ',
             'scentroid | wattfabricius.')
    }
    f <- cbind(f0=f0, f1=f1, f2=f2, f3=f3)
    if (m %in% 'bark') return(normBark(f))
    else if (m %in% 'mel') return(normMel(f))
    else if (m %in% 'log') return(normLog(f))
    else if (m %in% 'erb') return(normErb(f))
    else if (m %in% c('z','zscore','lobanov')) return(normLobanov(f, group))
    else if (m %in% c('logmean','nearey1')) return(normLogmean(f, group, ...))
    else if (m %in% c('shared','nearey2')) return(normSharedLogmean(f, group, ...))
    else {
        f <- as.matrix(cbind(f1=f1, f2=f2))
        return(normWattFabricius(f, vowel, group))
}	}


# INDIVIDUAL NORMALIZATION FUNCTIONS
#' @export
normBark <- function(f) {
    bark <- 26.81 * f / (1960 + f) - 0.53
    bark[bark < 2] <- bark[bark < 2] + 0.15 * (2 - bark[bark < 2])
    bark[bark > 20.1] <- bark[bark > 20.1] + 0.22 * (bark[bark > 20.1] - 20.1)
	return(bark)
}

#' @export
normLog <- function(f) {
    log10(f)
}

#' @export
normMel <- function(f) {
    2595*log10(1+f/700)
}

#' @export
normErb <- function(f) {
    21.4*log10(1+0.00437*f)
}

#' @export
normLobanov <- function(f, group=NULL) {
    if (is.null(group)) {
        return(scale(f))
    } else {
        f <- as.data.frame(f)
        groups <- split(f, group)
        scaled <- lapply(groups, function(x) as.data.frame(scale(x)))
        return(unsplit(scaled, group))
}	}

#' @export
normLogmean <- function(f, group=NULL, exp=FALSE, ...) {
    # AKA "Nearey1", what Adank confusingly calls "SingleLogmean".
    # Note that Adank et al 2004 (eq. 8) looks more like a shared logmean,
    # but in text she says it is applied to each formant separately.
    if (ncol(f) < 2) {
        stop("Missing values: normalization method 'normLogmean' ",
            "requires non-null values for f1 and f2).")
    }
    if (is.null(group)) {
        logmeans <- log(f) - rep(colMeans(log(f), ...), each=nrow(f))
    } else {
        f <- as.data.frame(f)
        groups <- split(f, group)
        logmeans <- lapply(groups,
                           function(x) log(x) - rep(apply(log(x), 2, mean),
                                                    each=nrow(x)))
    logmeans <- unsplit(logmeans, group)
	}
    if (exp) logmeans <- exp(logmeans)
    logmeans
}

#' @export
normNearey1 <- function(f, group=NULL, exp=FALSE, ...) {
    normLogmean(f, group=group, exp=exp, ...)
}

#' @export
normSharedLogmean <- function(f, group=NULL, exp=FALSE, ...) {
    # AKA "Nearey2"
    if (is.null(group)) {
        # this is my implementation of Nearey 1978's CLIH (eq. 3.1.10, page 95)
        # (cf. eqs. 1 & 2 of Morrison & Nearey 2006)
        logmeans <- log(f) - mean(log(unlist(f)), ...)
        # NOTE: Adank et al 2004 (eq. 9) would suggest this implementation:
        # logmeans <- log(f) - sum(colMeans(log(f), ...))
    } else {
        f <- as.data.frame(f)
        groups <- split(f, group)
        logmeans <- lapply(groups, function(x) log(x) - mean(log(unlist(x)), ...))
        # Adank:    lapply(groups, function(x) log(x) - sum(colMeans(log(x), ...)))
        logmeans <- unsplit(logmeans, group)
	}
    if (exp) logmeans <- exp(logmeans)
    logmeans
}

#' @export
normNearey2 <- function(f, group=NULL, exp=FALSE, ...) {
    normSharedLogmean(f, group=group, exp=exp, ...)
}

#' @export
normWattFabricius <- function(f, vowel, group=NULL) {
    if (ncol(f) != 2) {
        stop("Wrong dimensions: s-centroid normalization requires an Nx2 ",
             "matrix or data frame of F1 and F2 values.")
    }
    if (is.null(group)) group <- rep("g", nrow(f))
    subsets <- by(f, list(vowel, group), identity)  # 2D list (group x vowel) of lists (f1,f2)
    means <- matrix(sapply(subsets, function(i) ifelse(is.null(i),
                                                    data.frame(I(c(f1=NA, f2=NA))),
                                                    data.frame(I(colMeans(i))))
                        ), nrow=nrow(subsets), dimnames=dimnames(subsets))
    minima <- do.call(rbind, apply(means, 2, function(i) data.frame(t(apply(do.call(rbind, i), 2, min, na.rm=TRUE)))))
    maxima <- do.call(rbind, apply(means, 2, function(i) data.frame(t(apply(do.call(rbind, i), 2, max, na.rm=TRUE)))))
    min.id <- apply(means, 2, function(i) apply(do.call(rbind, i), 2, which.min))
    max.id <- apply(means, 2, function(i) apply(do.call(rbind, i), 2, which.max))
    if (length(unique(min.id['f1',]))>1) {
        warning("The vowel with the lowest mean F1 value (usually /i/)",
                "does not match across all speakers/groups. You'll ",
                "have to calculate s-centroid manually.")
        print(data.frame(minF1=minima["f1",],
                vowel=dimnames(means)[[1]][min.id["f1",]],
                group=dimnames(means)[[2]]))
        stop()
    } else if (length(unique(max.id["f1",]))>1) {
        warning("The vowel with the highest mean F1 value (usually /a/) ",
                "does not match across all speakers/groups. You'll ",
                "have to calculate s-centroid manually.")
        print(data.frame(maxF1=round(maxima["f1",]),
                vowel=dimnames(means)[[1]][max.id["f1",]],
                group=dimnames(means)[[2]]))
        stop()
    }
    lowvowf2 <- do.call(rbind, means[unique(max.id["f1",]),])[,"f2"]
    centroids <- t(rbind(f1=(2*minima$f1 + maxima$f1) / 3,
    				     f2=(minima$f2 + maxima$f2 + lowvowf2) / 3))
    f / centroids[group,]
}


# # # # # # # # # # # # # # # # #
# SECONDARY ANALYSIS FUNCTIONS  #
# # # # # # # # # # # # # # # # #

#' @export
vowelMeansPolygonArea <- function(f1, f2, vowel, poly.order, group=NULL) {
    if (length(poly.order) != length(unique(poly.order))) {
        warning("Duplicate entries in 'poly.order' detected; they will be ",
                "ignored.")
    }
    poly.order <- unique(as.character(poly.order)) # as.character in case factor
    v <- unique(as.character(vowel))
    if (length(setdiff(poly.order, v)) > 0) {
        warning("There are vowels in 'poly.order' that are not in ",
                "'vowel'; they will be ignored.")
        poly.order <- intersect(poly.order, v)
    }
    if (is.null(group))  group <- "all.points"
    df <- data.frame(f2=f2, f1=f1, v=factor(vowel, levels=poly.order), g=group)
    df <- df[order(df$v),]
    bygrouparea <- c(by(df, df$g, function(i) splancs::areapl(cbind(tapply(i$f2, i$v, mean),
                                                                    tapply(i$f1, i$v, mean)))))
}

#' @export
convexHullArea <- function(f1, f2, group=NULL) {
    if (is.null(group))  group <- "all.points"
    df <- data.frame(x=f2, y=f1, g=group, stringsAsFactors=FALSE)
    bygrouppts <- by(df, df$g, function(i) i[chull(i$x, i$y), c("x", "y")])
    bygrouparea <- sapply(bygrouppts, function(i) {
        splancs::areapl(as.matrix(data.frame(x=i$x, y=i$y, stringsAsFactors=FALSE)))
})  }

#' @export
repulsiveForce <- function(x, y, type, xform=log, exclude.inf=TRUE) {
    dmat <- as.matrix(dist(cbind(x, y)))
    if (exclude.inf) dmat[dmat == 0] <- min(dmat[dmat>0]) / 2
    force <- sapply(seq_along(type), function(i) {
        sum(1 / dmat[i, !(type %in% type[i])] ^ 2)
    })
    if (!is.null(xform)) force <- xform(force)
    force
}

#' @export
repulsiveForceHeatmap <- function(x, y, type=NULL, xform=log, exclude.inf=TRUE,
                                  resolution=10, colormap=NULL, fast=FALSE, ...) {
    # default to grayscale
    if (is.null(colormap)) colormap <- plotrix::color.scale(x=0:100, cs1=0, cs2=0,
                                                            cs3=c(25,100),
                                                            color.spec="hcl")
    # create grid encompassing vowel space
    gridlist <- createGrid(x, y, resolution)
    xgrid <- gridlist$x
    ygrid <- gridlist$y
    grid <- gridlist$g
    grid$z <- NA
    grid$v <- NA
    # if using fast method, pre-calculate force of vowel tokens
    force <- repulsiveForce(x, y, type, xform, exclude.inf)
    if (fast) vertices <- data.frame(x=x, y=y, z=force)[!is.na(x) & !is.na(y),]
    else      vertices <- data.frame(x=x, y=y, z=type)[!is.na(x) & !is.na(y),]
    if (fast) {
        # segregate grid points based on Delaunay triangulation of vowel space
        triang.obj <- with(vertices, deldir::deldir(x, y, z=z, suppressMsge=TRUE))
        triangs <- deldir::triang.list(triang.obj)
        sg.idxs <- lapply(triangs, function(i) splancs::inpip(grid, i[c("x", "y")],
                                                              bound=TRUE))
        subgrids <- lapply(sg.idxs, function(i) grid[i,])
        # ignore triangles too small to contain any grid points
        vacant <- sapply(subgrids, function(i) dim(i)[1] == 0)
        triangs <- triangs[!vacant]
        subgrids <- subgrids[!vacant]
        sg.idxs <- sg.idxs[!vacant]
        # assign a force value to each grid point
        sg.force <- mapply(function(i, j) fillTriangle(i[,"x"], i[,"y"],
                                                       j[, c("x", "y", "z")]),
                           subgrids, triangs, SIMPLIFY=FALSE)
        grid[do.call(c, sg.idxs), "z"] <- do.call(c, sg.force)
    } else {
        if (is.null(type)) stop("More accurate force calculation method (fast=FALSE) ",
                                "requires non-null values for 'type'.")
        hull <- vertices[chull(vertices),]  # polygon of hull
        sg.idxs <- splancs::inpip(grid[,1:2], hull)
        subgrid <- grid[sg.idxs,]
        # which vowel is closest to each grid point?
        dmat <- apply(subgrid, 1, function(i) apply(as.matrix(vertices[c('x','y')]),
                                                    1, function(j) dist(rbind(i, j))))
        # still need to exclude inf, in case duplicate pts have diff vowels
        if (exclude.inf) dmat[dmat == 0] <- min(dmat[dmat>0]) / 2
        # if there is a tie of which vowel is closest, pick first one (arbitrarily)
        which.min <- apply(dmat, 2, function(i) which(i == min(i))[1])
        #how.many.min <- apply(dmat, 2, function(i) length(which(i == min(i))))
        # TODO: sensible way to do tiebreaker? wherever how.many.min > 1, look at
        # next closest vowel iteratively until you find one that matches one of the
        # vowels tied as closest... could get ugly.
        sg.vowel <- vertices$z[which.min]
        sg.force <- sapply(seq_along(sg.vowel),
                           function(i) sum(1/dmat[!(vertices$z %in% sg.vowel[i]), i] ^ 2))
        if (!is.null(xform)) sg.force <- xform(sg.force)
        grid[sg.idxs, "z"] <- sg.force
    }
    image(xgrid, ygrid, matrix(grid$z, nrow=length(xgrid)),
          col=colormap, ...)
}

#' @export
repulsiveForceHeatmapLegend <- function (x, y, labels=c("low", "high"), pos=c(1, 3),
                                         colormap=NULL, smoothness=50, lend=2,
                                         lwd=12, ...) {
    if (is.null(colormap)) {  # default to grayscale
        colormap <- plotrix::color.scale(x=0:100, cs1=0, cs2=0, cs3=c(0,100),
                                         alpha=1, color.spec="hcl")
    }
    xvals <- seq(x[1], x[2], length.out=smoothness)
    yvals <- seq(y[1], y[2], length.out=smoothness)
    cols <- colormap[round(seq(1, length(colormap), length.out=smoothness))]
    invisible(plotrix::color.scale.lines(xvals, yvals, col=cols, lend=lend,
                                         lwd=lwd, ...))
    if (!is.null(labels)) {
        text(x, y, labels=labels, pos=pos, xpd=TRUE)
    }
}


# # # # # # # # # # # # # # # # # #
# UTILITY FUNCTIONS: NOT EXPORTED #
# # # # # # # # # # # # # # # # # #
prettyTicks <- function(lim) {
    axrange <- abs(diff(lim))
    step <- 10^(floor(log(axrange,10)))
    coef <- ifelse(axrange/step < 1, 0.1,
                   ifelse(axrange/step < 2, 0.2,
                          ifelse(axrange/step < 5, 0.5, 1)))
    step <- step*coef
    lims <- c(ceiling(max(lim)/step)*step, floor(min(lim)/step)*step)
    if (diff(lims) < 0) {step <- -step}
    seq(lims[1],lims[2],step)
}

ellipse <- function(mu, sigma, alpha=0.05, npoints=250, draw=TRUE, ...) {
    # adapted from the (now-defunct) mixtools package
    if (all(sigma == matrix(rep(0, 4), nrow=2))) return(rbind(mu, mu))
    es <- eigen(sigma)
    e1 <- es$vec %*% diag(sqrt(es$val))
    r1 <- sqrt(qchisq(1 - alpha, 2))
    theta <- seq(0, 2 * pi, len=npoints)
    v1 <- cbind(r1 * cos(theta), r1 * sin(theta))
    pts <- t(mu - (e1 %*% t(v1)))
    if (draw) {
        colnames(pts) <- c("x", "y")
        polygon(pts, ...)
    }
    invisible(pts)
}

makeTransparent <- function (color, opacity) {
    rgba <- t(col2rgb(color, alpha=TRUE))
    rgba[,4] <- round(255 * opacity)
    colnames(rgba) <- c("red", "green", "blue", "alpha")
    trans.color <- do.call(rgb, c(as.data.frame(rgba), maxColorValue=255))
}

# pineda's triangle filling algorithm
fillTriangle <- function(x, y, vertices) {
    x0 <- vertices[1,1]
    x1 <- vertices[2,1]
    x2 <- vertices[3,1]
    y0 <- vertices[1,2]
    y1 <- vertices[2,2]
    y2 <- vertices[3,2]
    z0 <- vertices[1,3]
    z1 <- vertices[2,3]
    z2 <- vertices[3,3]
    e0xy <- (x-x0)*(y1-y0)-(y-y0)*(x1-x0)
    e1xy <- (x-x1)*(y2-y1)-(y-y1)*(x2-x1)
    e2xy <- (x-x2)*(y0-y2)-(y-y2)*(x0-x2)
    e0x2 <- (x2-x0)*(y1-y0)-(y2-y0)*(x1-x0)
    e1x0 <- (x0-x1)*(y2-y1)-(y0-y1)*(x2-x1)
    e2x1 <- (x1-x2)*(y0-y2)-(y1-y2)*(x0-x2)
    f0 <- e0xy / e0x2
    f1 <- e1xy / e1x0
    f2 <- e2xy / e2x1
    z <- f0*z2 + f1*z0 + f2*z1
}

# create grid encompassing all pts (x, y) with `resolution` pts along short dimension
createGrid <- function(x, y, resolution) {
    vertices <- data.frame(x=x, y=y)
    vertices <- vertices[!is.na(vertices$x) & !is.na(vertices$y),]
    bounding.rect <- apply(vertices[c("x", "y")], 2, range, na.rm=TRUE)
    xr <- abs(diff(bounding.rect[,"x"]))
    yr <- abs(diff(bounding.rect[,"y"]))
    if (xr > yr) {
        xres <- round(resolution * xr / yr)
        yres <- resolution
    } else {
        xres <- resolution
        yres <- round(resolution * yr / xr)
    }
    xgrid <- seq(floor(bounding.rect[1, 1]), ceiling(bounding.rect[2, 1]),
                 length.out=xres)
    ygrid <- seq(floor(bounding.rect[1, 2]), ceiling(bounding.rect[2, 2]),
                 length.out=yres)
    grid <- expand.grid(x=xgrid, y=ygrid)
    list(g=grid, x=xgrid, y=ygrid)
}
