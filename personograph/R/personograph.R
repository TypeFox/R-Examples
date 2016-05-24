#' Generate personograph plots from data
#'
#' A personograph (Kuiper-Marshall plot) is a pictographic
#' representation of (relative) harm and benefit from an intervention. It is
#' similar to
#' \href{http://www.nntonline.net/visualrx/examples/}{Visual Rx (Cates
#' Plots)}. Each icon on the grid is colored to indicate whether that
#' percentage of people is harmed by the intervention, would benefit from the
#' intervention, has good outcome regardless of intervention, or bad outcome regardless of
#' intervention.
#' This terminology is similar to that of Uplift Modelling.
#'
#' The plot function \code{\link{personograph}} is implemented in such
#' a way that it's easy to just pass a named list of percentages,
#' colors, and an icon. Making it potentially useful for other use
#' cases as well.
#'
#' \if{html}{
#' The example code will generate the following graph if \code{higher_is_better=F}:
#'
#' \figure{green.png}{}
#'
#' }
#' \if{latex}{
#' The example code will generate the following graph if \code{higher_is_better=F}:
#'
#' \figure{green.pdf}{options: width=5in}
#' }
#'
#' \subsection{Funding & Acknowledgments}{
#' This software was commissioned and sponsored by \href{http://www.doctorevidence.com/}{Doctor Evidence}.
#' The Doctor Evidence mission is to improve clinical outcomes by
#' finding and delivering medical evidence to healthcare
#' professionals, medical associations, policy makers and
#' manufacturers through revolutionary solutions that enable anyone to
#' make informed decisions and policies using medical data that is
#' more accessible, relevant and readable.}
#'
#' \subsection{Source & Issues}{
#' Source code and issue tracker can be found on \href{https://github.com/joelkuiper/personograph}{Github}.
#' }
#'
#' @docType package
#' @name personograph-package
#' @seealso \code{\link{personograph}}
#' @seealso \code{\link{uplift}}
#' @import grid
#' @import grImport
#' @import grDevices
#' @importFrom utils tail
#' @examples
#' # Example data
#' data <- read.table(textConnection('
#'           name ev.trt n.trt ev.ctrl n.ctrl
#' 1     Auckland     36   532      60    538
#' 2        Block      1    69       5     61
#' 3        Doran      4    81      11     63
#' 4        Gamsu     14   131      20    137
#' 5     Morrison      3    67       7     59
#' 6 Papageorgiou      1    71       7     75
#' 7      Tauesch      8    56      10     71
#' '
#' ), header=TRUE)
#'
#' sm <- "RR" # The outcome measure (either Relative Risk or Odds Ratio)
#' if (requireNamespace("meta", quietly = TRUE)) { # use meta if available
#'     ## Calculate the pooled OR or RR point estimate
#'     m <- with(data, meta::metabin(ev.trt, n.trt, ev.ctrl, n.ctrl, sm=sm))
#'     point <- exp(m$TE.random) # meta returns random effects estimate on the log scale
#' } else {
#'     # Calculated Random Effects RR, using the meta package
#'     point <- 0.5710092
#' }
#'
#' # Approximate the Control Event Rates using a weighted median
#' cer <- w.approx.cer(data[["ev.ctrl"]], data[["n.ctrl"]])
#'
#' # Calculate the Intervention Event Rates (IER) from the CER and point estimate
#' ier <- calc.ier(cer, point, sm)
#'
#' # Calcaulte the "uplift" statistics
#' # Note that this depends on the direction of the outcome effect (higher_is_better)
#' u <- uplift(ier, cer, higher_is_better=FALSE)
#' plot(u, fig.title="Example", fig.cap="Example")
NULL

w.median <- function(x, w) {
    ## Lifted from cwhmisc, http://www.inside-r.org/packages/cran/cwhmisc/docs/w.median
    if (missing(w)) w <- rep(1,length(x))
    ok <- stats::complete.cases(x, w)
    x <- x[ok]
    w <- w[ok]
    ind <- sort.list(x)
    x <- x[ind]
    w <- w[ind]
    ind1 <- min(which(cumsum(w) / sum(w) >= 0.5))
    ind2 <- if((w[1] / sum(w)) > 0.5) {
        1
    } else {
        max(which(cumsum(w) / sum(w) <= 0.5))
    }
    max(x[ind1], x[ind2])
}

#' Calculate the CER (Control Event Rates)
#'
#' Calculates the CER from the data, this is a approximation of absolute
#' risk in the control population (from 0 to 1).
#'
#' By default it uses a weighted median of the indivdual control event rates. The weighted median has the benefit of always returning
#' an event rate that actually did occur. However, it is possible that this might return a CER of 0.
#' In this case we fall back to a weighted mean, and throw a warning.
#' If this too returns a CER of 0, it probably means that there was not enough data to estimate the control risk accurately.
#' In this case we recommend you obtain an estimate of the risk in the control group, for example from an observational study or expert opinion.
#'
#' @export
#' @param ev.ctrl Vector of event rates in the control group (/arm)
#' @param n.ctrl Vector of sample sizes in the control group (/arm)
#' @return Approximated Control Event Rates (CER)
w.approx.cer <- function(ev.ctrl, n.ctrl) {

    study_cer <- ev.ctrl / n.ctrl
    result <- w.median(study_cer, n.ctrl)
    if(result == 0) { # Weighted median was zero
        result <- sum(ev.ctrl) / sum(n.ctrl)
        warning("The estimated control arm risk from the studies using the weighted median was 0, using the weighted mean instead.", call.=F)
    }

    if(result == 0) { # Still zero!
        warning("The control arm risk was estimated as 0 from the studies. This probably indicates that the true control risk is a very small value greater than zero, but the studies were not able to estimate it accurately. Please note this may lead to unexpected results in the personograph chart. We recommend you obtain an estimate of the risk in the control group, for example from an observational study or expert opinion.", call.=F)
    }
    result
}

#' Calculate the IER (Intervention Event Rates)
#'
#' @export
#' @seealso \code{\link{w.approx.cer}}
#' @param cer Absolute risk with control (calculated; from 0 to 1)
#' @param point Relative risk with intervention (direct from meta-analysis)
#' @param sm The outcome measure, RR or OR as string
#' @return Absolute risk of intervention as Intervention Event Rates (IER)
calc.ier <- function(cer, point, sm) {
    if (sm == "RR") {
        return(cer * point)
    } else if(sm == "OR") {
        return(cer * (point / (1 - (cer * (1 - point)))))
    } else {
        stop("sm need to be OR (Odds Ratios) or RR (Relative Risk)")
    }
}

#' "Uplift" from IER and CER
#'
#' Calculates the percentage (from 0 to 1) of people who have an intervention benefit, intervention harm, bad outcome regardless, and good outcome regardless
#' from the Intervention Event Rates (IER) and Control Event Rates (CER).
#' Note that the result depends on the direction of the outcome measure,
#' e.g. \code{higher_is_better = T} (default) for intervention efficacy, \code{higher_is_better = F} for
#' adverse events.
#'
#' The adopted terminology is similar to that of Uplift modelling
#' \url{https://en.wikipedia.org/wiki/Uplift_modelling}
#'
#' @export
#' @param ier Intervention Event Rates
#' @param cer Control Event Rates
#' @param higher_is_better logical indicating the direction of the outcome measure, default TRUE
#' @return A list of S3 class \code{personograph.uplift} with the following elements:
#' \itemize{
#' \item{\code{good outcome}} {people who have a good outcome regardless of intervention}
#' \item{\code{bad outcome}} {people who have a bad outcome regradless of intervention}
#' \item{\code{intervention benefit}} {people who benefit from intervention}
#' \item{\code{intervention harm}} {people who are harmed by intervention}
#' }
#'
#' Can be plotted as a personograph with the S3 generic \code{plot}.
#' @examples
#' ier <- 0.06368133
#' cer <- 0.1115242
#' u <- uplift(ier, cer, higher_is_better=TRUE)
#' plot(u)
uplift <- function(ier, cer, higher_is_better=NULL) {
    if(is.null(higher_is_better)) {
        higher_is_better <- T
        warning("Setting higher_is_better as outcome direction to TRUE")
    }
    if (higher_is_better == F) {
        ## Always orient the numbers so that higher events represents a good outcome
        ier <- 1 - ier
        cer <- 1 - cer
    }

    ## [good outcome] people who have a good outcome no matter what intervention
    good <- min(ier, cer)

    ## [bad outcome] people who have bad outcome no matter what intervention
    bad <- 1 - max(ier, cer)

    ## [intervention benefit] people who would be benefit from the intervention
    benefit <- max(ier - cer, 0)

    ## [intervention harm] people who would harmed by intervention
    harm <- max(cer - ier, 0)

    result <- list("good outcome"=good,
                  "bad outcome"=bad,
                  "intervention harm"=harm,
                  "intervention benefit"=benefit)

    class(result) <- "personograph.uplift"
    result
}

as.colors <- function(lst, palette=gray.colors) {
    n <- names(lst)
    colors <- palette(length(n))
    sapply(n, function(name) { colors[[which(n == name)]]}, simplify = FALSE, USE.NAMES = TRUE)
}

round.standard <- function(x) {
    ## rounds numbers conventionally
    ## so that round.standard(0.5)==1
    floor(x + sign(x) * 0.5)
}

round.with.warn <- function(x, f=round.standard, name=NULL) {
    rounded <- f(x)
    if(x > 0 && rounded == 0) {
        warning(paste("truncating", ifelse(is.null(name), "a", name), "non-zero value of", x, "to 0"))
    }
    rounded
}

naturalfreq <- function(ar, denominator=100) {
    numerator <- ar
    if(numerator > 0 && numerator <0.5) {
        return(paste0("< 1/", denominator))
    } else {
        return(paste0(round.standard(numerator), "/", denominator))
    }
}

setColor <- function(icon, color) {
    for(i in seq_along(icon@paths)) { icon@paths[[i]]@rgb <- color}
    icon
}

#' Plots a personograph
#'
#' Plots a personograph from a named list with percentages (must sum to
#' 1). A personograph is a graphical represenation of relative benefit
#' or harm, using a grid of icons with different colors. Its intended
#' use is similar to that of Cates Plots (Visual Rx, Number Needed to
#' Treat visualization).
#' Although these could be seen as Kuiper-Marshall plots.
#'
#' \subsection{Supplying your own icon}{
#' You can supply your own icon by setting \code{icon} to a \code{grImport} \code{Picture}.
#' A \code{Picture} can be loaded with \code{grImport::readPicture} which requires a \code{grImport} XML file.
#' Obtaining this file from a standard SVG or PDF graphics file requires conversion.
#' The easiest way is to convert your original file to PDF and then to PostScript (PS) with the command-line \code{pdf2ps} tool, then tracing it with \code{grImport::PostScriptTrace}.
#' See the \code{grImport} package documentation for more details.}
#'
#' @export personograph
#' @param data A list of names to percentages (from 0 to 1)
#' @param icon.style A numeric from 1-11 indicating which of the included icons to use, they are mostly variations on the theme
#' @param icon A \code{grImport} \code{Picture} for the icon, overwrites \code{icon.style}
#' @param icon.dim The dimensions of icon as a vector \code{c(width, height)} as numerical. Calculated from the \code{dimensions} if not supplied
#' @param n.icons Number of icons to draw, defaults to 100
#' @param plot.width The percentage of width that the main plotting area should take (with respect to the frame)
#' @param dimensions A vector of \code{c(rows, columns)} for the dimensions of the grid
#' @param colors A vector of names to colors, must match the names in data. Uses \code{gray.colors} style if none supplied
#' @param fig.cap Figure caption
#' @param fig.title Figure title
#' @param draw.legend Logical if TRUE (default) will draw the legend
#' @param force.fill A character vector of 'ignore' (default), 'most', 'least', or one of the names from \code{data}.
#'     Defines the behaviour for cases when the rounding doesn't add
#'     up to \code{n.icons}. 'ignore' simply draws less icons, 'most' adds an
#'     icon to the largest group, 'least' to the smallest.
#'     If a name from \code{data} is supplied it will added to that element
#' @param fudge Fudge factor for the icon size, substracted from the \code{icon.size}
#' @param round.fn Function that is applied to round the percentages from \code{data} to \code{n.icons}. See also \code{force.fill}
#' @param legend.show.zeros Logical if TRUE indicating whether to show zero (0) values in the legend.
#' @return None.
#' @examples
#' data <- list(first=0.9, second=0.1)
#' personograph(data)
#' # With colors
#' personograph(data, colors=list(first="red", second="blue"))
#' # With different icon.style
#' personograph(data, icon.style=4) # numeric from 1-11
#' # Plot a thousand in a 20x50 grid
#' personograph(data, n.icons=1000, dimensions=c(20,50), plot.width=0.75)
personograph <- function(data,
                 fig.title=NULL,
                 fig.cap=NULL,
                 draw.legend=T,
                 icon=NULL,
                 icon.dim=NULL,
                 icon.style=1,
                 n.icons=100,
                 plot.width=0.75,
                 dimensions=ceiling(sqrt(c(n.icons, n.icons))),
                 fudge=0.0075,
                 legend.show.zeros=TRUE,
                 force.fill="ignore",
                 round.fn=round.standard,
                 colors=as.colors(data)) {

    stopifnot(sum(unlist(data)) == 1)

    devAskNewPage(FALSE)
    grid.newpage()

    fontfamily <- c("Helvetica", "Arial")

    if(is.null(icon)) {
        icon <- readPicture(system.file(paste0(icon.style, ".ps.xml"), package="personograph"))
    }

    master.rows <- sum(draw.legend, !is.null(fig.cap))
    master.heights <- c(0.1,
                       0.9 - (master.rows * 0.1),
                       ifelse(draw.legend, .1, 0),
                       ifelse(!is.null(fig.cap) || !draw.legend, .1, 0))

    masterLayout <- grid.layout(
        nrow    = 4,
        ncol    = 1,
        heights = unit(master.heights, rep("null", 4)))

    vp1 <- viewport(layout.pos.row=1, layout.pos.col = 1, name="title")
    vp2 <- viewport(layout.pos.row=2, layout.pos.col = 1, name="plot")
    vp3 <- viewport(layout.pos.row=3, layout.pos.col = 1, name="legend")
    vp4 <- viewport(layout.pos.row=4, layout.pos.col = 1, name="caption", just=c("centre", "top"))

    pushViewport(vpTree(viewport(layout = masterLayout, name = "master"), vpList(vp1, vp2, vp3, vp4)))

    if(!is.null(fig.title)) {
        seekViewport("title")
        grid.text(fig.title,
                  gp = gpar(fontsize = 18, fontfamily=fontfamily, fontface="bold"))
        popViewport()
    }

    rows <- dimensions[1]
    cols <- dimensions[2]

    if(is.null(icon.dim)) {
        icon.height <- 1 / rows
        icon.width <- 1 / cols
    } else {
        icon.height <- icon.dim[1]
        icon.width <- icon.dim[2]
    }

    data.names <- names(data)

    if(is.null(colors)) {
        colors <- as.colors(data)
    }

    ## round based on *cumulative* sum
    ## this means that icons will be aligned to the grid, and
    ## will avoid problems with having more icons than the
    ## grid allows due to rounding errors
    n.elements <- length(data)
    cum_data <- cumsum(c(data, 0))

    rounded <- round(cum_data * n.icons)
    rounded[2:n.elements] <- rounded[2:n.elements] - rounded[1:n.elements-1]
    counts <- round.with.warn(rounded, f=round.fn, name=name)

    if(sum(unlist(counts)) < n.icons) {
        ordered.names <- data.names[order(unlist(counts))]
        forceFill <- function(counts, name) {
            counts[[name]] <- counts[[name]] + 1
            warning(paste0("adding an extra icon to ", name, ", to fill to ", n.icons))
            counts
        }
        if(force.fill == "least") {
            counts <- forceFill(counts, utils::tail(ordered.names, n = 1))
        } else if(force.fill == "most") {
            counts <- forceFill(counts, ordered.names[[1]])
        } else if(force.fill == "ignore") {
            warning(paste0("rounded sum of icons does not add up to ", n.icons, ", drawing less icons"))
        } else {
            counts <- forceFill(counts, force.fill)
        }
    }

    flat <- unlist(lapply(data.names, function(name) { rep(name, counts[[name]])}))

    seekViewport("plot")
    pushViewport(viewport(width=unit(plot.width, "npc")))

    colorMatrix <- function(flat, colors) {
        m <- matrix(nrow=rows, ncol=cols)
        total <- 0
        for (i in rows:1) {
            for (j in 1:cols) {
                total <- total + 1
                if(total < length(flat) + 1) {
                    j_snake <- ifelse((i %% 2 == 1), j, cols - j + 1) # to group like icons together
                    m[i,j_snake] <- colors[[flat[[total]]]]
                }
            }
        }
        m
    }

    colorMask <- function(colorMatrix, color) {
        return(ifelse(colorMatrix == color, color, NA))
    }

    coordinates <- function(colorMask, width, height) {
        originalDim <- dim(colorMask)
        rows <- originalDim[1]
        cols <- originalDim[2]

        x <- matrix(seq((width/2), 1 - (width/2), by=width), nrow=rows, ncol=cols, byrow=T)
        y <- matrix(seq((height/2), 1 - (height/2), by=height), nrow=rows, ncol=cols, byrow=F)

        list(x=x[which(!is.na(colorMask), TRUE)], y=y[which(!is.na(colorMask), TRUE)])
    }

    colorM <- colorMatrix(flat, colors)

    for(name in data.names) {
        color <- colors[[name]]
        mask <- colorMask(colorM, color)
        coords <- coordinates(mask, icon.width, icon.height)
        if(length(coords$x) > 0 && length(coords$y) > 0) {
            icon <- setColor(icon, color)
            grid.symbols(icon, x=coords$x, y=coords$y, size=max(icon.height, icon.width) - fudge)
        }
    }
    popViewport(2)

    font <- gpar(fontsize=11, fontfamily)

    if(draw.legend) {
        seekViewport("legend")
        filtered.data.names <- if(legend.show.zeros) {data.names } else {data.names[which(data != 0)]}

        legendCols <- length(filtered.data.names)
        legendGrobs <- list()
        legendWidths <- list()
        for(name in filtered.data.names) {
            label <- paste(naturalfreq(counts[[name]], denominator=n.icons), name)
            grob <- textGrob(label, gp=font, just="left", x=-0)
            legendGrobs[[name]] <- grob
            legendWidths[[name]] <- widthDetails(grob)
        }

        legendWidths <- c(rbind(rep(unit(0.25, "inches"), legendCols), unlist(legendWidths)))

        pushViewport(viewport(
            clip   = F,
            width  = unit(0.8, "npc"),
            layout = grid.layout(ncol=legendCols * 2,
                                 nrow=1,
                                 widths=unit(legendWidths, "inches"),
                                 heights=unit(0.25, "npc"))))


        idx <- 0
        for(name in filtered.data.names)  {
            idx <- idx + 1
            pushViewport(viewport(layout.pos.row=1, layout.pos.col=idx))
            grid.circle(x=0.5, r=0.35, gp=gpar(fill=colors[[name]], col=NA))
            popViewport()

            idx <- idx + 1
            pushViewport(viewport(layout.pos.row=1, layout.pos.col=idx))
            grid.draw(legendGrobs[[name]])
            popViewport()
        }

        popViewport(2)
    }

    if(!is.null(fig.cap)) {
        seekViewport("caption")
        grid.text(fig.cap, gp = font)
        popViewport()
    }

    dev.flush()
    return(invisible(NULL))
}

#' @export
#' @method plot personograph.uplift
#' @seealso \code{\link{personograph}}
plot.personograph.uplift <- function(x, ...) {
    colors <- list("intervention harm"="firebrick3", "intervention benefit"="olivedrab3", "bad outcome"="azure4", "good outcome"="azure2")
    personograph(x, colors=colors, legend.show.zeros=F, ...)
}
