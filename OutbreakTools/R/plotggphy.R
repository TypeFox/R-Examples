## hack to remove the NOTE in R CMD check about:
## plotggphy: no visible binding for global variable ‘x.beg’
## plotggphy: no visible binding for global variable ‘x.end’
## plotggphy: no visible binding for global variable ‘y.beg’
## plotggphy: no visible binding for global variable ‘y.end’
## plotggphy: no visible binding for global variable ‘y’
## plotggphy: no visible binding for global variable ‘label’
if(getRversion() >= "2.15.1")  utils::globalVariables(c("x.beg","x.end","y.beg","y.end","y","label"))


#' Function to plot phylogenies using \code{\link{ggplot2}}
#'
#' @param x An \linkS4class{obkData} object.
#' @param which.tree Numeric. Specify the order of the tree to be plotted. Currently the function cannot plot multiple tree.
#' @param ladderize If \code{TRUE}, the phylogeny is ladderized
#' @param show.tip.label Logical. If \code{TRUE}, the labels of the tip are shown.
#' @param tip.label.size Numeric. Size of the tip labels.
#' @param build.tip.attribute Logical. If \code{TRUE}, then a data frame \code{tip.attribute} is constructed by merging the data frames \code{individuals} and \code{records}.
#' @param tip.color Character. Can be either the name of a color (e.g. "Black") or the name of a column of \code{tip.attribute}. In the first case, all tips have the specified color. In the second case, tips are color-codded according to the specified attribute.
#' @param tip.alpha Character (or Numeric). Can be either the value of transparency (between 0 and 1) or the name of a column of \code{tip.attribute}. In the first case, all tips have the specified transparency. In the second case, tips are transparency-codded according to the specified attribute.
#' @param tip.shape Character (or Numeric). Can be either the value of a shape (e.g. 16 correpond to filled circles) or the name of a column of \code{tip.attribute}. In the first case, all tips have the specified shape. In the second case, tips are shape-codded according to the specified attribute.
#' @param tip.size Character (or Numeric). Can be either the value of tip size or the name of a column of \code{tip.attribute}. In the first case, all tips have the specified size. In the second case, tips are size-codded according to the specified attribute.
#' @param branch.unit Character. The unit of the branch can be either "year", "month", "day" or "subst". If a time unit is provided, together with use.tip.dates, then the x-axis of the phylogeny is plotted in date format using standard POSIX specification.
#' @param tip.dates Character. If branch.unit is in unit of time, tip.dates indicates the name of the column of tip.attribute that contains the sampling dates of the tip. See also guess.tip.dates.from.labels.
#' @param guess.tip.dates.from.labels Logical. If \code{TRUE} then tip.dates are guessed from the tip labels using the information provided by 'set.guess'.
#' @param set.guess List. A list of three elements: prefix, order and from. For instance, if labels are formated like this: A/Shenzhen/40/2009_China_2009-06-09 then
#' set.guess = list(prefix="_",order=3,from="first") or set.guess = list(prefix="_",order=1,from="last").
#' @param axis.date.format Character. When x-axis is in date format, this argument allow to change the format of the tick labels. See strptime for more details.
#' @param major.breaks Character. Major x-axis breaks (only when x is in date format). Ex: "weeks", "15days", "months", etc.
#' @param minor.breaks Character. Minor x-axis breaks (only when x is in date format). Ex: "weeks", "15days", "months", etc.
#' @param color.palette Character. The palette for tip colors. Only palettes from the package \code{\link{RColorBrewer}} are available. See \code{\link{brewer.pal}} documentation for more details.
#' @param legend.position Character (or numeric). The position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
#' @export
#' @import ggplot2,RColorBrewer
#' @author Anton Camacho
#' @example ../misc/plotggphyExample.R
#'
plotggphy <- function(x, which.tree = 1, ladderize = TRUE, show.tip.label = NULL, tip.label.size = 3, build.tip.attribute = TRUE,
	tip.color = NULL, tip.alpha = NULL, tip.shape = NULL, tip.size = NULL, branch.unit = NULL, tip.dates = NULL,
	guess.tip.dates.from.labels = FALSE, set.guess = list(prefix = "_", order = 1, from = "last"), axis.date.format = NULL,
                      major.breaks = NULL, minor.breaks = NULL, color.palette = "Spectral", legend.position = "right") {

    ## STOP IF NOT OBKDATA OBJECT ##
    if (!inherits(x, "obkData")) {
        stop("argument x must be an object of class obkData")
    }


    ## GET PLOTTED TREE ##
    phylo <- get.trees(x)
    if (!length(phylo)) {
        stop("x doesn't contain any tree")
    }
    if (length(phylo) > 1) {
        warnings("plotggphy can't plot more than one tree at the moment. By default the first tree is plotted. This can be changed using the \"which.tree\" argument.")
    }

    phylo <- phylo[[which.tree]]

    ## handle default value for show.tip.label ##
    if(is.null(show.tip.label)){
        if(length(phylo$tip.label)<=50) {
            show.tip.label <- TRUE
        } else {
            show.tip.label <- FALSE
        }
    }


    ## BUILD DATA.FRAME OF TIP ATTRIBUTES IF NEEDED ##
    if(build.tip.attribute && is.null(x@dna@meta)) warning("Cannot build tip attribute without information on the sequences (@dna@meta is empty).")
    if(build.tip.attribute && !is.null(x@dna@meta)){
        tip.attribute <- make.tip.attributes(x, which.tree=which.tree)
    } else { # not building tip attributes

        tip.attribute <- NULL

        if (!is.null(tip.color)) {
            if (!is.color(tip.color)) {
                stop("tip.color must be either the name of a column of tip.attribute or a color name (e.g. \"red\").")
            }
        }
        if (!is.null(tip.shape)) {
            if (!is.numeric(tip.shape)) {
                stop("tip.shape must be either the name of a column of tip.attribute or a numeric value.")
            }
        }
        if (!is.null(tip.size)) {
            if (!is.numeric(tip.size)) {
                stop("tip.size must be either the name of a column of tip.attribute or a numeric value.")
            }
        }
        if (!is.null(tip.alpha)) {
            if (!is.numeric(tip.alpha)) {
                stop("tip.alpha must be either the name of a column of tip.attribute or a numeric value.")
            }
        }

    }

    if (!is.null(branch.unit)) {
        if (!branch.unit %in% c("subst", "year", "month", "day")) {
            stop("branch.unit should be either NULL or one of \"subst\",\"year\",\"month\",\"day\"")
        }
    }

    if (guess.tip.dates.from.labels) {
        if (any(names(set.guess) != c("prefix", "order", "from"))) {
            stop("set.guess must be a named list containing prefix, order and from, see documentation")
        }
        ##read dates
        tip.dates <- try(.process.Date(.extract.string(phylo$tip.label, set.guess[["prefix"]], set.guess[["order"]], set.guess[["from"]])))
        if (inherits(tip.dates, "try-error")) {
            stop("tip.dates are not in an unambiguous format when extracted from tip.label using the prefix and order provided by set.guess. See documentation of as.Date().")
        }

    } else {
        if(is.null(tip.dates)) tip.dates <- tip.attribute$date
    }

    ##ladderize?
    if (ladderize) {
        phylo <- ladderize(phylo)
    }

    ##transform the phylo object to plot it with ggplot2
    ggphy <- phylo2ggphy(phylo, tip.dates = tip.dates, branch.unit = branch.unit)


    ##TODO: allow edge and node attributes and merge with df.edge and df.node
    df.tip <- ggphy[[1]]
    df.node <- ggphy[[2]]
    df.edge <- ggphy[[3]]

    is.x.date <- inherits(df.edge$x.beg, "Date")

    if (!is.null(tip.attribute)) {
        ##merge df.tip with tip.attribute
        tip.attribute$label <- rownames(tip.attribute)
        tmp <- merge(df.tip, tip.attribute, by = "label", all.x = T)
        df.tip <- tmp
    }


    ##theme_set(theme_grey())
    theme.old <- theme_update(axis.ticks.y = element_blank(), axis.title.y = element_blank(), panel.grid.major.y = element_blank(),
                              panel.grid.minor.y = element_blank())

    p <- ggplot(data = df.edge)
    p <- p + geom_segment(data = df.edge, aes(x = x.beg, xend = x.end, y = y.beg, yend = y.end), lineend = "round")
    p <- p + scale_y_continuous("", breaks = NULL)

    if (is.x.date) {
        to_parse <- paste("scale_x_date(\"Time\"", ifelse(is.null(axis.date.format), "", ",labels=date_format(axis.date.format)"),
                          ifelse(is.null(major.breaks), "", ",breaks= date_breaks(major.breaks)"), ifelse(is.null(minor.breaks),
                                                              "", ",minor_breaks= minor.breaks"), ")")

        p <- p + eval(parse(text = to_parse))

    } else if (is.null(branch.unit)) {
        p <- p + scale_x_continuous("")
    } else {
        p <- p + scale_x_continuous("Number of substitutions per generation")
    }

    if (show.tip.label) {
        p <- p + geom_text(data = df.tip, aes(x = x, y = y, label = label), hjust = 0, size = tip.label.size)
    }


    tip.characteristic <- list(color = tip.color, alpha = tip.alpha, shape = tip.shape, size = tip.size)
    is.aes <- sapply(seq_along(tip.characteristic), function(i) {
        x.val <- tip.characteristic[[i]]
        x.name <- names(tip.characteristic)[i]

        if (is.null(x.val)) {
            return(NA)
        }
        if (x.val %in% names(tip.attribute)) {
            return(TRUE)
        }
        if (x.name == "color" & is.color(x.val)) {
            return(FALSE)
        }
        if (x.name %in% c("alpha", "shape", "size") & is.numeric(x.val)) {
            return(FALSE)
        }
    })
    names(is.aes) <- names(tip.characteristic)
    is.aes <- is.aes[!is.na(is.aes)]

    tip.aes <- names(is.aes)[is.aes]
    if (length(tip.aes)) {
        tip.aes.txt <- paste(",", paste0(tip.aes, " = tip.", tip.aes), collapse = "")
    } else {
        tip.aes.txt = ""
    }

    tip.fix <- names(is.aes)[!is.aes]
    if (length(tip.fix)) {
        tip.fix.txt <- paste(",", paste0(tip.fix, " = tip.", tip.fix), collapse = "")
    } else {
        tip.fix.txt = ""
    }

    if (length(tip.fix) + length(tip.aes)) {
        p <- p + eval(parse(text = paste("geom_point(data = df.tip, aes_string(x = \"x\", y = \"y\"", tip.aes.txt,
                            ")", tip.fix.txt, ")")))
    }

    if ("color" %in% tip.aes) {
        nMaxCol <- brewer.pal.info[color.palette, "maxcolors"]
        tmp <- df.tip[, tip.color]

        ## if (is.discrete(tmp)) {
        if (is.character(tmp) || is.factor(tmp)) {
            if (length(unique(tmp)) > nMaxCol) {
                warning(paste("too many tip colors for palette", color.palette, ": we use default palette instead"))
                x <- tryCatch(eval(parse(text='p + scale_color_discrete(tip.color)')))                
                if(!"error" %in% class(tryCatch(print(x), error= function(e) e))){                  
                  p <- p + scale_color_discrete(tip.color)
                }                                                                     
            } else {
              x <- tryCatch(eval(parse(text='p + scale_color_brewer(tip.color, palette = color.palette)')))                
              if(!"error" %in% class(tryCatch(print(x), error= function(e) e))){                
                p <- p + scale_color_brewer(tip.color, palette = color.palette)
              }else{
                warning(paste("To avoid error 'Continuous value supplied to discrete scale' resulting from", color.palette, "we use default palette instead"))
              }                  
            }
        } else {
          x <- tryCatch(eval(parse(text='p + scale_color_gradientn(tip.color, colours = brewer.pal(nMaxCol, color.palette))')))                
          if(!"error" %in% class(tryCatch(print(x), error= function(e) e))){            
            p <- p + scale_color_gradientn(tip.color, colours = brewer.pal(nMaxCol, color.palette))
          }else{
            warning(paste("To avoid error 'Discrete value supplied to continuous scale' resulting from", color.palette, "we use default palette instead"))
          }            
        }
    }

    p <- p + theme(legend.position = legend.position)

    print(p)

    theme_set(theme.old)

    return(p)
}

