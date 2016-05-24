## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

plot.settings_handler <- function(widget, playState)
{
    if (playState$is.lattice) {
        latticeSettingsGUI(playState)
    } else {
        basePlotSettingsGUI(playState)
    }
}

latticeSettingsGUI <- function(playState)
{
    dialog <- gwindow(title="Plot settings")
    wingroup <- ggroup(horizontal=FALSE, container = dialog)
    tabs <- gnotebook(container = wingroup)
    wid <- list()

    origCall <- playState$call
    makeScalesArgAList(playState)

    trell <- playState$trellis
    isLatt3D <- !is.null(trell$panel.args.common$scales.3d)
    isMulti <- (prod(dim(trell)) > 1)
    scaleNames <- c("x", "y")
    if (isLatt3D) scaleNames <- c("x", "y", "z")
    relations <- c("same", "free", "sliced")
    aspects <- c('"fill"', '"iso"', '"xy"', '0.5', '1', '2')
    indexConds <-
        c("",
          "function(x, y) median(x, na.rm=TRUE)",
          "function(x, y) median(y, na.rm=TRUE)",
          "function(x, y) mean(x, na.rm=TRUE)",
          "function(x, y) mean(y, na.rm=TRUE)",
          "function(x, y) mean(y - x, na.rm=TRUE)",
          "function(x, y) sum(complete.cases(x))",
          "function(x, y) sum(complete.cases(x, y))")

    setArg <- function(h, ...) {
        target <- parse(text=h$action)[[1]]
        val <- svalue(h$obj)
        ## NAs possible from coerce.with = as.numeric
        if ((length(val) == 1) && is.na(val))
            val <- NULL
        callArg(playState, target) <- val
    }

    parseArg <- function(txt) {
        val <- parse(text = txt)
        if (length(val) == 0) return(NULL)
        val[[1]]
    }

    setArgTitle <- function(h, ...) {
        nm <- h$action
        val <- svalue(wid[[nm]])
        if (val == "") val <- NULL
        isExpr <- svalue(wid[[paste(nm, "expr", sep=".")]])
        if (!is.null(val) && isExpr)
            val <- parse(text=val, srcfile=NULL)
        ## insert explicit NULL, do not just omit the argument
        ## (to over-ride default labels)
        if (is.null(val)) val <- expression(NULL)
        callArg(playState, nm) <- val
    }

    setArgLayout <- function(h, ...) {
        val <- c(svalue(wid$layout.cols),
                 svalue(wid$layout.rows))
        if (any(is.na(val))) val <- NULL
        callArg(playState, "layout") <- val
    }

    setArgTck <- function(h, ...) {
        w <- h$action
        target <- paste("scales", w, "tck", sep="$")
        target <- parse(text=target)[[1]]
        val <- svalue(wid[[w]]$tck)
        if (is.na(val)) val <- NULL
        if (!is.null(val) && !svalue(wid[[w]]$ticks.opp))
            val[2] <- 0
        callArg(playState, target) <- val
    }

    setAxisLabels <- function(h, ...) {
        w <- h$action
        ## get this from playState again; might have changed
        trell <- playState$trellis
        if (isLatt3D) {
            w.scales <- trell$panel.args.common$scales.3d
            w.scales <- w.scales[[paste(w, ".scales", sep = "")]]
        } else {
            w.scales <- trell[[paste(w, ".scales", sep = "")]]
        }
        w.limits <- trell[[paste(w, ".limits", sep="")]]
        isFactorScale <-
            if (is.list(w.limits)) {
                is.character(w.limits[[1]]) || is.expression(w.limits[[1]])
            } else {
                is.character(w.limits) || is.expression(w.limits)
            }
        ## callTxt will be source code to set objects 'labels' and
        ## optionally 'at', for the corresponding entries in 'scales'
        callTxt <- character()
        if (isFactorScale) {
            ## factor scale, only need 'labels'
            origVal <- w.limits
            curVal <- w.scales$labels
            if (is.logical(curVal)) {
                curVal <- origVal
            }
            callTxt <- c("## Default axis labels: ", #
                         paste("#",                  #
                               deparse(call("<-", quote(labels), origVal)) ),
                         "## Custom axis labels: ", #
                         deparse(call("<-", quote(labels), curVal)) )
        } else {
            ## numeric scale, need 'at' and optionally 'labels'
            curAt <- w.scales$at
            curLabs <- w.scales$labels
            prettyInside <- function(x, n) {
                p <- pretty(x, n)
                p[(min(x) <= p) & (p <= max(x))]
            }
            if (is.list(w.limits)) {
                origAt <- lapply(w.limits, prettyInside, n = w.scales$tick.number)
                origLabs <- lapply(origAt, as.character)
            } else {
                origAt <- prettyInside(w.limits, n = w.scales$tick.number)
                origLabs <- as.character(origAt)
            }
            if (is.logical(curAt)) {
                curAt <- origAt
                origAt <- NULL
            }
            if (!is.null(origAt))
                callTxt <- c("## Default axis tick/label locations: ", #
                             paste("#", #
                                   deparse(call("<-", quote(at), origAt)) ))
            callTxt <- c(callTxt,
                         "## Custom axis tick/label locations: ", #
                         deparse(call("<-", quote(at), curAt)) )
            if (!is.logical(curLabs)) {
                callTxt <- c(callTxt,
                             deparse(call("<-", quote(labels), curLabs)))
            } else {
                callTxt <- c(callTxt,
                             "## Corresponding labels, optional:", #
                             paste("#",                            #
                                   deparse(call("<-", quote(labels), origLabs)) ))
            }
        }
        callTxt <- paste(callTxt, collapse = "\n")
        repeat {
            newTxt <- guiTextInput(callTxt, title = "Edit axis labels",
                                   prompt = "To reset: delete everything and press OK",
                                   accepts.tab = FALSE)
            if (is.null(newTxt)) break
            callTxt <- newTxt
            tmp <- tryCatch(parse(text = callTxt), error = error_handler)
            ## check whether there was a syntax error
            if (!inherits(tmp, "error")) {
                ## if more than one call, wrap them in braces
                tmpEnv <- new.env()
                result <- tryCatch(eval(tmp, tmpEnv), error = error_handler)
                if (!inherits(result, "error")) {
                    atTarget <- substitute(scales$w$at, list(w = as.symbol(w)))
                    labTarget <- substitute(scales$w$labels, list(w = as.symbol(w)))
                    callArg(playState, atTarget) <- tmpEnv[["at"]]
                    callArg(playState, labTarget) <- tmpEnv[["labels"]]
                    ## redraw (preview)
                    playReplot(playState)
                    showingPreview <<- TRUE
                    break
                }
            }
        }
    }

    ## (NEW TAB)
    basicsTab <- ggroup(horizontal=FALSE, container = tabs,
                        label = "Basic settings")

    ## TITLES
    labgroup <- gframe("Titles", horizontal=FALSE, container=basicsTab)
    lay <- glayout(spacing = 2, container=labgroup)
    rownum <- 1
    titleNames <- c("main", "sub", "xlab", "ylab")
    if (isLatt3D) titleNames <- c(titleNames, "zlab")
    for (nm in titleNames) {
        ## lattice titles can be: grob / list / vector
        argVal <- trell[[nm]]
        if (inherits(argVal, "grob")) {
            argVal <- toString(argVal)
        }
        if (is.list(argVal)) {
            argValMatch <- match.call(function(label, ...) NA,
                                      as.call(c(quote(foo), argVal)))
            argVal <- argValMatch$label
        }
        if (is.character(argVal) && (length(argVal) > 1)) {
            argVal <- as.expression(argVal)
        }
        isExpr <- is.language(argVal)
        if (isExpr) {
            argVal <- as.expression(argVal)
            ## multiple expressions in this format are parsed correctly
            argVal <- paste(sapply(argVal, deparseOneLine), collapse="; ")
        }
        nm.expr <- paste(nm, "expr", sep=".")
        lay[rownum, 1] <- nm
        lay[rownum, 2] <- wid[[nm]] <-
            gedit(toString(argVal), width=60, container = lay,
                  handler = setArgTitle, action = nm)
        ## plotmath option
        lay[rownum, 3] <- wid[[nm.expr]] <-
            gcheckbox("plotmath", checked=isExpr, container = lay,
                      handler = setArgTitle, action = nm)
        rownum <- rownum + 1
    }
    visible(lay) <- TRUE

    if (isLatt3D) {
        tmp <- ggroup(horizontal = TRUE, container = basicsTab)
        glabel("Axis labels offset distance: ", container = tmp)
        val <- trell$panel.args.common$scales.3d$x.scales$distance
        wid$scales.distance <-
            gedit(toString(val), width = 4, container = tmp,
                  handler = setArg, coerce.with = as.numeric,
                  action = "scales$distance")
    }

    ## STRIPS
    stripgroup <- gframe("Strips", horizontal=FALSE, container=basicsTab)
    tmp <- ggroup(horizontal=TRUE, container = stripgroup)
    wid$strip <- gcheckbox("Show top strips", container = tmp,
                           checked = !identical(trell$strip, FALSE),
                           handler = setArg, action = "strip")
    wid$strip.left <- gcheckbox("Show left strips", container = tmp,
                                checked = !identical(trell$strip.left, FALSE),
                                handler = setArg, action = "strip.left")
    ## TODO: show names / show levels
    ## TODO: edit strip text... (only for one conditioning variable)
    ## strip.custom(factor.levels=c(...), strip.names=FALSE, strip.levels=TRUE)
    wid$strip.abbrev <- gcheckbox("Abbreviate text", container = stripgroup,
                                  checked = isTRUE(trell$par.strip.text$abbreviate),
                                  handler = setArg, action = "par.strip.text$abbreviate")

    ## LAYOUT
    if (isMulti) {
        layoutgroup <- gframe("Layout", horizontal = FALSE, container = basicsTab)
        tmp <- ggroup(horizontal = TRUE, container = layoutgroup)
        val <- if (!is.null(trell$layout)) trell$layout else c(1, 1, 1)
        glabel("Layout per page: ", container = tmp)
        wid$layout.cols <- gspinbutton(from = 0, to = 16, value = val[1],
                                       handler = setArgLayout, container = tmp)
        glabel("columns, ", container = tmp)
        wid$layout.rows <- gspinbutton(from = 1, to = 16, value = val[2],
                                       handler = setArgLayout, container = tmp)
        glabel("rows. ", container = tmp)
        wid$as.table <-
            gcheckbox("as table", container = tmp,
                      checked = trell$as.table,
                      handler = setArg, action = "as.table")
        ## index.cond
        tmp <- ggroup(horizontal=TRUE, container = layoutgroup)
        glabel("Order panels by (index.cond): ", container = tmp)
        wid$index.cond <-
            gdroplist(indexConds, container = tmp, selected = 0,
                      editable = TRUE, coerce.with = parseArg,
                      handler = setArg, action = "index.cond")
        val <- callArg(playState, "index.cond", eval = FALSE)
        if (!is.null(val))
            svalue(wid$index.cond) <- deparseOneLine(val)
    }

    ## ASPECT
    aspectgroup <- gframe("Aspect ratio", horizontal = FALSE, container = basicsTab)
    tmp <- ggroup(horizontal=TRUE, container = aspectgroup)
    glabel("Panel aspect ratio y/x: ", container = tmp)
    val <- callArg(playState, "aspect")
    if (isLatt3D || is.null(val)) {
        val <- round(trell$aspect.ratio, 3)
        if (trell$aspect.fill) val <- "fill"
    }
    wid$aspect <- gdroplist(aspects, container = tmp,
                            editable = TRUE,
                            handler = setArg, coerce.with = parseArg,
                            action = if (isLatt3D) "panel.aspect" else "aspect")
    svalue(wid$aspect) <- deparse(val)
    ## TODO: aspect 3D?


    ## (NEW TAB)
    scalesTab <- ggroup(horizontal=FALSE, container = tabs,
                        label="Scales")

    ## SCALES
    isFactorScale <- list()
    isFactorScale$x <-
        if (is.list(trell$x.limits)) {
            is.character(trell$x.limits[[1]])
        } else {
            is.character(trell$x.limits)
        }
    isFactorScale$y <-
        if (is.list(trell$y.limits)) {
            is.character(trell$y.limits[[1]])
        } else {
            is.character(trell$y.limits)
        }
    isFactorScale$z <- FALSE
    lay <- glayout(spacing = 1, container = scalesTab)
    col <- 2
    firstCol <- TRUE
    for (w in scaleNames) {
        wid[[w]] <- list()
        if (isLatt3D) {
            w.scales <- trell$panel.args.common$scales.3d
            w.scales <- w.scales[[paste(w, ".scales", sep = "")]]
        } else {
            w.scales <- trell[[paste(w, ".scales", sep = "")]]
        }
        row <- 1
        lay[1, col] <- glabel(paste("<b>", w, "-axis</b>", sep=""),
                              markup = TRUE, container = lay)
        if (isMulti) {
            row <- row + 1
            if (firstCol) lay[row, 1] <- "relation:"
            lay[row, col] <- wid[[w]]$relation <-
                gdroplist(relations, container = lay,
                          selected = which(w.scales$relation == relations),
                          handler = setArg,
                          action = paste("scales", w, "relation", sep="$"))
        }
        row <- row + 1
        lay[row, col] <- wid[[w]]$draw <-
            gcheckbox("draw axis", container = lay,
                      checked = w.scales$draw,
                      handler = setArg,
                      action = paste("scales", w, "draw", sep="$"))
        if (isLatt3D) {
            row <- row + 1
            lay[row, col] <- wid[[w]]$arrows <-
                gcheckbox("arrows only", container = lay,
                          checked = isTRUE(w.scales$arrows),
                          handler = setArg,
                          action = paste("scales", w, "arrows", sep="$"))
        }
        row <- row + 1
        lay[row, col] <- wid[[w]]$log <-
            gcheckbox("logarithmic", container = lay,
                      checked = !identical(w.scales$log, FALSE),
                      handler = setArg,
                      action = paste("scales", w, "log", sep="$"))
        row <- row + 1
        lay[row, col] <- wid[[w]]$axs <-
            gcheckbox("padding", container = lay,
                      checked = (w.scales$axs == "r"),
                      handler = function(h, ...) {
                          target <- parse(text=h$action)[[1]]
                          val <- if (svalue(h$obj)) "r" else "i"
                          callArg(playState, target) <- val
                      },
                      action = paste("scales", w, "axs", sep="$"))
        row <- row + 1
        if (firstCol) lay[row, 1] <- "~ num. ticks: "
        lay[row, col] <- wid[[w]]$tick.number <-
            gedit(toString(w.scales$tick.number), width = 4, container = lay,
                  coerce.with = as.numeric, handler = setArg,
                  action = paste("scales", w, "tick.number", sep="$"))
        enabled(wid[[w]]$tick.number) <- !isFactorScale[[w]]
        row <- row + 1
        if (firstCol) lay[row, 1] <- "tick length: "
        lay[row, col] <- wid[[w]]$tck <-
            gedit(toString(w.scales$tck[1]), width = 4, container = lay,
                  coerce.with = as.numeric, handler = setArgTck,
                  action = w)
        if (!isLatt3D) {
            row <- row + 1
            lay[row, col] <- wid[[w]]$ticks.opp <-
                gcheckbox("ticks on opp. side", container = lay,
                          checked = (w.scales$tck[2] != 0),
                          handler = setArgTck,
                          action = w)
        }
        row <- row + 1
        if (firstCol) lay[row, 1] <- "label side:"
        alternList <- list(none = 0, standard = 1, opposite = 2,
                           alternating = c(1, 2), both = 3)
        whichAltern <- which(sapply(alternList, identical, w.scales$alternating))
        if (length(whichAltern) == 0) whichAltern <- 0
        lay[row, col] <- wid[[w]]$alternating <-
            gdroplist(names(alternList), container = lay,
                      selected = whichAltern,
                      handler = setArg,
                      coerce.with = function(nm) alternList[[nm]],
                      action = paste("scales", w, "alternating", sep="$"))
        row <- row + 1
        if (firstCol) lay[row, 1] <- "rotation:"
        lay[row, col] <- wid[[w]]$rot <-
            gedit(as.numeric(w.scales$rot[1]), width = 4, container = lay,
                  coerce.with = as.numeric, handler = setArg,
                  action = paste("scales", w, "rot", sep="$"))
        row <- row + 1
        lay[row, col] <- wid[[w]]$abbreviate <-
            gcheckbox("abbreviate labels", container = lay,
                      checked = w.scales$abbreviate,
                      handler = setArg,
                      action = paste("scales", w, "abbreviate", sep="$"))
        row <- row + 1
        ## specify axis labels...
        lay[row, col] <- wid[[w]]$axis.labels <-
            gbutton("set labels...", container = lay,
                    handler = setAxisLabels, action = w)

        ## TODO: custom axis components?

        ## next column
        col <- col + 1
        firstCol <- FALSE
    }
    visible(lay) <- TRUE

    svalue(tabs) <- 1
    showingPreview <- FALSE

    ok_handler <- function(h, ...)
    {
        playReplot(playState)
        showingPreview <<- TRUE
        if (h$action == "apply") return()
        dispose(h$obj)
        playState$win$present()
    }

    buttgroup <- ggroup(container=wingroup)
    addSpring(buttgroup)
    okbutt <- gbutton("OK", handler=ok_handler,
                      action="ok", container=buttgroup)
    prebutt <- gbutton("Apply", handler=ok_handler,
                       action="apply", container=buttgroup)
    canbutt <- gbutton("Cancel", container=buttgroup,
                       handler = function(h, ...) {
                           dispose(h$obj)
                           playState$call <- origCall
                           if (showingPreview)
                               playReplot(playState)
                       })
    size(okbutt) <- size(prebutt) <- size(canbutt) <- c(80, 30)
}


basePlotSettingsGUI <- function(playState)
{
    dialog <- gwindow(title="Plot settings")
    wingroup <- ggroup(horizontal=FALSE, container = dialog)
    tabs <- gnotebook(container=wingroup)
    wid <- list()
    ## keeps track of changed text fields during editing
    needUpdate <- list()
    setNeedUpdate <- function(h, ...)
        needUpdate[[h$action]] <<- TRUE

    origCall <- playState$call

    setArg <- function(h, ...) {
        target <- parse(text=h$action)[[1]]
        val <- svalue(h$obj)
        if (is.na(val)) val <- NULL
        callArg(playState, target) <- val
    }

    setArgTitle <- function(h, ...) {
        nm <- h$action
        val <- svalue(wid[[nm]])
        if (val == "") val <- NULL
        isExpr <- svalue(wid[[paste(nm, "expr", sep=".")]])
        if (!is.null(val) && isExpr)
            val <- parse(text=val, srcfile=NULL)
        ## insert explicit NULL, do not just omit the argument
        ## (to over-ride default labels)
        if (is.null(val)) val <- expression(NULL)
        callArg(playState, nm) <- val
    }

    setMar <- function(h, ...) {
        val <- c(svalue(wid$mar.bottom),
                 svalue(wid$mar.left),
                 svalue(wid$mar.top),
                 svalue(wid$mar.right))
        par(mar = val)
    }

    setType <- function(h, ...) {
        val <- NULL
        if (svalue(wid$points)) {
            val <- "p"
            if (svalue(wid$lines)) val <- "b"
        } else {
            if (svalue(wid$lines)) val <- "l"
            if (svalue(wid$droplines)) val <- "h"
        }
        callArg(playState, "type") <- val
    }

    setLog <- function(h, ...) {
        val <- paste(c(if (svalue(wid$x$log)) "x",
                       if (svalue(wid$y$log)) "y"),
                     collapse="")
        callArg(playState, "log") <- val
    }

    setLab <- function(h, ...) {
        val <- c(svalue(wid$x$tick.number),
                 svalue(wid$y$tick.number),
                 par("lab")[3])
        callArg(playState, "lab") <- val
    }

    ## convenience extractor
    arg <- function(x) callArg(playState, x)

    ## (NEW TAB)
    basicsTab <- ggroup(horizontal=FALSE, container = tabs,
                        label = "Basic settings")

    ## TITLES
    labgroup <- gframe("Titles", horizontal=FALSE, container=basicsTab)
    lay <- glayout(spacing = 1, container=labgroup)
    rownum <- 1
    titleNames <- c("main", "sub", "xlab", "ylab")
    for (nm in titleNames) {
        ## (base graphics etc) just a named argument
        argVal <- callArg(playState, nm)
        isExpr <- is.language(argVal)
        if (isExpr) {
            argVal <- as.expression(argVal)
            ## multiple expressions in this format are parsed correctly
            argVal <- paste(sapply(argVal, deparseOneLine), collapse="; ")
        }
        nm.expr <- paste(nm, "expr", sep=".")
        lay[rownum, 1] <- nm
        lay[rownum, 2] <- wid[[nm]] <-
            gedit(toString(argVal), width=60, container = lay,
                  handler = setArgTitle, action = nm)
         ## plotmath option
        lay[rownum, 3] <- wid[[nm.expr]] <-
            gcheckbox("plotmath", checked=isExpr, container = lay,
                      handler = setArgTitle, action = nm)
        rownum <- rownum + 1
    }
    visible(lay) <- TRUE

    ## MARGINS
    marginsgroup <- gframe("Margins", horizontal=FALSE, container=basicsTab)
    lay <- glayout(spacing = 1, container=marginsgroup)
    mar <- par("mar")
    lay[2,1] <- "Margins (lines): "
    lay[1,2] <- "bottom"
    lay[2,2] <- wid$mar.bottom <-
            gedit(toString(mar[1]), width = 4, container = lay,
                  handler = setMar, coerce.with = as.numeric)
    lay[1,3] <- "left"
    lay[2,3] <- wid$mar.left <-
            gedit(toString(mar[2]), width = 4, container = lay,
                  handler = setMar, coerce.with = as.numeric)
    lay[1,4] <- "top"
    lay[2,4] <- wid$mar.top <-
            gedit(toString(mar[3]), width = 4, container = lay,
                  handler = setMar, coerce.with = as.numeric)
    lay[1,5] <- "right"
    lay[2,5] <- wid$mar.right <-
            gedit(toString(mar[4]), width = 4, container = lay,
                  handler = setMar, coerce.with = as.numeric)
    visible(lay) <- TRUE

    ## TYPE
    typegroup <- gframe("Plot type", horizontal=FALSE, container=basicsTab)
    arg_type <- callArg(playState, "type")
    hasPoints <- (is.null(arg_type) || any(c("p","b","o") %in% arg_type))
    hasLines <- any(c("l","b","o") %in% arg_type)
    hasDrops <- any("h" %in% arg_type)
    tmp <- ggroup(horizontal = TRUE, container=typegroup)
    wid$points <- gcheckbox("Points", checked=hasPoints, container = tmp,
                            handler = setType)
    wid$lines <- gcheckbox("Lines", checked=hasLines, container = tmp,
                           handler = setType)
    wid$droplines <- gcheckbox("Drop lines", checked=hasDrops, container = tmp,
                               handler = setType)


    ## (NEW TAB)
    scalesTab <- ggroup(horizontal=FALSE, container = tabs,
                        label="Scales")

    ## ASPECT
    aspectgroup <- gframe("Aspect ratio", horizontal = FALSE, container = scalesTab)
    tmp <- ggroup(horizontal=FALSE, container = aspectgroup)
    wid$asp1 <- gcheckbox("Isometric scales", container = aspectgroup,
                          checked = isTRUE(arg("asp") == 1),
                          handler = function(h, ...) {
                              val <- if (svalue(h$obj)) 1 else NULL
                              callArg(playState, "asp") <- val
                          })
    wid$ptys <- gcheckbox("Square plot region", container = aspectgroup,
                          checked = (par("pty") == "s"),
                          handler = function(h, ...) {
                              val <- if (svalue(h$obj)) "s" else "m"
                              par(pty = val)
                          })

    ## SCALES
    lay <- glayout(spacing = 1, container = scalesTab)
    col <- 2
    firstCol <- TRUE
    for (w in c("x", "y")) {
        wid[[w]] <- list()
        row <- 1
        lay[1, col] <- glabel(paste("<b>", w, "-axis</b>", sep=""),
                              markup = TRUE, container = lay)
        row <- row + 1
        w.axt <- paste(w, "axt", sep="")
        lay[row, col] <- wid[[w]]$draw <-
            gcheckbox("draw axis", container = lay,
                      checked = !(any(arg("axes") == FALSE) ||
                                  any(arg(w.axt) == "n")),
                      handler = function(h, ...) {
                          val <- if (svalue(h$obj)) NULL else "n"
                          callArg(playState, h$action) <- val
                      }, action = w.axt)
        row <- row + 1
        lay[row, col] <- wid[[w]]$log <-
            gcheckbox("logarithmic", container = lay,
                      checked = par(paste(w, "log", sep="")),
                      handler = setLog)
        row <- row + 1
        w.axs <- paste(w, "axs", sep="")
        lay[row, col] <- wid[[w]]$axs <-
            gcheckbox("padding", container = lay,
                      checked = !any(arg(w.axs) == "i"),
                      handler = function(h, ...) {
                          val <- if (svalue(h$obj)) NULL else "i"
                          callArg(playState, h$action) <- val
                      }, action = w.axs)
        row <- row + 1
        if (firstCol) lay[row, 1] <- "~num. ticks:"
        lay[row, col] <- wid[[w]]$tick.number <-
            gedit(toString(par("lab")[1]), width = 4, container = lay,
                  handler = setLab, coerce.with = as.numeric)
        row <- row + 1
        if (firstCol) lay[row, 1] <- "tick length:"
        lay[row, col] <- wid[[w]]$tcl <-
            gedit(toString(par("tcl")), width = 4, container = lay,
                  handler = setArg, coerce.with = as.numeric,
                  action = "tcl")
        ## next column
        col <- col + 1
        firstCol <- FALSE
    }
    visible(lay) <- TRUE

    svalue(tabs) <- 1
    showingPreview <- FALSE

    ok_handler <- function(h, ...)
    {
        playReplot(playState)
        showingPreview <<- TRUE
        if (h$action == "apply") return()
        dispose(h$obj)
        playState$win$present()
    }

    buttgroup <- ggroup(container=wingroup)
    addSpring(buttgroup)
    okbutt <- gbutton("OK", handler=ok_handler,
                      action="ok", container=buttgroup)
    prebutt <- gbutton("Apply", handler=ok_handler,
                       action="apply", container=buttgroup)
    canbutt <- gbutton("Cancel", container=buttgroup,
                       handler=function(h, ...) {
                           dispose(h$obj)
                           playState$call <- origCall
                           if (showingPreview)
                               playReplot(playState)
                       })
    size(okbutt) <- size(prebutt) <- size(canbutt) <- c(80, 30)
}

