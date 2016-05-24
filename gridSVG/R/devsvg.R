
# Functions to create an SVG graphics device object, complete with
# "methods" for performing all necessary graphical operations

# This is designed to foreshadow the time when graphics devices
# in R (or at least in grid) are R objects and graphics functions
# include the device as an argument (i.e., no longer have
# the notion of graphics always going to the "current device")

# This will not be called in this way yet (instead, I will just
# be running down the grid display list and calling appropriate
# methods from that, BUT I thought it was worth designing for the
# future anyway.

# In another forward-looking move, I will create the device class
# and methods for it using S4 methods

#################
# Utility functions
#################

# Any non-grid parameters is let through untouched
# BUT code later in svg.R will complain about things that are
# not SVG parameters
# NOTE that 'cex'/'lex' have been incorporated into 'fontsize'/'lwd'
# and removed by this point
devParNameToSVGStyleName <- function(name) {
  switch(name,
         col="stroke",
         colAlpha="stroke-opacity",
         fill="fill",
         fillAlpha="fill-opacity",
         fontweight="font-weight",
         fontfamily="font-family",
         fontstyle="font-style",
         fontsize="font-size",
         alpha="opacity",
         lty="stroke-dasharray",
         lwd="stroke-width",
         lineend="stroke-linecap",
         linejoin="stroke-linejoin",
         linemitre="stroke-miterlimit",
         name)
}

# R lwd is in points, pixels or 1/96 inches
# However, most (perhaps all?) devices use 1/96 for their
# definition of an 'lwd', so use that.
devLwdToSVG <- function(lwd, dev) {
    round(lwd/96 * dev@res, 2)
}

# An R lty has to become an SVG stroke-dasharray
# This is going to be imperfect (to say the least)
devLtyToSVG <- function(lty, lwd, dev) {
    # Convert lty to numeric vec
    numlty <- switch(lty,
                     solid=0,
                     # These numbers taken from ?par
                     dashed=c(4, 4),
                     dotted=c(1, 3),
                     dotdash=c(1, 3, 4, 3),
                     longdash=c(7, 3),
                     twodash=c(2, 2, 6, 2),
                     # Otherwise we're a hex string
                     as.numeric(as.hexmode(strsplit(lty, "")[[1]])))
    # Scale by lwd
    scaledlty <- numlty * lwd
    # Convert to SVG stroke-dasharray string
    paste(ifelse(scaledlty == 0,
                 "none",
                 round(scaledlty/96 * dev@res, 2)),
          collapse=",")
}

devColToSVG <- function(col) {
  if (length(col) > 1)
    warning("Only first colour used")
  if (is.numeric(col) && col == 0)
      col <- "transparent"
  # Handle "transparent" as a special case
  if (col == "transparent")
      "none"
  else
      paste("rgb(", paste(col2rgb(col), collapse=","), ")", sep="")
}

devColAlphaToSVG <- function(colAlpha) {
    round(colAlpha/255, 2)
}

devFontSizeToSVG <- function(fontsize, dev) {
    round(fontsize/72 * dev@res, 2)
}

devLineJoinToSVG <- function(linejoin, dev) {
    # Only need to change spelling of mitre, SVG takes american form
    if (linejoin == "mitre")
        "miter"
    else
        linejoin
}

devFontFaceToSVG <- function(fontface) {
    # CSS uses two different properties to configure the appearance of a font
    # Setting defaults to CSS defaults
    fontWeightCSS <- "normal"
    fontStyleCSS <- "normal"

    if (is.numeric(fontface)) {
        if (fontface == 1) {
            # plain
            fontWeightCSS <- "normal"
            fontStyleCSS <- "normal"
        }

        if (fontface == 2) {
            # bold
            fontWeightCSS <- "bold"
            fontStyleCSS <- "normal"
        }

        if (fontface == 3) {
            # italic
            fontWeightCSS <- "normal"
            fontStyleCSS <- "italic"
        }

        if (fontface == 4) {
            # bold italic
            fontWeightCSS <- "bold"
            fontStyleCSS <- "italic"
        }
    }

    if (is.character(fontface)) {
        if (fontface == "plain") {
            fontWeightCSS <- "normal"
            fontStyleCSS <- "normal"
        }

        if (fontface == "bold") {
            fontWeightCSS <- "bold"
            fontStyleCSS <- "normal"
        }

        if (fontface == "italic") {
            fontWeightCSS <- "normal"
            fontStyleCSS <- "italic"
        }

        if (fontface == "oblique") {
            fontWeightCSS <- "normal"
            fontStyleCSS <- "oblique"
        }

        if (fontface == "bold.italic") {
            fontWeightCSS <- "bold"
            fontStyleCSS <- "italic"
        }
    }

    list(fontweight=fontWeightCSS,
         fontstyle=fontStyleCSS)
}

getSVGFonts <- function() {
    get("gridSVG.fonts", envir = .gridSVGEnv)
}

setSVGFonts <- function(fontStacks) {
    if (! all(names(fontStacks) == c("sans", "serif", "mono")))
        stop("Font settings must have fonts available for 'sans', 'serif' and 'mono'.")

    # Need to ensure that basic font fallbacks are available and
    # are placed at the end of each of the font stacks.
    if (! "sans-serif" %in% fontStacks$sans) {
        fontStacks$sans <- c(fontStacks$sans, "sans-serif")
    } else if (tail(fontStacks$sans, n = 1) != "sans-serif") {
        ind <- which(fontStacks$sans == "sans-serif")
        cleanedSans <- fontStacks$sans[-ind]
        fontStacks$sans <- c(cleanedSans, "sans-serif")
    }
    if (! "serif" %in% fontStacks$serif) {
        fontStacks$serif <- c(fontStacks$serif, "serif")
    } else if (tail(fontStacks$serif, n = 1) != "serif") {
        ind <- which(fontStacks$serif == "serif")
        cleanedSerif <- fontStacks$serif[-ind]
        fontStacks$serif <- c(cleanedSerif, "serif")
    }
    if (! "monospace" %in% fontStacks$mono) {
        fontStacks$mono <- c(fontStacks$mono, "monospace")
    } else if (tail(fontStacks$mono, n = 1) != "monospace") {
        ind <- which(fontStacks$mono == "monospace")
        cleanedMono <- fontStacks$mono[-ind]
        fontStacks$mono <- c(cleanedMono, "monospace")
    }

    assign("gridSVG.fonts", fontStacks, envir = .gridSVGEnv)
}

# Setting default font stacks
sansFontStack <- c("Helvetica", "Arial", "FreeSans",
                   "Liberation Sans", "Nimbus Sans L", "sans-serif")
serifFontStack <- c("Times", "Times New Roman", "Liberation Serif",
                    "Nimbus Roman No9 L Regular", "serif")
monoFontStack <- c("Courier", "Courier New", "Nimbus Mono L", "monospace")
setSVGFonts(list(sans = sansFontStack,
                 serif = serifFontStack,
                 mono = monoFontStack))

fontStackFromFontFamily <- function(fontfamily, currentFonts) {
    if (fontfamily %in% c(currentFonts$sans, "sans"))
        "sans"
    else if (fontfamily %in% currentFonts$serif)
        "serif"
    else if (fontfamily %in% c(currentFonts$mono, "mono"))
        "mono"
    else 
        "unknown"
}

devFontFamilyToSVG <- function(fontfamily, dev) {
    currentFonts <- getSVGFonts()
    stackname <- fontStackFromFontFamily(fontfamily, currentFonts)

    if (stackname == "unknown") {
        if (nchar(fontfamily) > 0)
            # Assume font exists, but also assume sans-serif fallback
            fontstack <- c(fontfamily, currentFonts$sans)
        else 
            fontstack <- currentFonts$sans # Assuming a sans-serif font
    } else {
        fontstack <- currentFonts[[stackname]]
    }
    
    # Formatting the font stack for CSS
    fontStackCSS <- paste(fontstack, collapse=', ')

    # Returning the font stack
    fontStackCSS
}

devParToSVGPar <- function(name, par, dev) {
  if (is.null(par))
    "none"
  else {
      ifelse(is.na(par),
             "none",
             switch(name,
                    col=devColToSVG(par),
                    colAlpha=devColAlphaToSVG(par),
                    fill=devColToSVG(par),
                    fillAlpha=devColAlphaToSVG(par),
                    fontsize=devFontSizeToSVG(par, dev),
                    fontfamily=devFontFamilyToSVG(par, dev),
                    lwd=devLwdToSVG(par, dev),
                    linejoin=devLineJoinToSVG(par, dev),
                    # By default just pass through the actual value
                    # e.g., lty has already been converted at this point
                    par))
  }
}

devParToSVGStyle <- function(gp, dev) {
    if (is.null(gp))
        result <- svgStyle()
    else {
        result <- list()
        # convert "cex" into "fontsize"
        if ("cex" %in% names(gp)) {
            if ("fontsize" %in% names(gp))
                gp$fontsize <- (gp$fontsize * gp$cex)
            else
                gp$fontsize <- (get.gpar("fontsize")[[1]] * gp$cex)
            gp$cex <- NULL
        }
        # Do the same for "lex"
        if ("lex" %in% names(gp)) {
            if ("lwd" %in% names(gp))
                gp$lwd <- (gp$lwd * gp$lex)
            else
                gp$lwd <- (get.gpar("lwd")[[1]] * gp$lex)
            gp$lex <- NULL
        }
        # Just remove "lineheight" (this has already been incorporated
        # into text object information by this point)
        # Remove it so that it is not exported as SVG attribute
        gp$lineheight <- NULL
        # Scale lty by lwd
        if ("lty" %in% names(gp)) {
            if ("lwd" %in% names(gp)) {
                gp$lty <- devLtyToSVG(gp$lty, gp$lwd, dev)
            } else {
                gp$lty <- devLtyToSVG(gp$lty, 1, dev)
            }
        }
        # Font is an alias for fontface, set to fontface
        if ("font" %in% names(gp)) {
            gp$fontface <- gp$font
            gp$font <- NULL
        }
        # Split fontface into fontweight and fontstyle
        if ("fontface" %in% names(gp)) {
            svgFont <- devFontFaceToSVG(gp$fontface)
            gp$fontweight <- svgFont$fontweight
            gp$fontstyle <- svgFont$fontstyle
            gp$fontface <- NULL
        }
        for (i in names(gp))
            if (!is.na(devParNameToSVGStyleName(i)))
                result[[devParNameToSVGStyleName(i)]] <-
                    devParToSVGPar(i, gp[[i]], dev)
    }
    result
}

#################
# SVG Device Stuff
#################

setClass("svgDevice",
         representation("graphicsDevice",
                        res="numeric",
                        attrs="list",
                        links="character",
                        show="character",
                        # Object created by svgDevice() in svg.R
                        # has no S4 class yet
                        dev="ANY"))

setMethod("inchToDevX", signature(device="svgDevice"),
          function(x, device) {
            x * device@res
          })
          
setMethod("inchToDevY", signature(device="svgDevice"),
          function(x, device) {
            x * device@res
          })

setMethod("devArrow", signature(device="svgDevice"),
          function(arrow, gp, device) {
            # Angle is specified for the arrowhead in degrees, need radians
            ratAngle <- (pi / 180) * arrow$angle

            # We know the length, it is the hypotenuse, need to find the
            # length of the opposite line for the entire arrowhead, not
            # just one half
            midpoint <- sin(ratAngle) * arrow$length
            arrowWidth <- midpoint * 2
            xmult <- cos(ratAngle)
            arrowX <- xmult * arrow$length

            xs <- unit.c(unit(0, "inches"), arrowX, unit(0, "inches"))
            ys <- unit.c(unit(0, "inches"), midpoint, arrowWidth)
            x <- cx(xs, device)
            y <- cy(ys, device)

            svgMarker(x, y, arrow$type, arrow$ends, sign(xmult), arrow$name,
                      devParToSVGStyle(gp, device), device@dev)
          })
          
setMethod("devLines", signature(device="svgDevice"),
          function(lines, gp, device) {
            svgLines(lines$x, lines$y, lines$name, lines$arrow,
                     device@attrs, device@links, device@show,
                     devParToSVGStyle(gp, device), device@dev)
          })

setMethod("devPolygon", signature(device="svgDevice"),
          function(polygon, gp, device) {
            svgPolygon(polygon$x, polygon$y, polygon$name,
                       device@attrs, device@links, device@show,
                       devParToSVGStyle(gp, device), device@dev)
          })

setMethod("devPath", signature(device="svgDevice"),
          function(path, gp, device) {
            svgPath(path$x, path$y, path$rule, path$name,
                    device@attrs, device@links, device@show,
                    devParToSVGStyle(gp, device), device@dev)
          })

setMethod("devRaster", signature(device="svgDevice"),
          function(raster, gp, device) {
            svgRaster(raster$x, raster$y, raster$width, raster$height,
                      raster$angle, raster$datauri,
                      raster$name, raster$just, raster$vjust, raster$hjust,
                      listToSVGAttrib(raster$attributes), device@links,
                      device@show, devParToSVGStyle(gp, device), device@dev)
          })

setMethod("devRect", signature(device="svgDevice"),
          function(rect, gp, device) {
            svgRect(rect$x, rect$y, rect$width, rect$height, rect$angle,
                    rect$name,
                    device@attrs, device@links, device@show,
                    devParToSVGStyle(gp, device), device@dev)
          })

setMethod("devText", signature(device="svgDevice"),
          function(text, gp, device) {
              # SVG text will use fill, but fill has already been
              # set to col back in primToDev.text() in griddev.R
              svgText(text$x, text$y, text$text,
                      text$hjust, text$vjust, text$rot,
                      text$width, text$height, text$angle,
                      text$ascent, text$descent,
                      text$lineheight, text$charheight, text$fontheight,
                      text$fontfamily, text$fontface, text$name,
                      device@attrs, device@links, device@show,
                      devParToSVGStyle(gp, device), device@dev)
          })

setMethod("devCircle", signature(device="svgDevice"),
          function(circle, gp, device) {
            svgCircle(circle$x, circle$y, circle$r, circle$name,
                      device@attrs, device@links, device@show,
                      devParToSVGStyle(gp, device), device@dev)
          })

setMethod("devStartElement", signature(device="svgDevice"),
          function(element, gp, device) {
            # Ignore gp, complicates output
            svgStartElement(id = element$id,
                            classes = element$classes,
                            element = element$name,
                            attrs = element$attrs,
                            namespace = element$namespace,
                            namespaceDefinitions = element$namespaceDefinitions,
                            attributes = device@attrs,
                            links = device@links,
                            show = device@show,
                            svgdev = device@dev)
          })

setMethod("devEndElement", signature(device="svgDevice"),
          function(name, device) {
            svgEndElement(name, device@links, device@dev)
          })

setMethod("devTextNode", signature(device="svgDevice"),
          function(text, device) {
            svgTextNode(text$text, device@dev)
          })

setMethod("devStartClip", signature(device="svgDevice"),
          function(clip, gp, device) {
            svgClipPath(clip$name, clip$x, clip$y,
                        clip$width, clip$height, clip$angle,
                        device@dev)

            # Because of the fact that we never stop clipping until
            # we pop our current viewport, we need to store information
            # on how many times we have clipped.
            # This allows us to traverse back up the appropriate number
            # of SVG <g>s.
            cl <- get("contextLevels", envir = .gridSVGEnv)
            cl[length(cl)] <- cl[length(cl)] + 1
            assign("contextLevels", cl, envir = .gridSVGEnv)

            # Can hard-code 'clip' and 'coords' because we're always clipping
            # but we're not a viewport.
            # 'style' is always going to be NULL too.
            svgStartGroup(clip$name, clip=TRUE,
                          attributes=device@attrs,
                          links=device@links,
                          show=device@show,
                          style=devParToSVGStyle(gp, device),
                          coords = NULL,
                          classes = clip$classes,
                          svgdev=device@dev)
          })

setMethod("devStartClipPath", signature(device="svgDevice"),
          function(clippath, gp, device) {
            svgStartGrobClipPath(clippath$name, device@dev)
          })

setMethod("devEndClipPath", signature(device="svgDevice"),
          function(clippath, gp, device) {
            svgEndGrobClipPath(device@dev)
          })

setMethod("devStartClipPathGroup", signature(device="svgDevice"),
          function(clippath, gp, device) {
            svgStartGrobClipPathGroup(clippath$name, clippath$cp,
                                      clippath$classes, device@dev)

            # Because of the fact that we never stop clipping until
            # we pop our current viewport, we need to store information
            # on how many times we have clipped.
            # This allows us to traverse back up the appropriate number
            # of SVG <g>s.
            cl <- get("contextLevels", envir = .gridSVGEnv)
            cl[length(cl)] <- cl[length(cl)] + 1
            assign("contextLevels", cl, envir = .gridSVGEnv)
            # Also note the ID because we're pushing a context, makes it
            # easier to locate later
            assign("contextNames",
                   c(get("contextNames", envir = .gridSVGEnv), clippath$name),
                   envir = .gridSVGEnv)
          })

setMethod("devStartMask", signature(device="svgDevice"),
          function(mask, gp, device) {
            svgStartMask(mask$name, mask$x, mask$y, mask$width,
                         mask$height, device@dev)
          })

setMethod("devEndMask", signature(device="svgDevice"),
          function(mask, gp, device) {
            svgEndMask(device@dev)
          })

setMethod("devStartMaskGroup", signature(device="svgDevice"),
          function(mask, gp, device) {
            svgStartMaskGroup(mask$name, mask$mask,
                              mask$classes, device@dev)

            # Because of the fact that we never stop clipping until
            # we pop our current viewport, we need to store information
            # on how many times we have clipped.
            # This allows us to traverse back up the appropriate number
            # of SVG <g>s.
            cl <- get("contextLevels", envir = .gridSVGEnv)
            cl[length(cl)] <- cl[length(cl)] + 1
            assign("contextLevels", cl, envir = .gridSVGEnv)
            # Also note the ID because we're pushing a context, makes it
            # easier to locate later
            assign("contextNames",
                   c(get("contextNames", envir = .gridSVGEnv), mask$name),
                   envir = .gridSVGEnv)
          })

setMethod("devStartGroup", signature(device="svgDevice"),
          function(group, gp, device) {
              clip <- FALSE
              if (! is.null(group$clip)) {
                  if (group$clip) {
                      clip <- TRUE
                      svgClipPath(group$name, group$vpx, group$vpy,
                                  group$vpw, group$vph, group$angle,
                                  device@dev)
                  }
              }

            # If we're starting a VP, then allow for "contexts" to be
            # added to children of this VP. A context is a clip path
            # or mask. Coords are only present via VPs.
            if (! is.null(group$coords)) {

              assign("contextLevels",
                     c(get("contextLevels", envir = .gridSVGEnv), 0),
                     envir = .gridSVGEnv)
            }

            svgStartGroup(group$name, clip=clip,
                          attributes=device@attrs,
                          links=device@links,
                          show=device@show,
                          style=devParToSVGStyle(gp, device),
                          coords = group$coords,
                          classes = group$classes,
                          svgdev=device@dev)
          })

setMethod("devEndGroup", signature(device="svgDevice"),
          function(name, vp, device) {
            svgEndGroup(name, device@links, vp, device@dev)
          })

setMethod("devStartSymbol", signature(device="svgDevice"),
          function(pch, device) {
            svgStartSymbol(pch, device@dev)
          })

setMethod("devPoint", signature(device="svgDevice"),
          function(pch, device) {
            svgPoint(pch, device@dev)
          })

setMethod("devEndSymbol", signature(device="svgDevice"),
          function(device) {
            svgEndSymbol(device@dev)
          })

setMethod("devUseSymbol", signature(device="svgDevice"),
          function(point, gp, device) {
            svgUseSymbol(point$name, point$x, point$y, point$size, point$pch,
                         point$angle,
                         device@attrs, device@links, device@show,
                         devParToSVGStyle(gp, device), device@dev)
          })

setMethod("devClose", signature(device="svgDevice"),
          function(device) {
            svgClose(device@dev)
          })

#################
# User Functions
#################

openSVGDev <- function(name="Rplots.svg", width=6, height=6, res=NULL,
                       strict=TRUE, rootAttrs=NULL) {
    if (is.null(res))
        res <- par("cra")[1]/par("cin")[1]
        # par("cra")[2]/par("cin")[2]*height))
    
    new("svgDevice",
        width=width, height=height,
        res=res,
        dev=svgOpen(res*width, res*height, strict, rootAttrs))
}
                   


