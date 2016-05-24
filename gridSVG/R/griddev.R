
vpError <- function() {
  stop("vp should only be path")
}

# Functions to take a grid grob and call appropriate
# functions from dev.R to produce output on a device

# Each function has to convert locations and dimensions
# into device coordinates THEN call the dev.R function

# Convert a unit object to a value in "device" units
# The calls to convert*() are just to get 'valueOnly'
# from "inches" units.
cx <- function(x, dev) {
  inchToDevX(convertX(x, "inches", valueOnly=TRUE), dev)
}

cy <- function(x, dev) {
  inchToDevY(convertY(x, "inches", valueOnly=TRUE), dev)
}

cw <- function(x, dev) {
  inchToDevX(convertWidth(x, "inches", valueOnly=TRUE), dev)
}

ch <- function(x, dev) {
  inchToDevY(convertHeight(x, "inches", valueOnly=TRUE), dev)
}

# Convert a "distance" (e.g., a circle radius)
cd <- function(x, dev) {
  pmin(inchToDevX(convertWidth(x, "inches", valueOnly=TRUE), dev),
       inchToDevY(convertHeight(x, "inches", valueOnly=TRUE), dev))
}

# Create a full name for a sub-grob based on the name of a parent grob
subGrobName <- function(baseGrobName, subGrobName,
                        separator = getSVGoption("id.sep")) {
    paste(baseGrobName, subGrobName, sep=separator)
}

# Return the base grob name given the full name of a sub-grob
baseGrobName <- function(subGrobName, 
                         separator = getSVGoption("id.sep")) {
  splitName <- unlist(strsplit(subGrobName, separator, fixed = TRUE))
  grobName <- paste(splitName[-length(splitName)], collapse = separator)

  # Returning the base name
  grobName
}

prefixName <- function(name) {
  paste0(get("prefix", envir = .gridSVGEnv), name)
}

# Convert a gpar object to an device-neutral graphical parameter list
gparToDevPars <- function(gp) {
    # Split up col into col plus colAlpha
    if (!is.null(gp$col)) {
        if (is.numeric(gp$col) && gp$col == 0)
            gp$col <- "transparent"
        rgba <- col2rgb(gp$col, alpha=TRUE)
        gp$colAlpha <- rgba[4]
    }
    # Ditto fill
    if (!is.null(gp$fill)) {
        if (is.numeric(gp$fill) && gp$fill == 0)
            gp$fill <- "transparent"
        rgba <- col2rgb(gp$fill, alpha=TRUE)
        gp$fillAlpha <- rgba[4]
    }
    gp
}

# Repeats all elements in a gpar() so that it is fully defined for n values
expandGpar <- function(gp, n) {
    if (is.null(gp))
        return(gpar())
    # If there are actually gpar elements defined, repeat them
    if (length(gp) > 0) {
        for (i in 1:length(gp)) {
            gp[[i]] <- rep(gp[[i]], length.out = n)
        }
    }

    # Returning the gp
    gp
}

# Repeats all elements in an arrow() so that it is fully defined for n values
expandArrow <- function(arrow, n) {
    # If there is actually an arrow, repeat its components
    if (! is.null(arrow)) {
        for (i in 1:length(arrow)) {
            arrow[[i]] <- rep(arrow[[i]], length.out = n)
        }
    }

    # Returning the arrow
    arrow
}

# Converting locations and widths
locToInches <- function(x, y, dev) {
  # Convert x and y to inches
  x <- convertX(x, "inches", valueOnly=TRUE)
  y <- convertY(y, "inches", valueOnly=TRUE)
  # Transform to inches on device
  n <- max(length(x), length(y))
  loc <- cbind(rep(x, length=n),
               rep(y, length=n),
               rep(1, length=n)) %*% current.transform()
  x <- unit(loc[,1]/loc[,3], "inches")
  y <- unit(loc[,2]/loc[,3], "inches")
  list(x=x, y=y)
}

dimToInches <- function(w, h, dev) {
  # FIXME:  Doesn't handle rotated viewports!!
  w <- convertWidth(w, "inches")
  h <- convertHeight(h, "inches")
  list(w=w, h=h)
}

dToInches <- function(d, dev) {
  w <- convertWidth(d, "inches", valueOnly=TRUE)
  h <- convertHeight(d, "inches", valueOnly=TRUE)
  d <- unit(pmin(w, h), "inches")
  d
}

# Generate (left, bottom) from (x, y), (width, height), and justification
leftbottom <- function(x, y, width, height,
                       just, hjust, vjust, dev) {
    hjust <- resolveHJust(just, hjust)
    vjust <- resolveVJust(just, vjust)
    left <- unit(convertX(x, "inches", valueOnly=TRUE) -
                 convertWidth(hjust*width, "inches", valueOnly=TRUE),
                 "inches")
    bottom <- unit(convertY(y, "inches", valueOnly=TRUE) -
                   convertHeight(vjust*height, "inches", valueOnly=TRUE),
                   "inches")
    locToInches(left, bottom, dev)
}

# Generate hjust/vjust from just
justTohjust <- function(just) {
  if (length(just) > 1)
    just <- just[1]

  if (is.numeric(just)) {
    # Rounding to nearest of 0, 0.5, 1
    roundedJust <- round(2 * just) / 2
    # and clamped to 0 to 1
    if (roundedJust < 0)
        roundedJust <- 0
    if (roundedJust > 1)
        roundedJust <- 1
    switch(as.character(roundedJust),
           "0" = "left",
           "0.5" = "centre",
           "1" = "right")
  } else {
    if (is.na(match(just[1], c("left", "right"))))
      "centre"
    else
      just[1]
  }
}

justTovjust <- function(just) {
  if (length(just) > 1)
    just <- just[2]

  if (is.numeric(just)) {
    # Rounding to nearest of 0, 0.5, 1
    roundedJust <- round(2 * just) / 2
    # and clamped to 0 to 1
    if (roundedJust < 0)
        roundedJust <- 0
    if (roundedJust > 1)
        roundedJust <- 1
    switch(as.character(roundedJust),
           "0" = "bottom",
           "0.5" = "centre",
           "1" = "top")
  } else {
    if (is.na(match(just[1], c("top", "bottom"))))
      "centre"
    else
      just
  }
}

changedGPar <- function(startGP, endGP) {
    diffGP <- mapply(function(x, y) !isTRUE(all.equal(x, y)),
                     endGP, startGP)
    do.call("gpar", unclass(endGP)[diffGP])
}

# Enforce a 'vp' setting
# This could be a viewport (or vpTree or vpList or vpStack) OR a vpPath
# The general idea is to push or down the 'vp' slot
# THEN check how far down we have come
# IF we have come more than one level down then start a group for
# the appropriate number of parent viewports as well as the current viewport
# AND then do the corresponding number of end groups afterwards
startGroup <- function(vp, depth, dev) {
    if (depth > 1) {
        path <- upViewport(depth - 1, recording=FALSE)
        paths <- explode(path)
        for (i in paths) {
            parent <- current.viewport()
            parent$classes <- class(parent)
            devStartGroup(devGrob(parent, dev), gparToDevPars(parent$gp), dev)
            downViewport(i, recording=FALSE)
        }
    }
    vp$classes <- class(vp)
    devStartGroup(devGrob(vp, dev), gparToDevPars(vp$gp), dev)
}
enforceVP <- function(vp, dev) {
    depth <- 0
    if (!is.null(vp)) {
        if (!inherits(vp, "vpPath")) {
            pushViewport(vp, recording=FALSE)
            depth <- depth(vp)
        } else {
            depth <- downViewport(vp, recording=FALSE)
        }
        startGroup(current.viewport(), depth, dev)
    }
    depth
}
unwindVP <- function(vp, depth, dev) {
    if (depth > 0) {
        for (i in 1:depth)
            devEndGroup("", TRUE, dev)
        upViewport(depth, recording=FALSE)
    }
}

# Grob to SVG
# This mimics grid.draw()
# Push/down any viewports and then call primToDev() to produce SVG
grobToDev <- function(x, dev) {
  UseMethod("grobToDev", x)
}

grobToDev.default <- function(x, dev) {
  stop("We shouldn't be here!")
}

grobToDev.grob <- function(x, dev) {
  depth <- enforceVP(x$vp, dev)
  x$classes <- class(x)
  primToDev(x, dev)
  unwindVP(x$vp, depth, dev)
  progressStep("grob")
}

# grob to device grob
# This just converts a grid grob into a generic (bland) device grob
# (which is just a list of values)
devGrob <- function(x, dev) {
  UseMethod("devGrob")
}

devGrob.default <- function(x, dev) {
  list(name=x$name, classes=x$classes)
}

moveToGen <- function() {
    curx <- NA
    cury <- NA

    moveto <- function(x, dev) {
        loc <- locToInches(x$x, x$y, dev)
        curx <<- cx(loc$x, dev)
        cury <<- cy(loc$y, dev)
    }

    lineto <- function(x, dev) {
        loc <- locToInches(x$x, x$y, dev)
        lineArrow <- x$arrow
        if (! is.null(lineArrow)) {
            ends <- switch(as.character(lineArrow$ends),
                           "1" = "first",
                           "2" = "last",
                           "3" = "both")
            result <- list(x=c(curx, cx(loc$x, dev)),
                           y=c(cury, cy(loc$y, dev)),
                           arrow=list(ends = ends),
                           classes=x$classes,
                           name=x$name)
        } else {
            result <- list(x=c(curx, cx(loc$x, dev)),
                           y=c(cury, cy(loc$y, dev)),
                           classes=x$classes,
                           name=x$name)
        }
        curx <<- cx(loc$x, dev)
        cury <<- cy(loc$y, dev)
        result
    }

    list(moveto=moveto, lineto=lineto)
}

moveToFuns <- moveToGen()

devGrob.move.to <- moveToFuns$moveto
devGrob.line.to <- moveToFuns$lineto

devGrob.lines <- function(x, dev) {
  loc <- locToInches(x$x, x$y, dev)
  
  # Need to add in attributes to know where arrows
  # go if we have any
  lineArrow <- x$arrow
  if (! is.null(lineArrow)) {
      ends <- switch(as.character(lineArrow$ends),
                     "1" = "first",
                     "2" = "last",
                     "3" = "both")
      list(x=cx(loc$x, dev),
           y=cy(loc$y, dev),
           arrow=list(ends = ends),
           classes=x$classes,
           name=x$name)
  } else {
      list(x=cx(loc$x, dev),
           y=cy(loc$y, dev),
           classes=x$classes,
           name=x$name)
  }
}

devGrob.points <- function(x, dev) {
  loc <- locToInches(x$x, x$y, dev)
  list(name = x$name,
       x = cx(loc$x, dev),
       y = cy(loc$y, dev),
       size = cd(dToInches(x$size), dev),
       angle = current.rotation(),
       classes = x$classes,
       pch = x$pch)
}

devGrob.polygon <- function(x, dev) {
  loc <- locToInches(x$x, x$y, dev)
  list(x=cx(loc$x, dev),
       y=cy(loc$y, dev),
       classes=x$classes,
       name=x$name)
}

devGrob.pathgrob <- function(x, dev) {
    # The complication is converting the 'x', 'y', and 'id's
    # into lists
    if (is.null(x$id) && is.null(x$id.lengths)) {
        loc <- locToInches(x$x, x$y, dev)

        list(x=cx(loc$x, dev),
             y=cy(loc$y, dev),
             rule=x$rule,
             classes=x$classes,
             name=x$name)
    } else {
        if (is.null(x$id)) {
            n <- length(x$id.lengths)
            id <- rep(1L:n, x$id.lengths)
        } else {
            n <- length(unique(x$id))
            id <- x$id
        }
        listX <- split(x$x, id)
        listY <- split(x$y, id)
        listLoc <- mapply(locToInches, listX, listY, MoreArgs=list(dev),
                          SIMPLIFY=FALSE)
        list(x=lapply(listLoc,
               function(loc, dev) { cx(loc$x, dev) }, dev),
             y=lapply(listLoc,
               function(loc, dev) { cy(loc$y, dev) }, dev),
             rule=x$rule,
             classes=x$classes,
             name=x$name)
    }
}

devGrob.rastergrob <- function(x, dev) {
  lb <- leftbottom(x$x, x$y, x$width, x$height, x$just, x$hjust, x$vjust, dev)
  dim <- dimToInches(x$width, x$height, dev)

  list(x=cx(lb$x, dev),
       y=cy(lb$y, dev),
       width=cw(dim$w, dev),
       height=ch(dim$h, dev),
       angle=current.rotation(),
       datauri=x$datauri,
       classes=x$classes,
       name=x$name)
}

devGrob.rect <- function(x, dev) {
  lb <- leftbottom(x$x, x$y, x$width, x$height, x$just, x$hjust, x$vjust, dev)
  dim <- dimToInches(x$width, x$height, dev)
  list(x=cx(lb$x, dev),
       y=cy(lb$y, dev),
       width=cw(dim$w, dev),
       height=ch(dim$h, dev),
       angle=current.rotation(),
       classes=x$classes,
       name=x$name)
}

devGrob.text <- function(x, dev) {
  loc <- locToInches(x$x, x$y, dev)
  gp <- get.gpar()
  charHeight <- grobHeight(textGrob("M", gp = x$gp))
  # The R graphics engine does some crazy-ass calculations to
  # determine line height.  This does WAAAAY back so we just
  # have to swallow and follow along.
  # textLineHeight <-  ch(charHeight * gp$lineheight, dev)
  xcex <- if (is.null(x$gp$cex)) 1 else x$gp$cex
  textLineHeight <- ch(unit(gp$lineheight * gp$cex * xcex *
                            graphics::par("cin")[2], "inches"), dev)
  charHeight <- ch(charHeight, dev)
  
  # height of current font
  # This corresponds to lineheight in SVG terms,
  # which is defined to be font size
  # see http://www.w3.org/TR/SVG/propidx.html
  #   comment in row for 'baseline-shift' in the 'percentages' column
  # This is needed for positioning plotmath expressions
  # to anything close to the right place
  fontHeight <- ch(unit(gp$fontsize * gp$cex * xcex/ 72, "inches"), dev)

  # Width of the text/expression
  
  # MUST set x$vp to NULL before doing the following calculations
  # because x$vp has already been asserted and the calculation may
  # involve trying to assert it again!
  # (which would mean hidden error because viewport pushed twice OR
  #  visible error because try to "down" to viewport that does not exist)
  x$vp <- NULL
  
  width <- cw(grobWidth(x), dev)
  height <- ch(grobHeight(x), dev)
  ascent <- ch(grobAscent(x), dev)
  descent <- ch(grobDescent(x), dev)
  
  # Checking whether to use just or [h/v]just
  # Will convert numerics to strings in justTo_just function
  just <- rep(x$just, length.out = 2)
  just <- c(justTohjust(just[1]),
            justTovjust(just[2]))
  if (! is.null(x$hjust))
    just[1] <- justTohjust(x$hjust)
  if (! is.null(x$vjust))
    just[2] <-justTovjust(x$vjust)
  hjust <- just[1]
  vjust <- just[2]

  list(x=cx(loc$x, dev),
       y=cy(loc$y, dev),
       text=x$label,
       hjust=hjust,
       vjust=vjust,
       rot=x$rot,
       width=width,
       height=height,
       angle=current.rotation(),
       ascent=ascent,
       descent=descent,
       lineheight=textLineHeight,
       fontheight=fontHeight,
       charheight=charHeight,
       fontfamily=gp$fontfamily,
       fontface=switch(gp$font,
         "plain", "bold", "italic", "bold.italic"),
       classes=x$classes,
       name=x$name)  
}

devGrob.circle <- function(x, dev) {
  loc <- locToInches(x$x, x$y, dev)
  list(x=cx(loc$x, dev),
       y=cy(loc$y, dev),
       r=cd(dToInches(x$r), dev),
       classes=x$classes,
       name=x$name)
}

# Because viewports and grobs can be used many times, and each
# time we use one we start a group, we need a *unique* id for that
# group, otherwise things like clipping paths don't work correctly
#
# 'append' determines whether we add our ID to the usageTable. Useful
# not to in cases like animated grobs
getID <- function(name, type, append = TRUE) {
  # If this is a grob or ref, only modify if we're trying to ensure
  # uniqueness. We *really* need to do this for viewports though, so
  # viewports are a special case.
  if (type != "vp" && ! get("uniqueNames", envir = .gridSVGEnv))
      return(name)

  # Finding out how many times a VP or grob has been used so far
  ut <- get("usageTable", envir = .gridSVGEnv)
  suffix <- ut[ut$name == name, "suffix"]

  suffix <-
    if (length(suffix) == 0)
      1
    else
      max(suffix) + 1

  # Test if there are any existing names that might clash.
  # For example rect.1 has rect.1.1 children, test whether
  # these child names might clash.
  candidateName <- paste(name, suffix, sep = getSVGoption("id.sep"))
  while (length(ut[ut$name == candidateName, "suffix"])) {
    # Just increment the suffix number by 1 each time, should (eventually)
    # give us a unique number
    suffix <- suffix + 1
    candidateName <- paste(name, suffix, sep = getSVGoption("id.sep"))
  }

  if (append) {
    sel <- prefixName(escapeSelector(candidateName))
    xp <- prefixName(escapeXPath(candidateName))
    assign("usageTable",
           rbind(ut,
                 data.frame(name = name,
                            suffix = suffix,
                            type = type,
                            selector = sel,
                            xpath = xp,
                            stringsAsFactors = FALSE)),
           envir = .gridSVGEnv)
  }

  # Returning the new ID
  paste(name, suffix, sep = getSVGoption("id.sep"))
}

getCoordsInfo <- function(vp, tm, dev) {
  # Need to maintain x, y, xscale, yscale, transform
  # Units of particular interest, npc, native, inches
  # Keep inches as our baseline
  transloc <- c(0, 0, 1) %*% tm
  loc <- (transloc / transloc[3])[-3]

  coords <- list(x = round(cx(unit(loc[1], "inches"), dev), 2),
                 y = round(cy(unit(loc[2], "inches"), dev), 2),
                 width = round(cw(unit(1, "npc"), dev), 2),
                 height = round(ch(unit(1, "npc"), dev), 2),
                 angle = current.rotation(),
                 xscale = vp$xscale,
                 yscale = vp$yscale,
                 inch = round(cw(unit(1, "inches"), dev), 2))
  coords
}

devGrob.viewport <- function(x, dev) {
  vp <- x
  # Modify the path so that we can use a different separator
  if (get("use.vpPaths", envir = .gridSVGEnv)) {
    vpname <- as.character(current.vpPath())
    splitPath <- explode(vpname)
    vpname <- paste(splitPath, collapse = getSVGoption("vpPath.sep"))
  } else {
    vpname <- vp$name
  }
  coords <- getCoordsInfo(vp, current.transform(), dev)

  if (is.null(vp$clip)) {
      clip <- FALSE
      list(name=getID(vpname, "vp"), clip=clip,
           coords=coords, classes=x$classes)
  } else if (is.na(vp$clip)) {
      # Clipping has been turned OFF
      # FIXME:  CANNOT do this in SVG (enlarge the clip path)
      clip <- FALSE
      list(name=getID(vpname, "vp"), clip=clip,
           coords=coords, classes=x$classes)
  } else if (! vp$clip) {
      clip <- FALSE
      list(name=getID(vpname, "vp"), clip=clip,
           coords=coords, classes=x$classes)
  } else {
      clip <- TRUE
      list(vpx=coords$x,
           vpy=coords$y,
           vpw=coords$width,
           vph=coords$height,
           angle=current.rotation(),
           name=getID(vpname, "vp"),
           clip=clip,
           classes=x$classes,
           coords=coords)
  }
}

devGrob.vpPath <- function(x, dev) {
  vp <- current.viewport()
  tm <- current.transform()
  if (is.null(vp$clip)) {
      clip <- FALSE
      list(name=getID(vp$name, "vp"), clip=clip, classes=x$classes)
  } else if (is.na(vp$clip)) {
      # Clipping has been turned OFF
      # FIXME:  CANNOT do this in SVG (enlarge the clip path)
      clip <- FALSE
      list(name=getID(vp$name, "vp"), clip=clip, classes=x$classes)
  } else if (! vp$clip) {
      clip <- FALSE
      list(name=getID(vp$name, "vp"), clip=clip, classes=x$classes)
  } else {
      clip <- TRUE

      transloc <- c(0, 0, 1) %*% tm
      loc <- (transloc / transloc[3])[-3]
      
      list(vpx=cx(unit(loc[1], "inches"), dev),
           vpy=cy(unit(loc[2], "inches"), dev),
           vpw=cw(unit(1, "npc"), dev),
           vph=ch(unit(1, "npc"), dev),
           name=getID(vp$name, "vp"),
           classes=x$classes,
           clip=clip)
  }  
}

devGrob.clip <- function(x, dev) {
  # Should be similar to a rect in description, because this is a clipping rect
  lb <- leftbottom(x$x, x$y, x$width, x$height, x$just, x$hjust, x$vjust, dev)
  dim <- dimToInches(x$width, x$height, dev)
  list(x=cx(lb$x, dev),
       y=cy(lb$y, dev),
       width=cw(dim$w, dev),
       height=ch(dim$h, dev),
       angle=current.rotation(),
       classes=x$classes,
       name=getID(x$name, "grob"))
}

# Prim to Dev
# This generates SVG from the grob to reproduce the grob in SVG code
# General form:
#   startGroup
#   for i=1:n
#     dev&(devGrob(i))
#   endGroup
primToDev <- function(x, dev) {
  UseMethod("primToDev")
}

primToDev.grob <- function(x, dev) {
}

arrowAddName <- function(arrow, name) {
  list(angle = arrow$angle,
       length = arrow$length,
       ends = arrow$ends,
       type = arrow$type,
       name = name)
}

primToDev.clip <- function(x, dev) {
    devStartClip(devGrob(x, dev), NULL, dev)
}

primToDev.move.to <- function(x, dev) {
    devGrob(x, dev)
}

primToDev.line.to <- function(x, dev) {
  # NOTE:  MUST NOT evaluate devGrob() more than once
  #        because it has side-effects (within its closure)
  dgrob <- devGrob(x, dev)

  dgrob$name <- getID(dgrob$name, "grob")
  x$name <- getID(x$name, "grob")

  # Grouping the grob
  devStartGroup(dgrob, NULL, dev)

  # This is a bit of a special case where we know there is only one
  # actual graphical object that is being created, so we are simply
  # going to modify it's name in place.
  dgrob$name <- subGrobName(x$name, 1)

  if (! is.null(x$arrow))
    devArrow(arrowAddName(x$arrow, x$name), gparToDevPars(x$gp), dev)
  devLines(dgrob, gparToDevPars(x$gp), dev)

  # Ending the group
  devEndGroup(x$name, FALSE, dev)
}

primToDev.lines <- function(x, dev) {
  x$name <- getID(x$name, "grob")

  # Grouping the grob
  devStartGroup(devGrob(x, dev), NULL, dev)

  # This is a bit of a special case where we know there is only one
  # actual graphical object that is being created, so we are simply
  # going to modify it's name in place.
  oldname <- x$name
  x$name <- subGrobName(x$name, 1)

  if (! is.null(x$arrow))
    devArrow(arrowAddName(x$arrow, x$name), gparToDevPars(x$gp), dev)
  devLines(devGrob(x, dev), gparToDevPars(x$gp), dev)

  # Ending the group
  x$name <- oldname
  devEndGroup(x$name, FALSE, dev)
}

primToDev.polyline <- function(x, dev) {
  x$name <- getID(x$name, "grob")

  # If we only have one line
  if (is.null(x$id) && is.null(x$id.lengths)) {
      x$id <- rep(1L, length(x$x))
  }

  # Multiple lines exist
  if (is.null(x$id)) {
      n <- length(x$id.lengths)
      id <- rep(1L:n, x$id.lengths)
  } else {
      n <- length(unique(x$id))
      id <- x$id
  }

  # Each line has an id, grab corresponding positions
  listX <- split(x$x, id)
  listY <- split(x$y, id)
  n <- length(listX)

  # Gp needs to be defined for each sub-grob, as does arrow
  gp <- expandGpar(x$gp, n)
  arrows <- expandArrow(x$arrow, n)

  # Grouping each sub-grob
  devStartGroup(devGrob(x, dev), NULL, dev)

  # Now we want to create a new lineGrob for each line
  # Naming each line with the polyline name suffixed by its id
  for (i in 1:n) {
      lg <- linesGrob(x = listX[[i]],
                      y = listY[[i]],
                      gp = gp[i],
                      arrow = arrows[i],
                      default.units = x$default.units,
                      name = subGrobName(x$name, i))
      if (! is.null(lg$arrow))
          devArrow(arrowAddName(lg$arrow, lg$name), gparToDevPars(lg$gp), dev)
      devLines(devGrob(lg, dev), gparToDevPars(lg$gp), dev) 
  }

  # Ending the group
  devEndGroup(x$name, FALSE, dev)
}

# Any more efficient way of doing this?
# FIXME:  will lose any extra attributes of segments grob
primToDev.segments <- function(x, dev) {
  nx0 <- length(x$x0)
  nx1 <- length(x$x1)
  ny0 <- length(x$y0)
  ny1 <- length(x$y1)
  n <- max(nx0, nx1, ny0, ny1)

  # Gp needs to be defined for each sub-grob, as does arrow
  gp <- expandGpar(x$gp, n)
  arrows <- expandArrow(x$arrow, n)

  x$name <- getID(x$name, "grob")

  # Grouping each sub-grob
  devStartGroup(devGrob(x, dev), NULL, dev)

  for (i in 1:n) {
    lg <- linesGrob(unit.c(x$x0[(i-1) %% nx0 + 1],
                           x$x1[(i-1) %% nx1 + 1]),
                    unit.c(x$y0[(i-1) %% ny0 + 1],
                           x$y1[(i-1) %% ny1 + 1]),
                    arrow = arrows[i],
                    default.units = x$default.units,
                    gp = gp[i],
                    name = subGrobName(x$name, i))
    if (! is.null(lg$arrow))
      devArrow(arrowAddName(lg$arrow, lg$name), gparToDevPars(lg$gp), dev)
    devLines(devGrob(lg, dev), gparToDevPars(lg$gp), dev)
  }

  # Ending the group
  devEndGroup(x$name, FALSE, dev)
}

primToDev.polygon <- function(x, dev) {
  # If we have only one polygon
  if (is.null(x$id) && is.null(x$id.lengths)) {
      x$id <- rep(1L, length(x$x))
  }

  # If we have multiple polygons
  if (is.null(x$id)) {
      n <- length(x$id.lengths)
      id <- rep(1L:n, x$id.lengths)
  } else {
      n <- length(unique(x$id))
      id <- x$id
  }

  # Each polygon has an id, grab corresponding positions
  listX <- split(x$x, id)
  listY <- split(x$y, id)
  # May have id.length == 0 so use # of groups
  n <- length(listX)

  # Gp needs to be defined for each sub-grob
  gp <- expandGpar(x$gp, n)

  x$name <- getID(x$name, "grob")

  # Grouping each sub-grob
  devStartGroup(devGrob(x, dev), NULL, dev)

  # Now we want to create a new polygonGrob for each polygon
  # Naming each polygon with the polygon name suffixed by its id
  for (i in 1:n) {
      pg <- polygonGrob(x = listX[[i]],
                        y = listY[[i]],
                        gp = gp[i],
                        default.units = x$default.units,
                        name = subGrobName(x$name, i))
      devPolygon(devGrob(pg, dev), gparToDevPars(pg$gp), dev)
  }

  # Ending the group
  devEndGroup(x$name, FALSE, dev)
}

trim <- function(points) {
    n <- length(points$x)
    if (n > 2) {
        remove <- 1
        while (remove < n &&
               points$x[1] == points$x[1 + remove] &&
               points$y[1] == points$y[1 + remove]) {
            remove <- remove + 1
        }
        if (remove > 1) {
            points$x <- points$x[-(1:(remove - 1))]
            points$y <- points$y[-(1:(remove - 1))]            
        }
    }
    points
}

primToDev.xspline <- function(x, dev) {
  # Setting up function that turns an xspline into a series of points
  # which is then used to define a line or path
  splineToGrob <- function(spline) {
    splinePoints <- xsplinePoints(spline)
    if (spline$open) {
        linesGrob(x = splinePoints$x,
                  y = splinePoints$y,
                  gp = spline$gp,
                  arrow = spline$arrow,
                  default.units = spline$default.units,
                  name = spline$name)
    } else {
        pathGrob(x = splinePoints$x,
                 y = splinePoints$y,
                 gp = spline$gp,
                 default.units = spline$default.units,
                 name = spline$name)
    }
  }

  # 'grid' does not allow NAs in (x, y) for Xsplines
  if (any(is.na(x$x)) || any(is.na(x$y)))
      stop("non-finite control point in Xspline")
  
  # If we have only one spline
  if (is.null(x$id) && is.null(x$id.lengths)) {
      x$id <- rep(1L, length(x$x))
  }

  # If we're dealing with more than one spline
  if (is.null(x$id)) {
      n <- length(x$id.lengths)
      id <- rep(1L:n, x$id.lengths)
  } else {
      n <- length(unique(x$id))
      id <- x$id
  }

  # Each xspline has an id, grab corresponding positions
  listX <- split(x$x, id)
  listY <- split(x$y, id)
  n <- length(listX)

  # If x$shape is not defined for each point, repeat it for all points
  pointShapes <- rep(x$shape, length.out = length(x$x))
  listShape <- split(pointShapes, id)

  # Like x$shape, if the following attributes not defined for each grob id, repeat it
  splineOpen <- rep(x$open, length.out = n)
  splineEnds <- rep(x$repEnds, length.out = n)

  # Gp needs to be defined for each sub-grob, as does arrow
  gp <- expandGpar(x$gp, n)
  arrows <- expandArrow(x$arrow, n)

  x$name <- getID(x$name, "grob")

  # Grouping each sub-grob
  devStartGroup(devGrob(x, dev), NULL, dev)

  # Now we want to create a new xsplineGrob for each xspline
  # Naming each xspline with the xspline name suffixed by its id
  for (i in 1:n) {
      xsg <- xsplineGrob(x = listX[[i]],
                         y = listY[[i]],
                         open = x$open, # Could use splineOpen[i] but grid.xspline applies this for the entire group of grobs
                         shape = listShape[[i]],
                         default.units = x$default.units,
                         repEnds = splineEnds[i],
                         arrow = arrows[i],
                         gp = gp[i],
                         name = subGrobName(x$name, i))
      sg <- splineToGrob(xsg)
      if (inherits(sg, "pathgrob")) {
          devPath(devGrob(sg, dev), gparToDevPars(sg$gp), dev)
      } else {
          dg <- devGrob(sg, dev)
          if (! is.null(sg$arrow)) {
              devArrow(arrowAddName(sg$arrow, sg$name),
                       gparToDevPars(sg$gp), dev)
              # The arrow orientation is determined "auto"matically by
              # the SVG renderer, so we need to avoid identical values
              # at start or end of points (this has been done in
              # xsplinePoints(), but we need to do it again here because
              # we will be rounding to 2 dp for SVG output!)
              dgTrimFront <- trim(list(x=round(dg$x, 2), y=round(dg$y, 2)))
              dgTrimBack <- trim(list(x=rev(dgTrimFront$x),
                                      y=rev(dgTrimFront$y)))
              dg$x <- rev(dgTrimBack$x)
              dg$y <- rev(dgTrimBack$y)
          }
          devLines(dg, gparToDevPars(sg$gp), dev)
      }
  }

  # Ending the group
  devEndGroup(x$name, FALSE, dev)
}

primToDev.pathgrob <- function(x, dev) {
  x$name <- getID(x$name, "grob")

  # Grouping the grob
  devStartGroup(devGrob(x, dev), NULL, dev)

  # This is a bit of a special case where we know there is only one
  # actual graphical object that is being created, so we are simply
  # going to modify it's name in place.
  oldname <- x$name
  x$name <- subGrobName(x$name, 1)

  devPath(devGrob(x, dev), gparToDevPars(x$gp), dev)

  # Ending the group
  x$name <- oldname
  devEndGroup(x$name, FALSE, dev)
}

primToDev.rastergrob <- function(x, dev) {
  # Finding out how many rasters we're dealing with
  n <- max(length(x$x), length(x$y), length(x$width), length(x$height))
  # Repeating components as necessary
  xs <- rep(x$x, length.out = n)
  ys <- rep(x$y, length.out = n)

  # Finding the dimensions of the image, c(height, width)
  rasterDims <- dim(x$raster)
  rasterHeight <- rasterDims[1]
  rasterWidth <- rasterDims[2]

  # If we haven't been given any information about the h or w,
  # blow the image up to the full size but respect the aspect ratio
  x <- resolveRasterSize(x)
  
  # Use widthDetails() here (rather than grobWidth())
  # because drawing context already enforced
  widths <- rep(x$width, length.out = n)
  heights <- rep(x$height, length.out = n) 
  
  # Generating the filename of the raster
  fileloc <- tempfile(x$name, fileext = "png")

  # Because of issues regarding interpolation, it's best just to
  # store the raster with as large a dimension as possible.
  # OTOH, never want to REDUCE the size of the raw raster (?)
  rasterDims <- c(max(abs(rasterHeight), ch(max(heights), dev)),
                  max(abs(rasterWidth), cw(max(widths), dev)))

  olddev <- dev.cur()
  png(filename = fileloc, width = round(abs(rasterDims[2])),
      height = round(abs(rasterDims[1])), bg = "transparent")
      # Need to ensure that the raster is oriented correctly in the (more rare)
      # case of an xscale or yscale being big -> small
      # To do this, position natively in a new (temporary) viewport
      xscale <- if (rasterDims[2] < 0) 1:0
                else 0:1
      yscale <- if (rasterDims[1] < 0) 1:0
                else 0:1
      pushViewport(viewport(xscale = xscale, yscale = yscale),
                   recording = FALSE)
      # The raster stays the same and is only repeated for each appearance.
      # Given that we know the dimensions of the PNG, we can safely say that
      # the raster occupies the entireity of both the x and y dimensions.
      grid.raster(x$raster, width = 1, height = 1, interpolate = x$interpolate,
                  default.units = "native")
      popViewport(recording = FALSE)
  dev.off()
  dev.set(olddev)

  # base64 encoding the PNG so we can insert the image as a data URI
  base64Raster <- base64enc(fileloc)
  file.remove(fileloc)

  # Expand the gp such that it fully defines all sub-grobs
  gp <- expandGpar(x$gp, n)

  x$name <- getID(x$name, "grob")

  # Grouping each sub-grob
  devStartGroup(devGrob(x, dev), NULL, dev)

  for (i in 1:n) {
      rg <- rasterGrob(x$raster,
                       x = xs[i],
                       y = ys[i],
                       width = widths[i],
                       height = heights[i],
                       just = x$just,
                       hjust = x$hjust,
                       vjust = x$vjust,
                       default.units = x$default.units,
                       gp = gp[i], # Will be ignored, keeping anyway
                       name = subGrobName(x$name, i))
      rg$datauri <- base64Raster
      devRaster(devGrob(rg, dev), gparToDevPars(rg$gp), dev)
  }

  # Ending the group
  devEndGroup(x$name, FALSE, dev)
}

primToDev.rect <- function(x, dev) {
  # Finding out how many rects we're dealing with
  n <- max(length(x$x), length(x$y), length(x$width), length(x$height))
  # Repeating components as necessary
  xs <- rep(x$x, length.out = n)
  ys <- rep(x$y, length.out = n)
  widths <- rep(x$width, length.out = n)
  heights <- rep(x$height, length.out = n)

  # Expand the gp such that it fully defines all sub-grobs
  gp <- expandGpar(x$gp, n)

  x$name <- getID(x$name, "grob")

  # Grouping each sub-grob
  devStartGroup(devGrob(x, dev), NULL, dev)

  for (i in 1:n) {
      rg <- rectGrob(x = xs[i],
                     y = ys[i],
                     width = widths[i],
                     height = heights[i],
                     just = x$just,
                     hjust = x$hjust,
                     vjust = x$vjust,
                     default.units = x$default.units,
                     gp = gp[i],
                     name = subGrobName(x$name, i))
      devRect(devGrob(rg, dev), gparToDevPars(rg$gp), dev)
  }

  # Ending the group
  devEndGroup(x$name, FALSE, dev)
}

primToDev.text <- function(x, dev) {
  # Finding out how many pieces of text we're dealing with
  n <- max(length(x$x), length(x$y), length(x$label))
  # Repeating components as necessary
  textX <- rep(x$x, length.out = n)
  textY <- rep(x$y, length.out = n)
  textRot <- rep(x$rot, length.out = n)

  # If any given label is a vector of length 0, we don't want NA to appear
  if (length(x$label) == 0) {
    textLabel <- " "
    textLabel <- rep(textLabel, length.out = n)
  } else {
    # Checking that no element of label vector is empty
    if (!is.language(x$label)) {
        textLabel <- sapply(x$label, function(t) {
            if (nchar(t) == 0 | length(t) == 0)
                " "
            else
                t
        })
    }
    textLabel <- rep(x$label, length.out = n)
  }

  # Force fill to be col for text
  if (is.null(x$gp))
      x$gp <- gpar(fill = get.gpar()$col)
  else
      x$gp$fill <- if (! is.null(x$gp$col)) x$gp$col
                   else get.gpar()$col
  
  # Expand the gp such that it fully defines all sub-grobs
  gp <- expandGpar(x$gp, n)

  x$name <- getID(x$name, "grob")

  # Grouping each sub-grob
  devStartGroup(devGrob(x, dev), NULL, dev)

  for (i in 1:n) {
      tg <- textGrob(x = textX[i],
                     y = textY[i],
                     label = textLabel[i],
                     rot = textRot[i],
                     just = x$just,
                     hjust = x$hjust,
                     vjust = x$vjust,
                     default.units = x$default.units,
                     gp = gp[i],
                     name = subGrobName(x$name, i))
      devText(devGrob(tg, dev), gparToDevPars(tg$gp), dev)
  }

  # Ending the group
  devEndGroup(x$name, FALSE, dev)
}

primToDev.circle <- function(x, dev) {
  # Finding out how many circles we're dealing with
  n <- max(length(x$x), length(x$y), length(x$r))
  # Repeating components as necessary
  xs <- rep(x$x, length.out = n)
  ys <- rep(x$y, length.out = n)
  rs <- rep(x$r, length.out = n)

  # Expand the gp such that it fully defines all sub-grobs
  gp <- expandGpar(x$gp, n)

  x$name <- getID(x$name, "grob")

  # Grouping each sub-grob
  devStartGroup(devGrob(x, dev), NULL, dev)

  for (i in 1:n) {
      cg <- circleGrob(x = xs[i],
                       y = ys[i],
                       r = rs[i],
                       default.units = x$default.units,
                       gp = gp[i],
                       name = subGrobName(x$name, i))
      devCircle(devGrob(cg, dev), gparToDevPars(cg$gp), dev)
  }

  # Ending the group
  devEndGroup(x$name, FALSE, dev)
}

adjustSymbolSize <- function(pointSize, pgp) {
    # Points are affected by cex and fontsize but only if they are
    # char or lines, etc
    # Solution: push a viewport with new gps from the grob and can
    # therefore can convert unit safely to inches because grid's unit
    # conversion routines can handle when the *viewport* has the cex or
    # fontsize information but not when the *grob* has it.
    # Also, not recording on the DL because this viewport wasn't part
    # of the original vp tree.
    if (! is.null(pgp$cex) || ! is.null(pgp$fontsize)) {
        xscale <- current.viewport()$xscale
        yscale <- current.viewport()$yscale
        if (! (is.null(pgp$cex) & is.null(pgp$fontsize))) {
            pushViewport(viewport(xscale = xscale, yscale = yscale,
                                  gp = gpar(cex = pgp$cex,
                                      fontsize = pgp$fontsize)),
                         recording = FALSE)
        } else if (! is.null(pgp$cex)) {
            pushViewport(viewport(xscale = xscale, yscale = yscale,
                                  gp = gpar(cex = pgp$cex)),
                         recording = FALSE)
        } else {
            # if (! is.null(pgp$fontsize))
            pushViewport(viewport(xscale = xscale, yscale = yscale,
                                  gp = gpar(fontsize = pgp$fontsize)),
                         recording = FALSE)
        }
        pointSize <- convertWidth(pointSize, "inches") # Use width, matches grid
        popViewport(recording = FALSE)
    }
    pointSize
}

primToDev.points <- function(x, dev) {
    # Finding out how many grobs we're going to be dealing with
    # length of x and y already checked in grid.points
    n <- length(x$x)

    # Expand the gp such that it fully defines all sub-grobs
    gp <- expandGpar(x$gp, n)

    x$name <- getID(x$name, "grob")

    # Grouping each sub-grob
    devStartGroup(devGrob(x, dev), NULL, dev) 

    # For testing validity, convert to numerics
    chinds <- which(! as.character(x$pch) %in% as.character(c(0:25, 32:127)))
    pchtest <- x$pch
    if (length(chinds) > 0) {
        newpch <- integer(length(pchtest))
        newpch[chinds] <- as.numeric(sapply(pchtest[chinds],
                                            function(x) charToRaw(x)))
        newpch[!chinds] <- as.numeric(pchtest[!chinds])
        pchtest <- newpch
    }

    if (any(!pchtest %in% c(0:25, 32:127)))
        stop("Unsupported pch value")

    # These can differ for points
    pchs <- rep(pchtest, length.out = n)
    sizes <- rep(x$size, length.out = n)

    # Assume we need to define the point symbol
    createDef <- TRUE

    for (i in 1:n) {
        # Check whether the point symbol has been used yet
        pchUsageTable <- get("pchUsageTable", envir = .gridSVGEnv)
        # Update usages
        pchUsageTable[pchs[i] + 1, "used"] <- TRUE
        assign("pchUsageTable", pchUsageTable, envir = .gridSVGEnv)

        pgp <- gp[i]

        if (! is.unit(sizes[i]) && is.numeric(sizes[i])) {
            # Just a number -- convert to a unit
            pointSize <- unit(sizes[i], x$default.units)
        } else {
            # All other units
            pointSize <- sizes[i]
        }
        
        asciipch <- if (pchs[i] %in% 32:127)
                        rawToChar(as.raw(pchs[i]))
                    else
                        pchs[i]

        # Force a stroke-width
        pgp$lwd <- if (is.null(pgp$lwd)) get.gpar()$lwd
                   else pgp$lwd

        if (pchs[i] < 15)
            pgp$fill <- "transparent"

        if (pchs[i] %in% 15:20 | pchs[i] >= 32) {
            # 46 == "."
            # Don't do anything for a "." because we need a
            # stroke for it to be visible
            if (pchs[i] != 46) {
                pgp$fill <- if (is.null(pgp$col)) get.gpar()$col
                            else pgp$col
                if (pchs[i] %in% 15:18 | pchs[i] >= 32)
                    pgp$col <- "transparent"
            }
        }

        # Size is now relative to text so use text grob
        if (pchs[i] >= 32)
            pointSize <- grobWidth(textGrob(asciipch))

        # Enforce gp$cex or gp$fontsize
        pointSize <- adjustSymbolSize(pointSize, pgp)
        
        devUseSymbol(devGrob(pointsGrob(x$x[i], x$y[i],
                                        pch = asciipch,
                                        size = pointSize,
                                        default.units = x$default.units,
                                        name = subGrobName(x$name, i)), dev),
                     gparToDevPars(pgp), dev)
    }

    # Ending the group
    devEndGroup(x$name, FALSE, dev)
}

grobToDev.gTree <- function(x, dev) {
    depth <- enforceVP(x$vp, dev)
    if (!is.null(x$childrenvp)) {
        pushViewport(x$childrenvp, recording=FALSE)
        upViewport(depth(x$childrenvp), recording=FALSE)
    }
    primToDev(x, dev)
    unwindVP(x$vp, depth, dev)
    # Ignore wrapping gTree as it was not on the original DL
    if (x$name != "gridSVG")
        progressStep("grob")
}

primToDev.gTree <- function(x, dev) {
    if (x$name != "gridSVG") {
        x$name <- getID(x$name, "grob")
        x$classes <- class(x)
        children <- x$children[x$childrenOrder]
    } else {
        children <- x$children
    }
    devStartGroup(devGrob(x, dev), gparToDevPars(x$gp), dev)
    lapply(children, function(child) {
        # 'gridSVG' is a special case because it is just a wrapping gTree.
        # It is not useful for us to track the entire gPath as a result,
        # only the path *after* 'gridSVG'
        if (get("use.gPaths", envir = .gridSVGEnv) && x$name != "gridSVG")
            child$name <- paste(x$name, child$name, sep = getSVGoption("gPath.sep"))
        child$classes <- class(child)
        grobToDev(child, dev)
    })
    devEndGroup(x$name, FALSE, dev)
}

# Viewports (and vpPaths and downs and ups)
# on the display list get recorded as wrapped grobs
grobToDev.recordedGrob <- function(x, dev) {
    x <- x$list
    if (!is.null(x$vp)) { # recorded pushViewport
        enforceVP(x$vp, dev)
    } else if (!is.null(x$path)) { # recorded downViewport
        enforceVP(x$path, dev)
    } else if (!is.null(x$n)) { # recorded up or pop
        unwindVP(NULL, x$n, dev)
    }
}

# grid to SVG
# Given a gTree created by grid.grab()
gridToDev <- function(gTree, dev) {
  grobToDev(gTree, dev)
}


