xmlDecl <- function() {
  paste0('<?xml version="1.0" encoding="', localeToCharset()[1], '"?>\n')
}

htmlFile <- function(filename, svgdev) {
  # For viewing using Adobe SVG Viewer in IE
  # OR in Firefox 3 (native support)
  # create a "wrapper" html file
    # NOTE that for including plotmath output (as MathML), may
    # need to use the right sort of headers.
    # See ~/Research/Rstuff/SVG/PlotMath/README for notes from some
    # experiments AND email from David Scott that contains an
    # example from org-babel output 2011-11-01
  htmlfile <- paste0(filename, ".html")
    # NOTE that different browsers prefer different approaches
    # See email from David Scott 2011-11-03 for some sample code
  # The empty text node is so that we ensure the object tag is not
  # self-closing, i.e. there is an explicit closing tag written out
  obj <- newXMLNode("object",
                    attrs = list(data = filename,
                                 type = "image/svg+xml",
                                 width = paste0(ceiling(svgDevWidth(svgdev)), "px"),
                                 height = paste0(ceiling(svgDevHeight(svgdev)), "px")),
                    newXMLTextNode(""))
  fn <- saveXML(obj, file = htmlfile)
}

svgOpen <- function(width=200, height=200, strict=TRUE, rootAttrs=NULL) {
  # Ensure all vp contexts are now zero
  assign("contextLevels", 0, envir = .gridSVGEnv)
  svgdev <- svgDevice(width, height, strict, rootAttrs)
  svgHeader(width, height, svgdev)
  return(svgdev)
}

svgClose <- function(svgdev) {
  # Ensure all vp contexts are now zero
  assign("contextLevels", 0, envir = .gridSVGEnv)
  return(xmlRoot(svgDevParent(svgdev)))
}

svgJSUtils <- function(exportJS, svgfile, svgroot) {
  utilsFn <- paste0(svgfile, ".utils.js")
  utilsFile <- file(system.file("js/utils.js", package = "gridSVG"))
  utilsLines <- readLines(utilsFile)
  close(utilsFile)
  if (exportJS == "file") {
    destFile <- file(utilsFn)
    writeLines(utilsLines, destFile)
    close(destFile)
    newXMLNode("script", parent = svgroot, at = 0,
                      attrs = list(type = "application/ecmascript",
                                   "xlink:href" = utilsFn))
  }

  if (exportJS == "inline") {
    newXMLNode("script", parent = svgroot, at = 0,
               attrs = list(type = "application/ecmascript"),
               newXMLCDataNode(paste0(c("", utilsLines, ""), collapse = "\n")))
  }

  # When we don't want to write to a file we might want to retain some
  # info, thus just return the JS quietly
  invisible(paste(utilsLines, collapse = "\n"))
}

svgCoords <- function(exportCoords, svgfile, svgroot) {
  coordsJSON <- toJSON(get("vpCoords", envir = .gridSVGEnv))
  coordsJSON <- paste("var gridSVGCoords = ", coordsJSON, ";", sep = "")

  if (exportCoords == "file") {
    coordsFn <- paste0(svgfile, ".coords.js")
    coordsFile <- file(coordsFn, "w")
    cat(coordsJSON, "\n", file = coordsFile, sep = "")
    close(coordsFile)
    newXMLNode("script", parent = svgroot, at = 0,
               attrs = list(type = "application/ecmascript",
                            "xlink:href" = coordsFn))
  }

  if (exportCoords == "inline") {
    newXMLNode("script", parent = svgroot, at = 0,
               attrs = list(type = "application/ecmascript"),
               newXMLCDataNode(paste0(c("", coordsJSON, ""), collapse = "\n")))
  }

  # When we don't want to write to a file we might want to retain some
  # info, thus return coords info quietly
  invisible(get("vpCoords", envir = .gridSVGEnv))
}

svgMappings <- function(exportMappings, svgfile, svgroot) {
  usageTable <- get("usageTable", envir = .gridSVGEnv)

  if (exportMappings == "file") {
    mappingsFn <- paste0(svgfile, ".mappings.js")
    mappingsFile <- file(mappingsFn, "w")
    cat(exportMappings(usageTable), file = mappingsFile)
    close(mappingsFile)
    newXMLNode("script", parent = svgroot, at = 0,
               attrs = list(type = "application/ecmascript",
                            "xlink:href" = mappingsFn))
  }

  if (exportMappings == "inline") {
    newXMLNode("script", parent = svgroot, at = 0,
               attrs = list(type = "application/ecmascript"),
               newXMLCDataNode(exportMappings(usageTable)))
  }

  # When we don't want to write to a file we might want to retain some
  # info, thus return coords info quietly
  invisible(formatMappings(usageTable))
}

svgAnnotate <- function(svgRoot, callAttrs) {
    # The purpose of this function is to collate all the information
    # that gridSVG knows about as it being called. Provides us with a
    # method of potentially debugging output and version detection.
    #
    # We put all this information in a comment node so that the output
    # is not parsed by a viewer.
    # However, if we are able to get the *text* from the comment we want
    # to be able to *parse* the output.
    
    argNames <- names(callAttrs)
    argValues <- unname(unlist(callAttrs))

    # The call elements that we're going to be building up
    metadata <- newXMLNode("metadata", namespaceDefinitions =
                           c(gridsvg = "http://www.stat.auckland.ac.nz/~paul/R/gridSVG/"))

    # Using the package DESCRIPTION version instead of packageVersion
    # because packageVersion converts our versions from 1.0-0 to 1.0.0.
    # Ignoring timezone in Sys.time(), should be fine
    newXMLNode("generator",
               namespace = "gridsvg",
               attrs = c(name = "gridSVG",
                         version = packageDescription("gridSVG")$Version,
                         time = as.character(Sys.time())),
               parent = metadata)

    for (i in 1:length(callAttrs)) {
        newXMLNode("argument",
                   namespace = "gridsvg",
                   attrs = c(name = argNames[i], value = argValues[i]),
                   parent = metadata)
    }

    seps <- unlist(getSVGoptions())
    for (i in 1:length(seps)) {
        newXMLNode("separator",
                   namespace = "gridsvg",
                   attrs = c(name = names(seps[i]),
                             value = unname(seps[i])),
                   parent = metadata)
    }

    # at = 0 because we want this comment to be inserted directly after
    # the main <svg> element
    addChildren(svgRoot, metadata, at = 0)
}

svgComment <- function(comment, svgdev=svgDevice()) {
  # If this is a multi-line comment, to ensure comments have the same
  # indentation, prefix and suffix the comment with empty lines
  if (length(comment) > 1)
    comment <- paste0(c("", comment, ""), collapse="\n")
  newXMLCommentNode(comment, parent = svgDevParent(svgdev))
}

# <clipPath>, <rect>, <raster>, and <text> elements MAY
# have a rotation angle
svgAngleTransform <- function(x, y, angle) {
    if (!is.null(angle) && angle != 0) {
        paste0("rotate(", round(angle, 2), " ",
               round(x, 2), " ", round(y, 2), ")")
    } else {
        NULL
    }
}

svgClipPath <- function(id, vpx, vpy, vpw, vph, vpa,
                        svgdev=svgDevice()) {
  clipPathID <- prefixName(paste(id, "clipPath",
                                 sep = getSVGoption("id.sep")))

  # Correct w/h if necessary
  if (vpw < 0) {
    vpx <- vpx + vpw # shifts x to the left
    vpw <- abs(vpw)
  }

  if (vph < 0) {
    vpy <- vpy + vph # shifts y down
    vph <- abs(vph)
  }
  
  newXMLNode("defs", parent = svgDevParent(svgdev),
             newXMLNode("clipPath",
                        attrs = attrList(list(id = clipPathID,
                            transform=svgAngleTransform(vpx, vpy, vpa))),
                        newXMLNode("rect",
                                   attrs = list(x = round(vpx, 2),
                                                y = round(vpy, 2),
                                                width = round(vpw, 2),
                                                height = round(vph, 2),
                                                fill = "none",
                                                stroke = "none"))))
}

svgClipAttr <- function(id, clip) {
  if (clip)
    list("clip-path" = paste0("url(#", prefixName(id),
                              getSVGoption("id.sep"), "clipPath)"))
  else
    list()
}

svgMaskAttr <- function(id, mask) {
  if (mask)
    list("mask" = paste0("url(#", prefixName(id),
                         getSVGoption("id.sep"), "mask)"))
  else
    list()
}

svgStartElement <- function(id = NULL, classes = NULL, element = NULL, attrs = NULL,
                            namespace = NULL, namespaceDefinitions = NULL,
                            attributes=svgAttrib(), links=NULL, show = NULL,
                            svgdev = svgDevice()) {
  has.link <- hasLink(links[id])
  if (has.link)
    svgStartLink(links[id], show[id], svgdev)

  attrs$id <- prefixName(id)
  # If garnishing, clobber any existing attrs
  for (name in names(attributes))
      attrs[[name]] <- attributes[[name]]

  # Avoid clobbering "class" attribute if it exists
  # Instead, add to the list of classes available
  if (! is.null(attrs$class) && get("addClasses", envir = .gridSVGEnv)) {
      cls <- strsplit(attrs$class, "\\s")[[1]]
      cls <- cls[nzchar(cls)] # Get rid of whitespace
      classList <- svgClassList(unique(c(cls, classes)))
      attrs$class <- classList$class
  } else {
      classList <- svgClassList(classes)
      if (length(classList))
          attrs$class <- classList$class
  }

  attrs <- attrList(attrs)
  element <- newXMLNode(element, attrs = attrs,
                        namespace =
                            if (is.null(namespace))
                                character()
                            else 
                                namespace,
                        namespaceDefinitions =
                            if (is.null(namespaceDefinitions))
                                character()
                            else
                                namespaceDefinitions,
                        parent = svgDevParent(svgdev))
  svgDevChangeParent(element, svgdev)
}

# This is pretty much the same as svgEndGroup
svgEndElement <- function(id=NULL, links=NULL, svgdev=svgDevice()) {
  # In the case where we've got a link on our element, set the parent
  # one level up because we've got an "a" tag above the group
  has.link <- hasLink(links[id])
  if (has.link)
    svgEndLink(svgdev)

  svgDevChangeParent(xmlParent(svgDevParent(svgdev)), svgdev)
}

svgTextNode <- function(text, svgdev = svgDevice()) {
  newXMLTextNode(text, parent = svgDevParent(svgdev))
}

svgStartGroup <- function(id=NULL, clip=FALSE, mask=FALSE,
                          attributes=svgAttrib(), links=NULL, show=NULL,
                          style=svgStyle(), coords=NULL, classes = NULL,
                          svgdev=svgDevice()) {
  # If this is a viewport that we're starting a group for
  # we will have coordinate information, otherwise don't bother.
  if (! is.null(coords)) {
    currVpCoords <- get("vpCoords", envir = .gridSVGEnv)
    currId <- prefixName(getid(id, svgdev))
    currVpCoords[[currId]] <- coords
    assign("vpCoords", currVpCoords, envir = .gridSVGEnv)
  }

  has.link <- hasLink(links[id])
  if (has.link)
    svgStartLink(links[id], show[id], svgdev)

  attrlist <- c(list(id = prefixName(id)),
                svgClipAttr(id, clip),
                svgMaskAttr(id, mask),
                svgStyleAttributes(style, svgdev),
                svgAttribTxt(attributes, id, "g", svgdev))
  attrlist <- attrList(attrlist)

  # Avoid clobbering "class" attribute if it exists
  # Instead, add to the list of classes available
  if (! is.null(attrlist$class) && get("addClasses", envir = .gridSVGEnv)) {
      cls <- strsplit(attrlist$class, "\\s")[[1]]
      cls <- cls[nzchar(cls)] # Get rid of whitespace
      classList <- svgClassList(unique(c(cls, classes)))
      attrlist$class <- classList$class
  } else {
      classList <- svgClassList(classes)
      if (length(classList))
          attrlist$class <- classList$class
  }

  newparent <- newXMLNode("g", parent = svgDevParent(svgdev),
                          attrs = attrlist)
  svgDevChangeParent(newparent, svgdev)
}

svgEndGroup <- function(id=NULL, links=NULL, vp=FALSE, svgdev=svgDevice()) {
  # Handle case where clipGrobs, clipPath grobs and maskGrobs
  # have started groups. "pop" until we reach the appropriate group
  if (vp) {
    # In the case where we have reached something we know
    # is a viewport, then we don't need to unwind further
    parentIsVP <- function() {
      id <- xmlGetAttr(svgDevParent(svgdev), "id")
      ut <- get("usageTable", envir = .gridSVGEnv)
      ut <- ut[ut$type == "vp", ]
      baseGrobName(id) %in% ut$name
    }
    contextLevel <- tail(get("contextLevels", envir = .gridSVGEnv), 1)
    while (! parentIsVP() && contextLevel > 0) {
      svgDevChangeParent(xmlParent(svgDevParent(svgdev)), svgdev)
      contextLevel <- contextLevel - 1
    }
    # Remove latest vp from list of contexts
    assign("contextLevels",
           head(get("contextLevels", envir = .gridSVGEnv), -1),
           envir = .gridSVGEnv)
  } else {
    # In the case where we've got a link on our group, set the parent
    # one level up because we've got an "a" tag above the group.
    # Only doing this in the case where we're dealing with a grob.
    has.link <- hasLink(links[id])
    if (has.link)
      svgEndLink(svgdev)
  }

  svgDevChangeParent(xmlParent(svgDevParent(svgdev)), svgdev)
}

svgStartSymbol <- function(pch, svgdev = svgDevice()) {
  symbol <- newXMLNode("symbol", parent = svgDevParent(svgdev), at = 0,
                       attrs = list(id = prefixName(paste0("gridSVG.pch", pch)),
                                    viewBox = "-5 -5 10 10",
                                    overflow = "visible"))
  svgDevChangeParent(symbol, svgdev)
}

svgEndSymbol <- function(svgdev = svgDevice()) {
  # Close symbol
  svgDevChangeParent(xmlParent(svgDevParent(svgdev)), svgdev)
}

svgStartLink <- function(href="", show="", svgdev=svgDevice()) {
  linkAttrs <- list("xlink:href" = href)
  if (! is.null(show) && ! is.na(show) && nchar(show))
      linkAttrs$`xlink:show` <- show

  link <- newXMLNode("a",
                     parent = svgDevParent(svgdev),
                     attrs = linkAttrs) 
  svgDevChangeParent(link, svgdev)
}

svgEndLink <- function(svgdev=svgDevice()) {
  parent <- xmlParent(svgDevParent(svgdev))
  svgDevChangeParent(parent, svgdev)
}

svgAnimate <- function(attrib, values,
                       begin, interp, duration, rep, revert, id=NULL,
                       svgdev=svgDevice()) {
  n <- if (is.null(id)) 1 else length(unique(id))

  newXMLNode("animate", parent = svgDevParent(svgdev),
             attrs = list("xlink:href" = paste0("#", prefixName(getid(id, svgdev, n))),
                          attributeName = attrib,
                          begin = paste0(begin, "s"),
                          calcMode = interp,
                          dur = paste0(duration, "s"),
                          values = values,
                          repeatCount = if (is.numeric(rep)) rep else if (rep) "indefinite" else 1,
                          fill = if (revert) "remove" else "freeze"))
}

# Special case just for stroke-width
# values here is a vector of *numeric* values, not just
# a single element character vector (e.g. 'svgAnimate')
svgAnimatePointSW <- function(values,
                              begin, interp, duration, rep, revert,
                              id=NULL, svgdev=svgDevice()) {
  n <- if (is.null(id)) 1 else length(unique(id))
  keyTimes <- round(seq(from = 0, to = 1, length.out = length(values)), 2)
  # Change the spline depending on whether we're increasing
  # the "size" of the stroke width or decreasing
  keySplines <- -diff(values)
  keySplines <- sapply(keySplines, function(x) {
                         if (x >= 0)
                           "0 1" # point is growing
                         else
                           "1 0" # point is shrinking
                       })
  keySplines <- paste(keySplines, "1 1", collapse = ";")
  keyTimes <- paste0(round(keyTimes, 2), collapse = ";")
  values <- paste0(round(values, 2), collapse = ";")

  newXMLNode("animate", parent = svgDevParent(svgdev),
             attrs = list("xlink:href" = paste0("#", prefixName(getid(id, svgdev, n))),
                          attributeName = "stroke-width",
                          begin = paste0(begin, "s"),
                          calcMode = "spline",
                          dur = paste0(duration, "s"),
                          values = values,
                          repeatCount = if (is.numeric(rep)) rep else if (rep) "indefinite" else 1,
                          fill = if (revert) "remove" else "freeze",
                          keyTimes = keyTimes,
                          keySplines = keySplines))
}

# This and svgAnimateY are untested with id != NULL
# and I have a strong suspicion there may be problems
# because tapply returns a list -- see svgAnimatePoints
# for ideas for a possible solution (esp. the lpaste function)
svgAnimateXYWH <- function(attrib, values,
                           begin, interp, duration, rep, revert,
                           id=NULL,
                           svgdev=svgDevice()) {
  svgAnimate(attrib,
             paste(round(values, 2), collapse=";"),
             begin, interp, duration, rep, revert, id, svgdev)
}

# DON'T call this with a list of length < 2!
old.lpaste <- function(alist, collapse) {
  n <- length(alist)
  if (n == 2)
    result <- paste(alist[[1]], alist[[2]])
  else
    result <- paste(alist[[n]], lpaste(alist[1:(n-1)], collapse))
  paste(result, collapse=collapse)
}

lpaste <- function(alist, collapse) {
  n <- length(alist)
  result <- alist[[1]]
  for (i in 2:n)
    result <- paste(result, alist[[i]])
  paste(result, collapse=collapse)
}

svgAnimatePoints <- function(xvalues, yvalues, timeid,
                             begin, interp, duration, rep, revert,
                             id=NULL,
                             svgdev=svgDevice()) {
  if (is.null(id))
    warning("Only one point to animate")
  else
    svgAnimate("points",
                paste(lapply(split(paste(round(xvalues, 2),
                                         round(yvalues, 2), sep=","),
                                   timeid),
                             paste, collapse=" "),
                      collapse=";"),
               begin, interp, duration, rep, revert, id, svgdev)
}

svgAnimatePath <- function(xvalues, yvalues, pathid, timeid,
                           begin, interp, duration, rep, revert,
                           id=NULL,
                           svgdev=svgDevice()) {
  if (is.null(id))
    warning("Not sure what this animation means?")
  else {
      # Split into time segments
      x <- split(xvalues, timeid)
      y <- split(yvalues, timeid)
      pid <- split(pathid, timeid)
      d <- mapply(function(xtime, ytime, pid) {
                      # Split into path components
                      xx <- split(xtime, pid)
                      yy <- split(ytime, pid)
                      txt <- mapply(function(x, y) {
                                        paste(paste(c("M",
                                                      rep("L", length(x) - 1)),
                                                    round(x, 2), round(y, 2),
                                                    collapse=" "),
                                              "Z")
                                    }, xx, yy)
                      paste(unlist(txt), collapse=" ")
                  }, x, y, pid)
      svgAnimate("d", paste(d, collapse=";"),
                 begin, interp, duration, rep, revert, id, svgdev)
  }
}

svgAnimateTransform <- function(attrib, values,
                                begin, interp, duration, rep, revert,
                                additive = "replace",
                                id=NULL,
                                svgdev=svgDevice()) {
  n <- if (is.null(id)) 1 else length(unique(id))
  newXMLNode("animateTransform", parent = svgDevParent(svgdev),
             attrs = list("xlink:href" = paste0("#", prefixName(getid(id, svgdev, n))),
                          attributeName = "transform",
                          type = attrib,
                          begin = paste0(begin, "s"),
                          calcMode = interp,
                          dur = paste0(duration, "s"),
                          values = values,
                          additive = additive,
                          repeatCount = if (is.numeric(rep)) rep else if (rep) "indefinite" else 1,
                          fill = if (revert) "remove" else "freeze"))
}

svgAnimateTranslation <- function(xvalues, yvalues,
                                  begin, interp, duration, rep, revert,
                                  additive = "replace",
                                  id=NULL,
                                  svgdev=svgDevice()) {
  svgAnimateTransform("translate",
                      paste(round(xvalues, 2),
                            round(yvalues, 2),
                            sep=",", collapse=';'),
                      begin, interp, duration, rep, revert,
                      additive, id, svgdev)
}

svgAnimateRotation <- function(angle, xvalues, yvalues,
                               begin, interp, duration, rep, revert,
                               additive = "replace",
                               id=NULL,
                               svgdev=svgDevice()) {
  svgAnimateTransform("rotate",
                      paste(round(angle, 2),
                            round(xvalues, 2),
                            round(yvalues, 2),
                            sep=" ", collapse=';'),
                      begin, interp, duration, rep, revert,
                      additive, id, svgdev)
}

svgAnimateScale <- function(xvalues, yvalues,
                            begin, interp, duration, rep, revert,
                            additive = "replace",
                            id=NULL,
                            svgdev=svgDevice()) {
  svgAnimateTransform("scale",
                      paste(round(xvalues, 2),
                            round(yvalues, 2),
                            sep=",", collapse=';'),
                      begin, interp, duration, rep, revert,
                      additive, id, svgdev)
}

svgLines <- function(x, y, id=NULL, arrow = NULL,
                     attributes=svgAttrib(), links=NULL, show=NULL,
                     style=svgStyle(), svgdev=svgDevice()) {

  # Never fill a line
  style$fill <- "none"

  has.link <- hasLink(links[id])
  if (has.link)
    svgStartLink(links[id], show[id], svgdev)

  # Handle NA values in (x, y)
  xylist <- splitOnNA(x, y)
  N <- length(xylist)
  if (N > 0) {
      if (N > 1) {
          # assume no more than 26 NAs in (x, y)!
          alpha <- letters[1:N]
      } else {
          alpha <- ""
      }
      for (i in 1:N) {
          # Grabbing arrow info for marker element references
          # Arrows NOT drawn at NA splits
          if (! is.null(arrow$ends)) {
              if (arrow$ends == "both" && i == 1 && N == 1) 
                  lineMarkerTxt <- markerTxt("both", id)
              else if ((arrow$ends == "both" || arrow$ends == "first") &&
                       i == 1 && is.finite(x[1]) && is.finite(y[1]))
                  lineMarkerTxt <- markerTxt("first", id)
              else if ((arrow$ends == "both" || arrow$ends == "last") &&
                       i == N &&
                       is.finite(x[length(x)]) && is.finite(y[length(y)]))
                  lineMarkerTxt <- markerTxt("last", id)
              else 
                  lineMarkerTxt <- NULL
          } else {
              lineMarkerTxt <- NULL
          }
          attrlist <- c(list(id = prefixName(paste0(id, alpha[i])),
                             points = paste0(round(xylist[[i]]$x, 2), ",",
                                 round(xylist[[i]]$y, 2),
                                 collapse=" ")),
                        lineMarkerTxt,
                        svgStyleAttributes(style, svgdev),
                        svgAttribTxt(attributes, id, "polyline", svgdev))
          attrlist <- attrList(attrlist)
          newXMLNode("polyline", parent = svgDevParent(svgdev),
                     attrs = attrlist)
      }
  }
  
  if (has.link)
    svgEndLink(svgdev)
}

svgMarker <- function(x, y, type, ends, direction, name,
                      style=svgStyle(), svgdev=svgDevice()) {
    width <- abs(max(x) - min(x))
    height <- abs(max(y) - min(y))
    if (length(x) != length(y))
        stop("x and y must be same length")
    if (is.atomic(x)) {
        if (is.atomic(y)) {
            x <- list(x)
            y <- list(y)
        } else {
            stop("'x' and 'y' must both be lists or both be atomic")
        }
    }

    d <- mapply(
                function(subx, suby) {
                    openPath <- paste(c("M",
                                      rep("L", length(subx) - 1)),
                                      round(subx, 2), round(suby, 2),
                                      collapse=" ")
                    if (type == 2) # Closed arrow
                      paste(openPath, "Z")
                    else
                      openPath
                }, x, y)

    # If the arrow is open, we don't want to fill it
    if (type == 1)
        style$fill <- "none"

    # [[1]] and [1]: markerStart
    # [[2]] and [2]: markerEnd
    # pathattrs is simply a list where each element
    # is a list that we can simply pass in as attrs
    # to newXMLNode
    ids <- markerName("both", name)
    refXs <- direction * round(c(-width, width), 2)
    refYs <- round(c(-height / 2, height / 2), 2)
    pathlist <- attrList(c(list(d = d), svgStyleAttributes(style, svgdev)))
    # It is possible for width to be 0, i.e. when angle=90.
    # Ensure that the marker is always at least as wide as the
    # stroke width that it is given.
    mwidth <- max(as.numeric(pathlist$`stroke-width`), width)
    mheight <- max(as.numeric(pathlist$`stroke-width`), height)
    pathattrs <- list(pathlist, pathlist)
    pathattrs[[1]]$transform <- "rotate(180)"

    newXMLNode("defs", parent = svgDevParent(svgdev),
               newXMLNode("marker",
                          attrs = list(id = ids[1],
                                       refX = refXs[1],
                                       refY = refYs[1],
                                       overflow = "visible",
                                       markerUnits = "userSpaceOnUse",
                                       markerWidth = round(mwidth, 2),
                                       markerHeight = round(mheight, 2),
                                       orient = "auto"),
                          newXMLNode("path", attrs = pathattrs[[1]])),
               newXMLNode("marker",
                          attrs = list(id = ids[2],
                                       refX = refXs[2],
                                       refY = refYs[2],
                                       overflow = "visible",
                                       markerUnits = "userSpaceOnUse",
                                       markerWidth = round(mwidth, 2),
                                       markerHeight = round(mheight, 2),
                                       orient = "auto"),
                          newXMLNode("path", attrs = pathattrs[[2]])))
}

markerTxt <- function(ends, name) {
    mname <- markerName(ends, name)

    if (ends == "first")
        lmt <- list("marker-start" = paste0("url(#", mname, ")"))
    if (ends == "last")
        lmt <- list("marker-end" = paste0("url(#", mname, ")"))
    if (ends == "both")
        lmt <- list("marker-start" = paste0("url(#", mname[1], ")"),
                    "marker-end" = paste0("url(#", mname[2], ")"))
    lmt
}

markerName <- function(ends, name) {
    if (ends == "first")
        mname <- paste(name, getSVGoption("id.sep"), "markerStart", sep="")
    if (ends == "last")
        mname <- paste(name, getSVGoption("id.sep"), "markerEnd", sep="")
    if (ends == "both")
        mname <- c(paste(name, getSVGoption("id.sep"), "markerStart", sep=""),
                   paste(name, getSVGoption("id.sep"), "markerEnd", sep=""))
    prefixName(mname)
}

svgPolygon <- function(x, y, id=NULL,
                       attributes=svgAttrib(), links=NULL, show=NULL,
                       style=svgStyle(), svgdev=svgDevice()) {
  if (length(x) != length(y))
    stop("x and y must be same length")

  has.link <- hasLink(links[id])
  if (has.link)
    svgStartLink(links[id], show[id], svgdev)

  # Handle NA values in (x, y)
  xylist <- splitOnNA(x, y)
  N <- length(xylist)
  if (N > 0) {
      if (N > 1) {
          # assume no more than 26 NAs in (x, y)!
          alpha <- letters[1:N]
      } else {
          alpha <- ""
      }
      for (i in 1:N) {
          attrlist <- c(list(id = prefixName(paste0(id, alpha[i])),
                             points = paste0(round(xylist[[i]]$x, 2), ",",
                                 round(xylist[[i]]$y, 2),
                                 collapse=" ")),
                        svgStyleAttributes(style, svgdev),
                        svgAttribTxt(attributes, id, "polygon", svgdev))
          attrlist <- attrList(attrlist)
          newXMLNode("polygon", parent = svgDevParent(svgdev),
                     attrs = attrlist)
      }
  }

  if (has.link)
    svgEndLink(svgdev)
}

# Differs from polygon because it can have sub-paths
svgPath <- function(x, y, rule, id=NULL,
                    attributes=svgAttrib(), links=NULL, show=NULL,
                    style=svgStyle(), svgdev=svgDevice()) {
    if (length(x) != length(y))
        stop("x and y must be same length")
    if (is.atomic(x)) {
        if (is.atomic(y)) {
            # NAs allowed;  turn into sub-paths
            xy <- splitOnNA(x, y)
            x <- lapply(xy, "[[", "x")
            y <- lapply(xy, "[[", "y")
        } else {
            stop("'x' and 'y' must both be lists or both be atomic")
        }
        # If the path contains multiple sub-paths then NAs are NOT allowed
        # (to follow 'grid' behaviour)
        bad <- !all(sapply(x, function(z) { all(is.finite(z)) })) ||
            !all(sapply(y, function(z) { all(is.finite(z)) }))
        if (bad)
            stop("non-finite x or y in graphics path")
    }
    n <- length(x)
    d <- mapply(function(subx, suby) {
                    paste(paste(c("M",
                                  rep("L", length(subx) - 1)),
                                round(subx, 2), round(suby, 2),
                                collapse=" "),
                          "Z")
                }, x, y)

    tmpattr <- c(list(id = prefixName(id),
                      d = paste(unlist(d), collapse = " "),
                      "fill-rule" = switch(rule, winding="nonzero", "evenodd")),
                 svgStyleAttributes(style, svgdev),
                 svgAttribTxt(attributes, id, "path", svgdev))
    tmpattr <- attrList(tmpattr)


    has.link <- hasLink(links[id])
    if (has.link)
        svgStartLink(links[id], show[id], svgdev)

    newXMLNode("path", parent = svgDevParent(svgdev),
               attrs = tmpattr)

    if (has.link)
        svgEndLink(svgdev)
}

svgRaster <- function(x, y, width, height, angle=0, datauri, id=NULL,
                      just, vjust, hjust,
                      attributes=svgAttrib(), links=NULL, show=NULL,
                      style=svgStyle(), svgdev=svgDevice()) {
    # Draw nothing if non-finite location or size
    if (!(is.finite(x) && is.finite(y) &&
              is.finite(width) && is.finite(height)))
        return()
    has.link <- hasLink(links[id])
    if (has.link)
        svgStartLink(links[id], show[id], svgdev)

    if (width < 0) {
        x <- x + width # shifts x to the left
        width <- abs(width)
    }
    if (height < 0) {
        y <- y + height # shifts y down
        height <- abs(height)
    }

    rx <- round(x, 2)
    ry <- round(y, 2)

    transform <- paste0("translate(", rx, ", ", round(height + y, 2), ")")
    angleTransform <- svgAngleTransform(rx, ry, angle)
    if (!is.null(angleTransform)) {
        transform <- paste(angleTransform, transform)
    }
    attrlist <- c(list(id = prefixName(id),
                       transform = transform),
                  svgStyleAttributes(style, svgdev),
                  svgAttribTxt(attributes, id, "g", svgdev))
    attrlist <- attrList(attrlist)
    newXMLNode("g", parent = svgDevParent(svgdev),
               attrs = attrlist,
               newXMLNode("g",
                          attrs = list(id = paste(prefixName(id), "scale",
                                           sep = getSVGoption("id.sep")),
                              transform = paste0("scale(",
                                  round(width, 2), ", ",
                                  round(-height, 2), ")")),
                          newXMLNode("image",
                                     # Suppress the namespace warning because
                                     # we know in this specific case it is
                                     # a spurious warning
                                     suppressNamespaceWarning = TRUE,
                                     attrs = list(x = 0,
                                         y = 0,
                                         width = 1,
                                         height = 1,
                                         "xlink:href" = datauri,
                                         preserveAspectRatio = "none"))))

    if (has.link)
        svgEndLink(svgdev)
}

svgRect <- function(x, y, width, height, angle=0, id=NULL,
                    attributes=svgAttrib(), links=NULL, show=NULL,
                    style=svgStyle(), svgdev=svgDevice()) {
    # Draw nothing if non-finite location or size
    if (!(is.finite(x) && is.finite(y) &&
              is.finite(width) && is.finite(height)))
        return()
    has.link <- hasLink(links[id])
    if (has.link)
        svgStartLink(links[id], show[id], svgdev)
    
    if (width < 0) {
        x <- x + width # shifts x to the left
        width <- abs(width)
    }
    
    if (height < 0) {
        y <- y + height # shifts y down
        height <- abs(height)
    }
    
    rx <- round(x, 2)
    ry <- round(y, 2)
    
    attrlist <- c(list(id = prefixName(id),
                       x = rx,
                       y = ry,
                       width = round(width, 2),
                       height = round(height, 2),
                       transform = svgAngleTransform(rx, ry, angle)), 
                  svgStyleAttributes(style, svgdev),
                  svgAttribTxt(attributes, id, "rect", svgdev))
    attrlist <- attrList(attrlist)
    newXMLNode("rect", parent = svgDevParent(svgdev),
               attrs = attrlist)
    
    if (has.link)
        svgEndLink(svgdev)
}

svgTextSplitLines <- function(text, id, lineheight, charheight,
                              vjust, svgdev) {
    # Splitting based on linebreaks
    splitText <- strsplit(text, "\n")
    # If text is "", produces character(0), so fix that
    if (length(splitText[[1]]) == 0)
        splitText[[1]] <- ""

    n <- length(splitText[[1]])

    # Need to adjust positioning based on vertical justification.
    # Horizontal justification is done for us.
    # Only the first line needs to be modified, the rest are all
    # just one line below the previous line
    if (vjust %in% c("centre", "center"))
        firstDelta <- - ((lineheight * (n - 1) - charheight) / 2)
    if (vjust == "bottom")
        firstDelta <- - (n - 1) * lineheight
    if (vjust == "top")
        firstDelta <- charheight
    lineheight <- c(firstDelta, rep(lineheight, n - 1))

    textContent <- splitText[[1]]
    # Note that x=0 here so that we push it to the left, hjust
    # is worked out automatically from there
    for (i in 1:n) {
        newXMLNode("tspan", parent = svgDevParent(svgdev),
                   attrs = list(id = paste(id, "tspan", i,
                                    sep=getSVGoption("id.sep")),
                       dy = round(lineheight[i], 2),
                       x = 0),
                   newXMLTextNode(textContent[i]))
    }
}

svgTextElement <- function(text, id, rot, hjust, vjust,
                           lineheight, charheight, style, svgdev=svgDevice()) {
    # Rotation in SVG goes clockwise from +ve x=axis
    transform <- if (rot != 0)
                   list(transform = paste0("rotate(", round(-rot, 2), ")"))
                 else
                   NULL
    attrlist <- c(list(x = 0,
                       y = 0,
                       id = paste(id, "text", sep=getSVGoption("id.sep"))),
                  transform,
                  textAnchor(hjust),
                  svgStyleAttributes(style, svgdev))
    attrlist <- attrList(attrlist)
    newpar <- newXMLNode("text", parent = svgDevParent(svgdev),
                         attrs = attrlist)
    # Set parent of all <tspan>s to be the <text> el
    svgDevChangeParent(newpar, svgdev)
    # Write each of the lines here
    svgTextSplitLines(text, id, lineheight, charheight, vjust, svgdev)
    # Resetting parent
    svgDevChangeParent(xmlParent(newpar), svgdev)
}

# NOTE that the precise placement of math is even less likely to work
# than normal text.  Besides the problem of the browser using a
# different font (which is more likely because a math expression
# typically uses multiple fonts), the web browser will be using
# a different formula layout engine compared to R so things like
# the spacing between operators will be different.
# One particular problem is that R justifies math formulas
# relative to the bounding box of the formula, whereas it
# appears that Firefox at least justifies relative to the formula
# baseline (just from observation).
# The code below tries to do something rational by making use
# of finer detail metric information for the formula
# to mimic R's vertical justification.
svgMathElement <- function(text, id, rot, hjust, vjust,
                           width, height, ascent, descent,
                           lineheight, charheight, fontheight,
                           fontfamily, fontface, style,
                           svgdev=svgDevice()) {
    # Determine x/y based on width/height and hjust/vjust
    if (hjust %in% c("centre", "center"))
        x <- -width/2
    if (hjust == "left")
        x <- 0
    if (hjust == "right")
        x <- -width
    if (vjust %in% c("centre", "center"))
        y <- -(max(ascent, fontheight) + descent)/2
    if (vjust == "bottom")
        y <- -(max(ascent, fontheight) + descent)
    if (vjust == "top") {
        if (fontheight > ascent)
            y <- -(fontheight - ascent)
        else
            y <- (ascent - fontheight)
    }

    tmpattr <- c(list(x = round(x, 2),
                      y = round(y, 2),
                      id = paste(id, "mathtext", sep=getSVGoption("id.sep")),
                      width = round(3*width, 2),
                      height = round(3*height, 2)),
                 svgStyleAttributes(style, svgdev))
    if (rot != 0)
        tmpattr$transform <- paste0("rotate(", round(-rot, 2), ")")

    switch <- newXMLNode("switch", parent = svgDevParent(svgdev))
    foreignObj <- newXMLNode("foreignObject", parent = switch,
                             attrs = attrList(tmpattr))
    svgDevChangeParent(foreignObj, svgdev)
    expr2mml(text, fontfamily, fontface, svgdev)
    svgDevChangeParent(xmlParent(switch), svgdev)
}

svgText <- function(x, y, text, hjust="left", vjust="bottom", rot=0,
                    width=1, height=1, angle=0, ascent=1, descent=0,
                    lineheight=1, charheight=.8, fontheight=1,
                    fontfamily="sans", fontface="plain",
                    id=NULL, attributes=svgAttrib(), links=NULL, show=NULL,
                    style=svgStyle(), svgdev=svgDevice()) {
    # Draw nothing if x/y non-finite
    if (!(is.finite(x) && is.finite(y)))
        return()
    has.link <- hasLink(links[id])
    if (has.link)
        svgStartLink(links[id], show[id], svgdev)

    if (width < 0) {
        x <- x + width # shifts x to the left
        width <- abs(width)
    }
    if (height < 0) {
        y <- y + height # shifts y down
        height <- abs(height)
    }

    rx <- round(x, 2)
    ry <- round(y, 2)
    
    topattrs <- list()
    topattrs$id <- prefixName(id)
    angleTransform <- svgAngleTransform(rx, ry, angle)
    topattrs$transform <- paste0("translate(", rx, ", ", ry, ")")
    if (!is.null(angleTransform)) {
        topattrs$transform <- paste(angleTransform, topattrs$transform)
    }
    topattrs$`stroke-width` <- "0.1"
    topattrs <- c(topattrs, svgAttribTxt(attributes, id, "g", svgdev))

    # Flip the y-direction again so that text is drawn "upright"
    # Do the flip in a separate <g> so that can animate the
    # translation easily
    # Use a tspan to do the vertical alignment
    topg <- newXMLNode("g", parent = svgDevParent(svgdev),
                       attrs = topattrs)
    sec <- newXMLNode("g", parent = topg,
                      attrs = list(id = paste(prefixName(id), "scale",
                                       sep = getSVGoption("id.sep")),
                          transform = "scale(1, -1)"))

    # Let all child <tspan> elements or MathML fragments be
    # located under the *second* <g>
    svgDevChangeParent(sec, svgdev)

    if (is.language(text)) {
        svgMathElement(text, prefixName(id), rot, hjust, vjust,
                       width, height, ascent, descent,
                       lineheight, charheight, fontheight,
                       fontfamily, fontface, style,
                       svgdev)
    } else {
        svgTextElement(text, prefixName(id), rot, hjust, vjust,
                       lineheight, charheight, style,
                       svgdev)
    }

    # Reset parent to parent of entire text "grob"
    svgDevChangeParent(xmlParent(topg), svgdev)

    if (has.link)
        svgEndLink(svgdev)
}

svgCircle <- function(x, y, r, id=NULL,
                      attributes=svgAttrib(), links=NULL, show=NULL,
                      style=svgStyle(), svgdev=svgDevice()) {
    # Draw nothing if non-finite location or size
    if (!(is.finite(x) && is.finite(y) && is.finite(r)))
        return()
    has.link <- hasLink(links[id])
    if (has.link)
        svgStartLink(links[id], show[id], svgdev)
    
    tmpattr <- c(list(id = prefixName(id),
                      cx = round(x, 2),
                      cy = round(y, 2),
                      r = round(r, 2)),
                 svgStyleAttributes(style, svgdev),
                 svgAttribTxt(attributes, id, "circle", svgdev))
    tmpattr <- attrList(tmpattr)
    has.link <- hasLink(links[id])
    newXMLNode("circle", parent = svgDevParent(svgdev),
             attrs = tmpattr)
    
    if (has.link)
        svgEndLink(svgdev)
}

svgScript <- function(body, href, type="application/ecmascript",
                      id=NULL, svgdev=svgDevice()) {
  tmpattr <- list(type = type,
                  id = prefixName(getid(id, svgdev, 1)))
  if (nchar(href) > 0)
    tmpattr$`xlink:href` <- href

  script <- newXMLNode("script", parent = svgDevParent(svgdev),
                       attrs = tmpattr)

  if (nchar(body) > 0) {
    # "body" adds newlines because otherwise the CDATA delimiters are part
    # of the first and last line of text, break it apart to look nicer
    newXMLCDataNode(paste0("\n", body, "\n"),
                    parent = script)
  }
}

# Beginning of definition of all PCH elements
# Note that these definitions come ported from
# R's /src/main/engine.c
# Note in particular that radius is defined to be 0.375 * size
# so that width is 0.75 of the specified size. Most of the time
# this means we have a computed radius of 3.75

svgUseSymbol <- function(id, x, y, size, pch, angle=0,
                         attributes=svgAttrib(), links=NULL, show=NULL,
                         style=svgStyle(), svgdev=svgDevice()) {

    # Draw nothing if non-finite location or size
    if (!(is.finite(x) && is.finite(y) && is.finite(size)))
        return()
    has.link <- hasLink(links[id])
    if (has.link)
        svgStartLink(links[id], show[id], svgdev)

    # Ensure the "dot" is only 1px wide
    if (pch == ".")
        size <- 1

    # Ensure we refer to the correct <symbol> id
    numpch <- if (is.character(pch))
                  as.numeric(charToRaw(pch))
              else
                  pch
    
    rx <- round(x, 2)
    ry <- round(y, 2)
    
    tmpattr <- list(id = prefixName(id),
                    "xlink:href" =
                        paste0("#", prefixName(paste0("gridSVG.pch", numpch))),
                    x = rx, y = ry,
                    width = round(size, 2),
                    height = round(size, 2))

    # centering adjustment
    r <- round(-size / 2, 2)
    tmpattr$transform <- paste0("translate(", r, ",", r, ")")
    angleTransform <- svgAngleTransform(rx, ry, angle)
    if (!is.null(angleTransform)) {
        tmpattr$transform <- paste(angleTransform, tmpattr$transform)
    }
  
    # Preserve order
    tmpattr <- c(tmpattr,
                 svgStyleAttributes(style, svgdev),
                 svgAttribTxt(attributes, id, "use", svgdev))
    
    # Need to scale the stroke width otherwise for large points
    # we also have large strokes
    sw <- as.numeric(tmpattr$`stroke-width`)
    scalef <- size / 10 # 10 is the point viewBox size
    sw <- sw / scalef
    tmpattr$`stroke-width` <- round(sw, 2)

    # For pch outside 0-25 or characters
    if (is.character(pch) || (is.numeric(pch) && pch > 25)) {
        # When we have a "." we have a special case
        if ((is.character(pch) && pch == ".") ||
            (is.numeric(pch) && pch == 46)) {
            # Strip unnecessary attribs
            fsind <- which(names(tmpattr) == "font-size")
            if (length(fsind) > 0)
                tmpattr <- tmpattr[-fsind]
            # Because we really want just a dot, use crispEdges
            # as anti-aliasing isn't really necessary
            tmpattr$`shape-rendering` <- "crispEdges"
        } else {
            # Make the s-w small so we see a stroke just barely
            tmpattr$`stroke-width` <- "0.1"
            # Set the font-size, otherwise it's going to mess with our scaling.
            # 10px so it's the size of the point definition
            tmpattr$`font-size` <- "10"
        }
    }
    
    newXMLNode("use", parent = svgDevParent(svgdev),
               attrs = attrList(tmpattr))

    if (has.link)
        svgEndLink(svgdev)
}

# Dispatching function, simply following a naming scheme,
# somewhat nasty but works fine
svgPoint <- function(pch, svgdev = svgDevice()) {
  textpch <- FALSE
  if (is.character(pch)) {
    if (pch == ".")
      fnname <- "svgPointDot"
    else {
      fnname <- "svgPointChar"
      textpch <- TRUE
    }
  } else {
    fnname <- paste0("svgPoint", pch)
  }
  do.call(fnname, if (textpch) list(pch = pch, svgdev = svgdev)
                  else list(svgdev = svgdev))
}

# Special point, the dot
svgPointDot <- function(svgdev = svgDevice()) {
  newXMLNode("rect", parent = svgDevParent(svgdev),
             attrs = list(x = -0.5, y = -0.5,
                          width = 1, height = 1))
}

# Actual point character
svgPointChar <- function(pch, svgdev = svgDevice()) {
  # Transform to "flip" the text back
  newXMLNode("text", parent = svgDevParent(svgdev),
             attrs = list(x = 0, y = 0,
                          fontsize = 7.5,
                          transform = "scale(1, -1)",
                          "text-anchor" = "middle",
                          "baseline-shift" = "-25%"),
             newXMLTextNode(pch))
}

# S square
svgPoint0 <- function(svgdev = svgDevice()) {
    newXMLNode("rect", parent = svgDevParent(svgdev),
               attrs = list(x = -3.75, y = -3.75,
                            width = 7.5, height = 7.5))
}

# S octahedron (circle)
svgPoint1 <- function(svgdev = svgDevice()) {
    newXMLNode("circle", parent = svgDevParent(svgdev),
               attrs = list(cx = 0, cy = 0,
                            r = 3.75))
}

# S triangle - point up
svgPoint2 <- function(svgdev = svgDevice()) {
    TRC0 <- sqrt(4 * pi/(3 * sqrt(3)))
    TRC1 <- TRC0 * sqrt(3) / 2
    TRC2 <- TRC0 / 2
    r <- TRC0 * 3.75
    xc <- TRC1 * 3.75
    yc <- TRC2 * 3.75
    linexs <- round(c(0, xc, -xc, 0), 2)
    lineys <- round(c(r, -yc, -yc, r), 2)
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(linexs, lineys,
                                           sep = ",", collapse = " ")))
}

# S plus
svgPoint3 <- function(svgdev = svgDevice()) {
    xc <- sqrt(2) * 3.75
    yc <- sqrt(2) * 3.75
    l1xs <- round(c(-xc, xc), 2)
    l1ys <- c(0, 0)
    l2xs <- c(0, 0)
    l2ys <- round(c(-yc, yc), 2)
    # Horizontal
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(l1xs, l1ys,
                                           sep = ",", collapse = " ")))
    # Vertical
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(l2xs, l2ys,
                                           sep = ",", collapse = " ")))
}

# S times
svgPoint4 <- function(svgdev = svgDevice()) {
    xc <- 3.75
    yc <- 3.75
    l1xs <- c(-xc, xc)
    l1ys <- c(-yc, yc)
    l2xs <- c(-xc, xc)
    l2ys <- c(yc, -yc)
    # /
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(l1xs, l1ys,
                                           sep = ",", collapse = " ")))
    # \
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(l2xs, l2ys,
                                           sep = ",", collapse = " ")))
}

# S diamond
svgPoint5 <- function(svgdev = svgDevice()) {
    xc <- sqrt(2) * 3.75
    yc <- sqrt(2) * 3.75
    linexs <- round(c(-xc, 0, xc, 0, -xc), 2)
    lineys <- round(c(0, yc, 0, -yc, 0), 2)
    newXMLNode("polygon", parent = svgDevParent(svgdev),
               attrs = list(points = paste(linexs, lineys,
                                           sep = ",", collapse = " ")))
}

# S triangle - point down
svgPoint6 <- function(svgdev = svgDevice()) {
    TRC0 <- sqrt(4 * pi/(3 * sqrt(3)))
    TRC1 <- TRC0 * sqrt(3) / 2
    TRC2 <- TRC0 / 2
    r <- TRC0 * 3.75
    xc <- TRC1 * 3.75
    yc <- TRC2 * 3.75
    linexs <- round(c(0, xc, -xc, 0), 2)
    lineys <- round(c(-r, yc, yc, -r), 2)
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(linexs, lineys,
                                           sep = ",", collapse = " ")))
}

# S square and times superimposed
svgPoint7 <- function(svgdev = svgDevice()) {
    svgPoint0(svgdev)
    svgPoint4(svgdev)
}

# S plus and times superimposed
svgPoint8 <- function(svgdev = svgDevice()) {
    svgPoint3(svgdev)
    svgPoint4(svgdev)
}

# S diamond and plus superimposed
svgPoint9 <- function(svgdev = svgDevice()) {
    svgPoint3(svgdev)
    svgPoint5(svgdev)
}

# S hexagon (circle) and plus superimposed
svgPoint10 <- function(svgdev = svgDevice()) {
    newXMLNode("circle", parent = svgDevParent(svgdev),
               attrs = list(cx = 0, cy = 0,
                            r = 3.75))
    l1xs <- c(-3.75, 3.75)
    l1ys <- c(0, 0)
    l2xs <- c(0, 0)
    l2ys <- c(-3.75, 3.75)
    # Horizontal
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(l1xs, l1ys,
                                           sep = ",", collapse = " ")))
    # Vertical
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(l2xs, l2ys,
                                           sep = ",", collapse = " ")))
}

# S superimposed triangles
svgPoint11 <- function(svgdev = svgDevice()) {
    TRC0 <- sqrt(4 * pi/(3 * sqrt(3)))
    TRC1 <- TRC0 * sqrt(3) / 2
    TRC2 <- TRC0 / 2
    xc <- 3.75
    r <- TRC0 * xc
    yc <- TRC2 * xc
    yc <- 0.5 * (yc + r)
    xc <- TRC1 * xc

    # Pointing down
    linexs <- round(c(0, xc, -xc, 0), 2)
    lineys <- round(c(-r, yc, yc, -r), 2)
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(linexs, lineys,
                                           sep = ",", collapse = " ")))

    # Pointing up
    linexs <- round(c(0, xc, -xc, 0), 2)
    lineys <- round(c(r, -yc, -yc, r), 2)
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(linexs, lineys,
                                           sep = ",", collapse = " ")))
}

# S square and plus superimposed
svgPoint12 <- function(svgdev = svgDevice()) {
    svgPoint0(svgdev)
    l1xs <- c(-3.75, 3.75)
    l1ys <- c(0, 0)
    l2xs <- c(0, 0)
    l2ys <- c(-3.75, 3.75)
    # Horizontal
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(l1xs, l1ys,
                                           sep = ",", collapse = " ")))
    # Vertical
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(l2xs, l2ys,
                                           sep = ",", collapse = " ")))
}

# S octagon (circle) and times superimposed
svgPoint13 <- function(svgdev = svgDevice()) {
    svgPoint1(svgdev)
    svgPoint4(svgdev)
}

# S square and point-*down* triangle superimposed
# Note: R source refers to this as being point-up
svgPoint14 <- function(svgdev = svgDevice()) {
    r <- 3.75
    xs <- c(0, r, -r, 0)
    ys <- c(-r, r, r, -r)
    newXMLNode("polyline", parent = svgDevParent(svgdev),
               attrs = list(points = paste(xs, ys, sep = ",",
                                           collapse = " ")))
    newXMLNode("rect", parent = svgDevParent(svgdev),
               attrs = list(x = -r, y = -r,
                            width = 2*r, height = 2*r))
}

# S filled square
svgPoint15 <- function(svgdev = svgDevice()) {
    svgPoint0(svgdev)
}

# S filled octagon (circle)
svgPoint16 <- function(svgdev = svgDevice()) {
    svgPoint1(svgdev)
}

# S filled point-up triangle
svgPoint17 <- function(svgdev = svgDevice()) {
    svgPoint2(svgdev)
}

# S filled diamond
svgPoint18 <- function(svgdev = svgDevice()) {
    svgPoint5(svgdev)
}

# R filled circle
svgPoint19 <- function(svgdev = svgDevice()) {
    svgPoint1(svgdev)
}

# R `Dot' (small circle)
svgPoint20 <- function(svgdev = svgDevice()) {
    newXMLNode("circle", parent = svgDevParent(svgdev),
               attrs = list(cx = 0, cy = 0,
                            r = 2.5))
}

# circles
svgPoint21 <- function(svgdev = svgDevice()) {
    svgPoint1(svgdev)
}

# squares
svgPoint22 <- function(svgdev = svgDevice()) {
    r <- round(sqrt(pi / 4) * 3.75, 2)
    newXMLNode("rect", parent = svgDevParent(svgdev),
               attrs = list(x = -r, y = -r,
                            width = 2*r, height = 2*r))
}

# diamonds
svgPoint23 <- function(svgdev = svgDevice()) {
    r <- 3.75 * sqrt(pi / 4) * sqrt(2)
    xs <- round(c(-r, 0, r, 0, -r), 2)
    ys <- round(c(0, r, 0, -r, 0), 2)
    newXMLNode("polygon", parent = svgDevParent(svgdev),
               attrs = list(points = paste(xs, ys,
                                           sep = ",", collapse = " ")))
}

# triangle (point up)
svgPoint24 <- function(svgdev = svgDevice()) {
    svgPoint2(svgdev)
}

# triangle (point down)
svgPoint25 <- function(svgdev = svgDevice()) {
    svgPoint6(svgdev)
}

#############
# Internal functions
#############

# SVG Devices
# A device is an environment so that we can modify values
# stored within it.
# Store a list of transformation functions for
# x, y, width, and height;  this will allow viewports
# to be defined within user coordinates (see svgPushViewport
# and svgPopViewport)

svgDevice <- function(width=200, height=200, strict=TRUE, rootAttrs) {
  dev <- new.env(FALSE, emptyenv())
  assign("width", width, envir=dev)
  assign("height", height, envir=dev)
  assign("parent", NULL, envir=dev)
  assign("strict", strict, envir=dev)
  assign("rootAttrs", rootAttrs, envir=dev)
  assign("id", 1, envir=dev)
  return(dev)
}

svgDevWidth <- function(svgdev) {
  get("width", envir=svgdev)
}

svgDevHeight <- function(svgdev) {
  get("height", envir=svgdev)
}

svgStrict <- function(svgdev) {
    get("strict", envir=svgdev)
}

svgRootAttrs <- function(svgdev) {
    get("rootAttrs", envir=svgdev)
}

svgDevParent <- function(svgdev) {
  get("parent", envir=svgdev)
}

svgDevChangeParent <- function(newpar, svgdev) {
  assign("parent", newpar, envir=svgdev)
}

getid <- function(id, svgdev, n=1) {
  if (is.null(id))
    svgID(svgdev) + (1:n - 1)
  else {
    if (n > 1)
      paste(id, 1:n, sep="")
    else
      id
  }
}

svgID <- function(svgdev) {
  get("id", envir=svgdev)
}

hasLink <- function(link) {
  ! (is.null(link) || is.na(link))
}

incID <- function(svgdev, n=1) {
  assign("id", get("id", envir=svgdev) + n, envir=svgdev)
}

svgHeader <- function(width, height, svgdev=svgDevice()) {
    # This header tested on standalone SVG file in Firefox 3
    # FIXME:  add default xmlns for animation and scripts too?
    attrs <- svgRootAttrs(svgdev)
    # Resolve any conflict between rootAttrs[c("width", "height")]
    # and 'width' and 'height' passed in (from device size)
    if (!"width" %in% names(attrs)) {
        attrs$width <- paste0(round(width, 2), "px")
    }
    if (!"height" %in% names(attrs)) {
        attrs$height <- paste0(round(height, 2), "px")
    }
    if (!"viewBox" %in% names(attrs)) {
        attrs$viewBox <- paste(0, 0, round(width, 2), round(height, 2))
    }
    if (!"version" %in% names(attrs)) {
        attrs$version <- "1.1"
    }
    # Give the <svg> element an ID only if there is a prefix AND
    # 'rootAttrs' does not specify an ID
    if (!"id" %in% names(attrs)) {
        if (nzchar(get("prefix", envir = .gridSVGEnv))) {
            attrs$id <- get("prefix", envir = .gridSVGEnv)
        }
    }
    if (svgStrict(svgdev)) {
        checkAttrs(attrs, "svg")
    }
    svgdoc <-
        newXMLDoc(namespaces = list("http://www.w3.org/2000/svg",
                      xlink = "http://www.w3.org/1999/xlink"),
                  node =  newXMLNode("svg", attrs = attrList(attrs),
                      namespaceDefinitions = list("http://www.w3.org/2000/svg",
                          xlink = "http://www.w3.org/1999/xlink")))
    # Invert the y-axis so that y and height values measure "up"
    rootg <- newXMLNode("g",
                        parent = xmlRoot(svgdoc),
                        attrs = list(transform = paste0("translate(0, ",
                                         round(svgDevHeight(svgdev), 2),
                                         ") scale(1, -1)")))
    svgDevChangeParent(rootg, svgdev)
}

# SVG attributes
svgAttrib <- function(...) {
  temp <- list(...)
  if (length(temp) == 0)
    list()
  else if (is.null(temp[[1]]))
    list()
  else
    temp
}

# Removes NULL values and flattens our attrib list
# so we can include lists as elements in "alist"
# and arrive at a flattened list
# ALSO removes names from attr values so they don't corrupt attr names
# e.g., avoid list(a=c(b=1)) becoming c(a.b=1)
attrList <- function(alist) {
  as.list(unlist(lapply(alist, unname)))
}

listToSVGAttrib <- function(alist) {
  alist
}

emptyAttrib <- function(attributes) {
  length(attributes) == 0
}

svgClassList <- function(classes) {
    if (is.null(classes) || ! get("addClasses", envir = .gridSVGEnv))
        list()
    else
        list(class = paste0(unique(classes), collapse = " "))
}

checkAttrs <- function(attrList, eltName) {
    attrTable <- get("genAttrTable", .gridSVGEnv)
    validAttrs <- attrTable$attribute[attrTable$element == eltName]
    names <- names(attrList)
    svgnames <- names %in% validAttrs
    if (!all(svgnames)) {
        warning(paste("Removing non-SVG attribute name(s):",
                      paste(names[!svgnames], collapse=", ")))
        attrList[names[svgnames]]
    } else {
        attrList
    }
}

# Only use the attributes that are for this 'id'
svgAttribTxt <- function(attributes, id, eltName, svgdev) {
    if (emptyAttrib(attributes)) {
        list()
    } else {
        attributes <- lapply(attributes,
                             function(attr, id) {
                                 kept <- attr[names(attr) == id]
                                 if (length(kept) == 0)
                                     NULL
                                 else
                                     kept
                             },
                             id)
        # Drop NULL attributes
        attributes <- attributes[!sapply(attributes, is.null)]
        # Need to wipe out names because it messes things up when we
        # need to create an attribute list for nodes
        if (length(attributes) > 0) {
            attrs <- lapply(attributes,
                            function(x) {
                                names(x) <- NULL
                                x
                            })
            if (svgStrict(svgdev)) {
                checkAttrs(attrs, eltName)
            } else {
                attrs
            }
        } else {
            list()
        }
    }
}

# SVG styling
svgStyle <- function(...) {
  list(...)
}

listToSVGStyle <- function(alist) {
  alist
}

emptyStyle <- function(svgstyle) {
  length(svgstyle) == 0
}

svgStyleCSS <- function(svgstyle) {
  if (emptyStyle(svgstyle)) {
    ""
  } else {
      paste('style="',
            do.call("paste",
                    c(mapply(function(name, style) {
                        paste(name, ":", style, sep="")
                    }, names(svgstyle), svgstyle),
                      list(sep="; "))),
            '"', sep="")
    # paste('style="', paste(names(svgstyle), ":",
    #                        paste(svgstyle), sep="", collapse="; "),
    #       '"', sep="")
  }
}

# SVG version of gpar()

# NOTE using SVG presentation attributes
# RATHER THAN CSS style attribute
# BECAUSE ...
# - can modify single presentation attribute without affecting
#   other presentation attributes (if CSS style then have to
#   reset the entire thing) and can do this from JavaScript.
# - presentation attributes have lower priority than CSS style
#   SO this allows overriding by specifying CSS style later.
#   Can also override with general style sheet later.
svgStyleAttributes <- function(svgstyle, dev) {
    if (emptyStyle(svgstyle)) {
        list()
    } else {
        # Remove non-SVG style attributes (with warning)
        names <- names(svgstyle)
        svgnames <- names %in% get("SVGParList", envir=.gridSVGEnv)
        if (svgStrict(dev) && any(!svgnames)) {
            warning(paste("Removing non-SVG style attribute name(s):",
                          paste(names[!svgnames], collapse=", ")))
            svgstyle <- svgstyle[names[svgnames]]
        }
        if (any(sapply(svgstyle, length) > 1))
            stop("All SVG style attribute values must have length 1")
        svgstyle
    }
}

# Specifying text justification
textAnchor <- function(hjust) {
  list("text-anchor" =
       switch(hjust,
              left="start",
              center="middle",
              centre="middle",
              right="end",
              "start"))
}

dominantBaseline <- function(vjust) {
  list("dominant-baseline" =
       switch(vjust,
              bottom="auto",
              center="middle",
              centre="middle",
              top="text-top",
              "baseline"))
}

baselineShift <- function(vjust) {
  list("baseline-shift" =
       switch(vjust,
              bottom="0%",
              center="-50%",
              centre="-50%",
              top="-100%",
              "0%"))
}

alignmentBaseline <- function(vjust) {
  list("alignment-baseline" =
       switch(vjust,
              baseline="baseline",
              bottom="bottom",
              center="middle",
              centre="middle",
              top="top",
              "baseline"))
}

splitOnNA <- function(x, y) {
    # Check for non-finite (rather than just NA)
    nas <- !is.finite(x) | !is.finite(y)
    wnas <- which(nas)
    N <- max(length(x), length(y))
    start <- c(1, pmin(wnas + 1, N))
    end <- c(pmax(wnas - 1, 1), N)
    xylist <- mapply(function(s, e) {
                         if (e > s) list(x=x[s:e], y=y[s:e])
                     },
                     start, end, SIMPLIFY=FALSE)
    xylist[!sapply(xylist, is.null)]
}

