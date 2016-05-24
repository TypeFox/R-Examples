grid.patternFill <- function(path, pattern = NULL, label = NULL,
                             alpha = 1, group = TRUE, redraw = FALSE,
                             strict = FALSE, grep = FALSE, global = FALSE) {
    if (is.null(label) & is.null(pattern)) {
        stop("At least one of 'label' or 'pattern' must be supplied")
    } else if (is.null(label)) {
        label <- getNewLabel("gridSVG.patternFill")
        registerPatternFill(label, pattern)
        pattern <- NULL # use the ref from now on
    } else if (is.null(pattern)) {
        checkForDefinition(label)
    } else {
        checkExistingDefinition(label)
        registerPatternFill(label, pattern)
        pattern <- NULL # use the ref from now on
    }

    grobApply(path, function(path) {
        grid.set(path, patternFillGrob(grid.get(path), pattern = pattern,
                                       label = label, alpha = alpha,
                                       group = group),
                 redraw = redraw)
    }, strict = strict, grep = grep, global = global)

    invisible()
}

patternFillGrob <- function(x, pattern = NULL, label = NULL,
                            alpha = 1, group = TRUE) {
    if (is.null(label) & is.null(pattern)) {
        stop("At least one of 'label' or 'pattern' must be supplied")
    } else if (is.null(label)) {
        label <- getNewLabel("gridSVG.patternFill")
        registerPatternFill(label, pattern)
    } else if (is.null(pattern)) {
        checkForDefinition(label)
    } else {
        checkExistingDefinition(label)
        registerPatternFill(label, pattern)
    }

    if (length(alpha) != length(label))
        alpha <- rep(alpha, length.out = length(label))

    x$referenceLabel <- c(x$referenceLabel, label)
    # Attribs to be garnished *at draw time*. In particular needs to be
    # done because the label ID is not known until then, because of things
    # like prefixes and separators.
    x$patternFillLabel <- label
    x$patternFillAlpha <- alpha
    x$patternFillGroup <- group
    class(x) <- unique(c("patternFilled.grob", class(x)))
    x
}

pattern <- function(grob = NULL,
                    x = unit(0, "npc"), y = unit(0, "npc"),
                    width = unit(0.1, "npc"), height = unit(0.1, "npc"),
                    default.units = "npc",
                    just = "centre", hjust = NULL, vjust = NULL,
                    dev.width = 7, dev.height = 7) {
    if (! is.unit(x))
        x <- unit(x, default.units)
    if (! is.unit(y))
        y <- unit(y, default.units)
    if (! is.unit(width))
        width <- unit(width, default.units)
    if (! is.unit(height))
        height <- unit(height, default.units)

    pattern <- list(grob = grob,
                    x = x, y = y,
                    width = width, height = height,
                    just = just, hjust = hjust, vjust = vjust,
                    dev.width = dev.width, dev.height = dev.height)
    class(pattern) <- "pattern"
    pattern
}

registerPatternFill <- function(label, pattern = NULL, ...) {
    checkExistingDefinition(label)
    refDefinitions <- get("refDefinitions", envir = .gridSVGEnv)

    if (is.null(pattern)) {
        pattern <- gridSVG::pattern(...)
    } else if (! inherits(pattern, "pattern")) {
        stop("'pattern' must be a 'pattern' object")
    }

    if (is.null(pattern$grob))
        stop("A grob must be given for a fill pattern definition")

    # Now convert *at time of definition* to absolute units (inches)
    loc <- leftbottom(pattern$x, pattern$y, pattern$width, pattern$height,
                      pattern$just, pattern$hjust, pattern$vjust, NULL)
    x <- loc$x
    y <- loc$y
    width <- convertWidth(pattern$width, "inches")
    height <- convertHeight(pattern$height, "inches")

    # ID will be overwritten later, because we might change
    # the separator used for "id.sep"
    defList <- list(
        label = label,
        id = getID(label, "ref"),
        grob = pattern$grob,
        vp = getAbsoluteVp(),
        x = x,
        y = y,
        width = width,
        height = height,
        dev.width = pattern$dev.width,
        dev.height = pattern$dev.height
    )

    class(defList) <- "patternFillDef"

    refDefinitions[[label]] <- defList
    assign("refDefinitions", refDefinitions, envir = .gridSVGEnv)
    assign("refUsageTable",
           rbind(get("refUsageTable", envir = .gridSVGEnv),
                 data.frame(label = label, used = FALSE,
                            stringsAsFactors = FALSE)),
           envir = .gridSVGEnv)

    # Return NULL invisibly because we don't actually care what the
    # definition looks like until gridSVG tries to draw it. 
    invisible()
}

registerPatternFillRef <- function(label, refLabel, pattern = NULL, ...) {
    checkExistingDefinition(label)
    refDefinitions <- get("refDefinitions", envir = .gridSVGEnv)
    if (! refLabel %in% names(refDefinitions))
        stop(paste("The reference labelled", sQuote(label), "does not exist."))

    if (is.null(pattern)) {
        pattern <- gridSVG::pattern(...)
    } else if (! inherits(pattern, "pattern")) {
        stop("'pattern' must be a 'pattern' object")
    }

    # Now convert *at time of definition* to absolute units (inches)
    offsets <- getAbsoluteOffset()
    loc <- leftbottom(pattern$x, pattern$y,
                      pattern$width, pattern$height,
                      pattern$just, pattern$hjust, pattern$vjust, NULL)
    x <- loc$x + offsets[1]
    y <- loc$y + offsets[2]
    width <- convertWidth(pattern$width, "inches")
    height <- convertHeight(pattern$height, "inches")

    defList <- list(
        label = label,
        refLabel = refLabel,
        id = getID(label, "ref"),
        x = x,
        y = y,
        width = width,
        height = height
    )

    class(defList) <- "patternFillRefDef"

    refDefinitions[[label]] <- defList
    assign("refDefinitions", refDefinitions, envir = .gridSVGEnv)
    assign("refUsageTable",
           rbind(get("refUsageTable", envir = .gridSVGEnv),
                 data.frame(label = label, used = FALSE,
                            stringsAsFactors = FALSE)),
           envir = .gridSVGEnv)

    # Return NULL invisibly because we don't actually care what the
    # definition looks like until gridSVG tries to draw it. 
    invisible()
}

primToDev.patternFilled.grob <- function(x, dev) {
    setLabelUsed(x$referenceLabel)
    label <- getLabelID(x$patternFillLabel)
    # Allowing fill-opacity to be set by a garnish because
    # grid only knows about a colour and its opacity. If we use a
    # reference instead of a then nothing is known about the opacity.
    # We want to ensure that we can still set it, so use the garnish
    # to overwrite it.
    pg <- garnishGrob(x, fill = paste0("url(#", label, ")"),
                      "fill-opacity" = x$patternFillAlpha,
                      group = x$patternFillGroup)
    # Now need to remove all pattern fill appearances in the class list.
    # This is safe because repeated pattern filling just clobbers existing
    # attributes.
    cl <- class(pg)
    class(pg) <- cl[cl != "patternFilled.grob"]
    primToDev(pg, dev)
}

drawDef.patternFillDef <- function(def, dev) {
    svgdev <- dev@dev

    # Convert grid coords to SVG coords
    x <- round(cx(def$x, dev), 2)
    y <- round(cy(def$y, dev), 2)
    width <- round(cw(def$width, dev), 2)
    height <- round(ch(def$height, dev), 2)

    # Checking for flipped scales
    if (width < 0) {
        x <- x + width # shifts x to the left
        width <- abs(width)
    }
    if (height < 0) {
        y <- y + height # shifts y down
        height <- abs(height)
    }

    # Attempt to use a known-safe prefix
    # If the prefix is safe, then it will *always* be safe
    # because the names are known *after* content is drawn
    # and the referenced labels must be unique
    prefix <- paste0("gridSVG.patternFill.", def$id)

    # There is a little bit of replication going on from
    # 'grid.export' but it avoids some problems.
    # We could use 'grid.export' recursively but we lose the ability to
    # definitely generate unique names if that is the case because usage
    # tables would be wiped.
    # A viewport and gTree are forced to ensure everything is unique because
    # we want paths to be used.
    # We do not care at this point whether it is strictly necessary to
    # perform all of this because we just want unique IDs.
    olddev <- dev.cur()
    pdf(file = NULL, width = def$dev.width, height = def$dev.height)
        newdev <- openSVGDev("", res = dev@res,
                             width = def$dev.width, height = def$dev.height)
        pushViewport(viewport(name = getID(prefix, "vp")))
        grid.draw(gTree(name = getID(prefix, "grob"),
                  children = gList(grid.force(def$grob)),
                  gp = get.gpar())) # Force gp to ensure correct styling
        grid.force(redraw = FALSE)
        gt <- grid.grab(name = "gridSVG", wrap = TRUE)
        gridToDev(gt, newdev)
        newroot <- devClose(newdev)
        viewBox <- xmlGetAttr(newroot, "viewBox")
        gridSVGNode <- prefixName("gridSVG")
        # Clone this node so that when the pattern definition refers to
        # namespaces (e.g. rasterGrobs have 'xlink:href'), the namespaces
        # are not destroyed when this temporary device is closed.
        gridSVGNode <-
            xmlClone(getNodeSet(newroot,
                                paste0("//*[@id='", gridSVGNode, "']"))[[1]])
    dev.off()
    dev.set(olddev)

    # Creating the pattern element
    pattern <- newXMLNode("pattern",
                          attrs = list(id = def$id, x = x, y = y,
                                       width = width, height = height,
                                       viewBox = viewBox,
                                       patternUnits = "userSpaceOnUse"),
                          parent = svgDevParent(svgdev))
    # Assigning its children
    xmlChildren(pattern) <- xmlChildren(gridSVGNode)
}

drawDef.patternFillRefDef <- function(def, dev) {
    svgdev <- dev@dev

    # Convert grid coords to SVG coords
    x <- round(cx(def$x, dev), 2)
    y <- round(cy(def$y, dev), 2)
    width <- round(cw(def$width, dev), 2)
    height <- round(ch(def$height, dev), 2)

    # Checking for flipped scales
    if (width < 0) {
        x <- x + width # shifts x to the left
        width <- abs(width)
    }
    if (height < 0) {
        y <- y + height # shifts y down
        height <- abs(height)
    }

    # Creating the pattern element
    pattern <- newXMLNode("pattern",
        attrs = list(id = def$id, x = x, y = y,
                     width = width, height = height,
                     "xlink:href" =
                       paste0("#", getLabelID(def$refLabel))),
        parent = svgDevParent(svgdev))
}

# Ensure the patterns are retained on a forced grob
forceGrob.patternFilled.grob <- function(x) {
    y <- NextMethod()
    if (inherits(y, "forcedgrob")) {
        y$referenceLabel <- x$referenceLabel
        y$patternFillLabel <- x$patternFillLabel
        y$patternFillAlpha <- x$patternFillAlpha
        y$patternFillGroup <- x$patternFillGroup
        class(y) <- unique(c("patternFilled.grob", class(y)))
    }
    y
}
