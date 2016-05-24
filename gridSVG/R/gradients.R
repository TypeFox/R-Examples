# High level functions for applying gradients as fills to grobs
grid.gradientFill <- function(path, gradient = NULL, label = NULL,
                              alpha = 1, group = TRUE, redraw = FALSE,
                              strict = FALSE, grep = FALSE, global = FALSE) {
    if (is.null(gradient) & is.null(label)) {
        stop("At least one of 'gradient' or 'label' must be supplied")
    } else if (is.null(label)) {
        label <- getNewLabel("gridSVG.gradientFill")
        registerGradientFill(label, gradient)
        gradient <- NULL # use the ref from now on
    } else if (is.null(gradient)) {
        checkForDefinition(label)
    } else {
        checkExistingDefinition(label)
        registerGradientFill(label, gradient)
        gradient <- NULL # use the ref from now on
    }

    grobApply(path, function(path) {
        grid.set(path, gradientFillGrob(grid.get(path), gradient = gradient,
                                        label = label, alpha = alpha,
                                        group = group),
                 redraw = redraw)
    }, strict = strict, grep = grep, global = global)

    invisible()
}

gradientFillGrob <- function(x, gradient = NULL, label = NULL,
                             alpha = 1, group = TRUE) {
    if (is.null(gradient) & is.null(label)) {
        stop("At least one of 'gradient' or 'label' must be supplied")
    } else if (is.null(label)) {
        label <- getNewLabel("gridSVG.gradientFill")
        registerGradientFill(label, gradient)
    } else if (is.null(gradient)) {
        checkForDefinition(label)
    } else {
        checkExistingDefinition(label)
        registerGradientFill(label, gradient)
    }

    if (length(alpha) != length(label))
        alpha <- rep(alpha, length.out = length(label))

    x$referenceLabel <- c(x$referenceLabel, label)
    x$gradientFillLabel <- label
    x$gradientFillAlpha <- alpha
    x$gradientFillGroup <- group
    class(x) <- unique(c("gradientFilled.grob", class(x)))
    x
}

linearGradient <- function(col = c("black", "white"),
                           stops = seq(0, 1, length.out = length(col)),
                           gradientUnits = c("bbox", "coords"),
                           x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                           y0 = unit(0, "npc"), y1 = unit(1, "npc"),
                           default.units = "npc",
                           spreadMethod = c("pad", "reflect", "repeat")) {
    # Vectorising colours & stops
    nstops <- max(length(col), length(stops))
    col <- rep(col, length.out = nstops)
    stops <- rep(stops, length.out = nstops)

    offset <- round(stops, 2)
    stopCol <- sapply(col, function(x) devColToSVG(x), USE.NAMES = FALSE)
    stopOpacity <- devColAlphaToSVG(col2rgb(col, alpha = TRUE)[4, ])

    gradientUnits <- match.arg(gradientUnits)
    spreadMethod <- match.arg(spreadMethod)
    
    if (! is.unit(x0))
        x0 <- unit(x0, default.units)
    if (! is.unit(x1))
        x1 <- unit(x1, default.units)
    if (! is.unit(y0))
        y0 <- unit(y0, default.units)
    if (! is.unit(y1))
        y1 <- unit(y1, default.units)

    # Convert gradientUnits to SVG values
    gradientUnits <- switch(gradientUnits,
                            bbox = "objectBoundingBox",
                            coords = "userSpaceOnUse")

    # Need to get npc-like values from units
    if (gradientUnits == "objectBoundingBox") {
        # Convert to npc 
        x0 <- convertX(x0, "npc", valueOnly = TRUE)
        x1 <- convertX(x1, "npc", valueOnly = TRUE)
        y0 <- convertY(y0, "npc", valueOnly = TRUE)
        y1 <- convertY(y1, "npc", valueOnly = TRUE)
    }

    grad <- list(element = "linearGradient",
                 gradientUnits = gradientUnits,
                 x1 = x0, x2 = x1,
                 y1 = y0, y2 = y1,
                 spreadMethod = spreadMethod,
                 offset = offset, stopCol = stopCol,
                 stopOpacity = stopOpacity)
    class(grad) <- c("linear.gradient", "gradient")
    grad
}

radialGradient <- function(col = c("black", "white"),
                           stops = seq(0, 1, length.out = length(col)),
                           gradientUnits = c("bbox", "coords"),
                           x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                           r = unit(0.5, "npc"),
                           fx = unit(0.5, "npc"), fy = unit(0.5, "npc"),
                           default.units = "npc",
                           spreadMethod = c("pad", "reflect", "repeat")) {
    # Vectorising colours & stops
    nstops <- max(length(col), length(stops))
    col <- rep(col, length.out = nstops)
    stops <- rep(stops, length.out = nstops)

    offset <- round(stops, 2)
    stopCol <- sapply(col, function(x) devColToSVG(x), USE.NAMES = FALSE)
    stopOpacity <- devColAlphaToSVG(col2rgb(col, alpha = TRUE)[4, ])

    gradientUnits <- match.arg(gradientUnits)
    spreadMethod <- match.arg(spreadMethod)
    if (is.null(stops))
        stops <- list()
    
    if (! is.unit(x))
        x <- unit(x, default.units)
    if (! is.unit(y))
        y <- unit(y, default.units)
    if (! is.unit(r))
        r <- unit(r, default.units)
    if (! is.unit(fx))
        fx <- unit(fx, default.units)
    if (! is.unit(fy))
        fy <- unit(fy, default.units)

    # Convert gradientUnits to SVG values
    gradientUnits <- switch(gradientUnits,
                            bbox = "objectBoundingBox",
                            coords = "userSpaceOnUse")

    # Need to get npc-like values from units
    if (gradientUnits == "objectBoundingBox") {
        x <- convertX(x, "npc", valueOnly = TRUE)
        y <- convertY(y, "npc", valueOnly = TRUE)

        rw <- convertWidth(r, "npc", valueOnly = TRUE)
        rh <- convertHeight(r, "npc", valueOnly = TRUE)
        r <- pmin(rw, rh)

        fx <- convertX(fx, "npc", valueOnly = TRUE)
        fy <- convertY(fy, "npc", valueOnly = TRUE)
    }

    grad <- list(element = "radialGradient",
                 gradientUnits = gradientUnits,
                 cx = x, cy = y, r = r,
                 fx = fx, fy = fy,
                 spreadMethod = spreadMethod,
                 offset = offset, stopCol = stopCol,
                 stopOpacity = stopOpacity)
    class(grad) <- c("radial.gradient", "gradient")
    grad
}

print.gradient <- function(x, ...) {
    prln <- function(label, value) {
        cat(sprintf(paste0(label, ": %s\n"), value))
    }
    prln("Type", x$element)
    n <- length(x$offset)
    prln("Number of stops", n)
    cat("\n")
    prln("Gradient stops", "")
    for (i in 1:n) {
        cat("  ")
        cat("Offset:", x$offset[i])
        cat("  ")
        cat("Colour:", x$stopCol[i])
        cat("  ")
        cat("Opacity:", x$stopOpacity[i])
        cat("\n")
    }
    invisible(x)
}

flattenLinearGradient <- function(gradient) {
    # Flatten all locations here
    if (gradient$gradientUnits == "userSpaceOnUse") {
        offsets <- getAbsoluteOffset()
        width <- convertWidth(gradient$x2 - gradient$x1, "inches",
                              valueOnly = TRUE)
        height <- convertHeight(gradient$y2 - gradient$y1, "inches",
                                valueOnly = TRUE)
        gradient$x1 <- convertX(gradient$x1, "inches") + offsets[1]
        gradient$x2 <- convertX(gradient$x2, "inches") + offsets[1]
        gradient$y1 <- convertY(gradient$y1, "inches") + offsets[2]
        gradient$y2 <- convertY(gradient$y2, "inches") + offsets[2]
    }
    gradient
}

flattenRadialGradient <- function(gradient) {
    # Flatten all locations here
    if (gradient$gradientUnits == "userSpaceOnUse") {
        offsets <- getAbsoluteOffset()
        gradient$cx <- convertX(gradient$cx, "inches") + offsets[1]
        gradient$cy <- convertY(gradient$cy, "inches") + offsets[2]
        gradient$r <- dToInches(gradient$r, NULL)
        gradient$fx <- convertX(gradient$fx, "inches") + offsets[1]
        gradient$fy <- convertY(gradient$fy, "inches") + offsets[2]
    }
    gradient
}

registerGradientFill <- function(label, gradient) {
    checkExistingDefinition(label)
    
    # Flattening all locations
    gradient <-
        if (inherits(gradient, "radial.gradient"))
            flattenRadialGradient(gradient)
        else
            flattenLinearGradient(gradient)

    gradient$label <- label
    gradient$id <- getID(label, "ref")
    class(gradient) <- "gradientDef"

    refDefinitions <- get("refDefinitions", envir = .gridSVGEnv)
    refDefinitions[[label]] <- gradient
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

svgLinearGradient <- function(def, dev) {
    svgdev <- dev@dev

    # Convert grid coords to SVG coords if we are using coordinates
    # rather than the bounding box of the referring object
    if (def$gradientUnits == "userSpaceOnUse") {
        def$x1 <- cx(def$x1, dev)
        def$x2 <- cx(def$x2, dev)
        def$y1 <- cy(def$y1, dev)
        def$y2 <- cy(def$y2, dev)
    }

    gradient <- newXMLNode("linearGradient",
        parent = svgDevParent(svgdev),
        attrs = list(id = def$id,
                     x1 = round(def$x1, 2), x2 = round(def$x2, 2),
                     y1 = round(def$y1, 2), y2 = round(def$y2, 2),
                     gradientUnits = def$gradientUnits,
                     spreadMethod = def$spreadMethod))

    svgDevChangeParent(gradient, svgdev)
}

svgRadialGradient <- function(def, dev) {
    svgdev <- dev@dev

    # Convert grid coords to SVG coords if we are using coordinates
    # rather than the bounding box of the referring object
    if (def$gradientUnits == "userSpaceOnUse") {
        def$cx <- cx(def$cx, dev)
        def$cy <- cy(def$cy, dev)
        def$r <- cd(def$r, dev)
        def$fx <- cx(def$fx, dev)
        def$fy <- cy(def$fy, dev)
    }

    gradient <- newXMLNode("radialGradient",
        parent = svgDevParent(svgdev),
        attrs = list(id = def$id,
                     cx = round(def$cx, 2), cy = round(def$cy, 2),
                     r = round(def$r, 2),
                     fx = round(def$fx, 2), fy = round(def$fy, 2),
                     gradientUnits = def$gradientUnits,
                     spreadMethod = def$spreadMethod))

    svgDevChangeParent(gradient, svgdev)
}

primToDev.gradientFilled.grob <- function(x, dev) {
    setLabelUsed(x$referenceLabel)
    label <- getLabelID(x$gradientFillLabel)
    # Allowing fill-opacity to be set by a garnish because
    # grid only knows about a colour and its opacity. If we use a
    # reference instead of a then nothing is known about the opacity.
    # We want to ensure that we can still set it, so use the garnish
    # to overwrite it.
    gf <- garnishGrob(x, fill = paste0("url(#", label, ")"),
                      "fill-opacity" = x$gradientFillAlpha,
                      group = x$gradientFillGroup)
    # Now need to remove all gradient fill appearances in the class list.
    # This is safe because repeated gradient filling just clobbers existing
    # attributes.
    cl <- class(gf)
    class(gf) <- cl[cl != "gradientFilled.grob"]
    primToDev(gf, dev)
}

drawDef.gradientDef <- function(def, dev) {
    svgdev <- dev@dev

    if (def$element == "linearGradient")
        svgLinearGradient(def, dev)
    else
        svgRadialGradient(def, dev)

    # Adding the gradient stops
    for (i in 1:length(def$offset)) {
        newXMLNode("stop",
                   attrs = list(offset = def$offset[i],
                                "stop-color" = def$stopCol[i],
                                "stop-opacity" = def$stopOpacity[i]),
                   parent = svgDevParent(svgdev))
    }

    # Going back up from the stops to the parent of the gradient
    svgDevChangeParent(xmlParent(svgDevParent(svgdev)), svgdev)
}

# Ensure the gradient fill is retained on a forced grob
forceGrob.gradientFilled.grob <- function(x) {
    y <- NextMethod()
    if (inherits(y, "forcedgrob")) {
        y$referenceLabel <- x$referenceLabel
        y$gradientFillLabel <- x$gradientFillLabel
        y$gradientFillAlpha <- x$gradientFillAlpha
        y$gradientFillGroup <- x$gradientFillGroup
        class(y) <- unique(c("gradientFilled.grob", class(y)))
    }
    y
}
