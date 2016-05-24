# This file is concerned with the use of objects that need to be
# referenced, e.g. pattern fills, gradient fills, filters, etc.

drawDef <- function(def, dev) {
    UseMethod("drawDef")
}

# This function ensures that if we change the ID separator value between
# the time of definition and draw time, we still get the expected IDs.
assignRefIDs <- function() {
    refdefs <- get("refDefinitions", envir = .gridSVGEnv)
    for (i in seq_along(refdefs))
        refdefs[[i]]$id <- getLabelID(refdefs[[i]]$label)
    assign("refDefinitions", refdefs, envir = .gridSVGEnv)

    # Because the separators might have changed, ensure that the
    # usageTable has the correct (escaped) values for selectors and xpath
    ut <- get("usageTable", envir = .gridSVGEnv)
    inds <- which(ut$type == "ref")
    for (i in inds) {
        fullName <- paste(ut[i, "name"],
                          ut[i, "suffix"],
                          sep = getSVGoption("id.sep"))
        sel <- prefixName(escapeSelector(fullName))
        xp <- prefixName(escapeXPath(fullName))
        ut[i, "selector"] <- sel
        ut[i, "xpath"] <- xp
    }
    assign("usageTable", ut, envir = .gridSVGEnv)
}

totalDefinitions <- function() {
    # This function should only be called for calculating progress bar
    # length. In addition, it will only be called if definitions are
    # required to be flushed so this should always return a non-zero
    # result.
    sum(get("refUsageTable", envir = .gridSVGEnv)[, "used"])
}

flushDefinitions <- function(dev) {
    svgdev <- dev@dev

    # Keep copies of old tables because they will be modified when
    # we draw any children
    usageTable <- get("usageTable", envir = .gridSVGEnv)
    vpCoords <- get("vpCoords", envir = .gridSVGEnv)
    use.vpPaths <- get("use.vpPaths", envir = .gridSVGEnv)
    use.gPaths <- get("use.gPaths", envir = .gridSVGEnv)
    uniqueNames <- get("uniqueNames", envir = .gridSVGEnv)

    # Set required options -- we don't care about what the user
    # has specified at this point because they shouldn't be
    # touching any reference definitions
    assign("use.vpPaths", TRUE, envir = .gridSVGEnv)    
    assign("use.gPaths", TRUE, envir = .gridSVGEnv)
    assign("uniqueNames", TRUE, envir = .gridSVGEnv)
    
    # Begin creating definitions
    # First ensure we're under #gridSVG
    rootID <- prefixName("gridSVG")
    gridSVGNode <- getNodeSet(xmlRoot(svgDevParent(svgdev)),
                              paste0("//*[@id='", rootID, "']"))[[1]]
    svgDevChangeParent(gridSVGNode, svgdev)
    refDefinitions <- get("refDefinitions", envir = .gridSVGEnv)
    pchUsageTable <- get("pchUsageTable", envir = .gridSVGEnv)
    if (! length(refDefinitions) && ! any(pchUsageTable[, "used"]))
        return() # fast path for leaving early, avoids creating a <defs>
    defs <- newXMLNode("defs", parent = svgDevParent(svgdev), at = 0)
    svgDevChangeParent(defs, svgdev)

    # Check whether we have any dependent references, e.g. have a pattern
    # fill by reference in use but not the pattern itself. We need to ensure
    # that both are written out.
    n <- length(refDefinitions)
    rut <- get("refUsageTable", envir = .gridSVGEnv)
    for (i in seq_len(n)) {
        def <- refDefinitions[[i]]
        if (isLabelUsed(def$label) &&
            ! is.null(def$refLabel) && ! isLabelUsed(def$refLabel))
            rut[rut$label == def$refLabel, "used"] <- TRUE
    }

    # Now trying to find out if there are trees of referenced content.
    for (i in seq_len(n)) {
        used <- labelsUsed(refDefinitions[[i]])
        if (is.null(used))
            next
        flaggedLabels <- used %in% rut$label
        if (any(flaggedLabels))
            rut[flaggedLabels, "used"] <- TRUE
    }

    assign("refUsageTable", rut, envir = .gridSVGEnv) 
    # Now can work out how many defs to flush
    ndefs <- totalDefinitions()
    progressInit("defs", ndefs)

    # Now try drawing
    for (i in seq_len(n)) {
        def <- refDefinitions[[i]]
        upViewport(0)
        if (! is.null(def$vp))
            pushViewport(def$vp)
        if (isLabelUsed(def$label))
            drawDef(def, dev)
        upViewport(0)
        progressStep("defs", ndefs)
    }

    # Resetting to original values
    assign("vpCoords", vpCoords, envir = .gridSVGEnv)
    assign("use.vpPaths", use.vpPaths, envir = .gridSVGEnv)
    assign("use.gPaths", use.gPaths, envir = .gridSVGEnv)
    assign("uniqueNames", uniqueNames, envir = .gridSVGEnv)

    # All of the points that have been used in the image will now be flushed.
    # This is done after any of the other references primarily because
    # a pattern could use a pch but we want the definition of the pch to
    # appear beforehand.
    flushPchs(dev)

    # Reset ref usage table
    rut <- get("refUsageTable", envir = .gridSVGEnv)
    rut$used <- logical(nrow(rut))
    assign("refUsageTable", rut, envir = .gridSVGEnv)
    # Again for usage table 
    assign("usageTable", usageTable, envir = .gridSVGEnv)

    # Get out of defs
    svgDevChangeParent(xmlParent(defs), svgdev)
}

flushPchs <- function(dev) {
    pchUsageTable <- get("pchUsageTable", envir = .gridSVGEnv)
    if (! any(pchUsageTable[, "used"]))
        return()
    usedPchs <- pchUsageTable[pchUsageTable[, "used"] > 0, "pch"]
    progressInit("pch", length(usedPchs))
    # Reversing so that when we insert at the start of the <defs>
    # the pchs are ordered from small to big, not big to small.
    # This is purely cosmetic.
    for (pch in rev(usedPchs)) {
        asciipch <- if (pch %in% 32:127)
                        rawToChar(as.raw(pch))
                    else
                        pch
        devStartSymbol(pch, dev)
        devPoint(asciipch, dev)
        devEndSymbol(dev)
        progressStep("pch")
    }
}

anyRefsDefined <- function() {
    ut <- get("usageTable", envir = .gridSVGEnv)
    nrow(ut) > 0 && any(ut$type == "ref")
}

# Methods used for grabbing the list of references used by a definition.
# Particularly useful in the case of patterns where it could contain
# content which also references other content. In other words, it allows
# us to be able to get a tree of dependencies, rather than just a flat
# list.
labelsUsed <- function(x) {
    UseMethod("labelsUsed")
}

labelsUsed.patternFillRefDef <- function(x) {
    NULL
}

labelsUsed.filterDef <- function(x) {
    NULL
}

labelsUsed.gradientDef <- function(x) {
    NULL
}

labelsUsed.patternFillDef <- function(x) {
    labelsUsed(x$grob)
}

labelsUsed.maskDef <- function(x) {
    labelsUsed(x$grob)
}

labelsUsed.clipPathDef <- function(x) {
    labelsUsed(x$grob)
}

labelsUsed.grob <- function(x) {
    x$referenceLabel
}

labelsUsed.gTree <- function(x) {
    c(x$referenceLabel, unlist(lapply(x$children, labelsUsed)))
}

# Used for knowing whether to write out a definition.
# If a definition has not been used we do not write it out.
# If it has been used more than once we do not repeat the
# definition.
isLabelUsed <- function(label) {
    rut <- get("refUsageTable", envir = .gridSVGEnv)
    any(rut[rut$label %in% label, "used"])
}

setLabelUsed <- function(label) {
    if (! is.null(label) && length(label)) {
        rut <- get("refUsageTable", envir = .gridSVGEnv)
        if (any(rut$label %in% label)) {
            rut[rut$label %in% label, "used"] <- TRUE
            assign("refUsageTable", rut, envir = .gridSVGEnv)
            # Need to ensure that nested dependencies are also handled.
            # e.g. a mask definition that is filtered needs to trigger a
            #      filter to be drawn.
            refdefs <- get("refDefinitions", envir = .gridSVGEnv)
            for (i in seq_along(label)) {
                def <- refdefs[[label[i]]]
                setLabelUsed(labelsUsed(def))
            }
        } else {
            stop("An attempt was made to reference content that no longer exists.")
        }
    }
}

# Convenience function to list all referenced content definitions
listSVGDefinitions <- function(print = TRUE) {
    refdefs <- get("refDefinitions", envir = .gridSVGEnv)
    n <- length(refdefs)
    if (!n)
        return(invisible())
    defs <- data.frame(label = character(n),
                       type = character(n),
                       refLabel = character(n),
                       stringsAsFactors = FALSE)

    for (i in 1:n) {
        curdef <- refdefs[[i]]
        defs$label[i] <- curdef$label
        if (! is.null(curdef$refLabel))
            defs$refLabel[i] <- curdef$refLabel
        defs$type[i] <- switch(class(curdef)[1],
                               clipPathDef = "Clipping Path",
                               filterDef = "Filter Effect",
                               gradientDef = "Gradient Fill",
                               maskDef = "Mask",
                               patternFillDef = "Pattern Fill",
                               patternFillRefDef = "Pattern Fill Reference",
                               "")
    }

    if (print) {
        orderedTypes <- sort(unique(defs$type))
        indent <- "  "
        cat("Reference Definitions\n")
        for (i in 1:length(orderedTypes)) {
            typesub <- defs[defs$type == orderedTypes[i], ]
            cat("\n", orderedTypes[i], "s\n", sep = "")
            for (j in 1:nrow(typesub)) {
                cat(indent, typesub$label[j], sep = "")
                # If this is a pattern fill, show us what we're referencing
                if (nchar(typesub$refLabel[j]))
                    cat(" ", paste0("(referencing ", typesub$refLabel[j], ")"),
                        "\n", sep = "")
                else
                    cat("\n")
            }
        }
    }

    invisible(defs)
}

checkForDefinition <- function(label) {
    if (! all(label %in% names(get("refDefinitions", envir = .gridSVGEnv))))
        stop("A reference definition must be created before using its label")
}

checkExistingDefinition <- function(label) {
    if (any(label %in% names(get("refDefinitions", envir = .gridSVGEnv))))
        stop(paste("A label already exists as a reference definition")) 
}

# When we need to generate a temporary label (i.e. when specifying a
# gradient fill directly on a grob with no label), we supply a prefix and
# return a new label that is going to be unique (among labels).
# getID() will perform the task of ensuring uniqueness among IDs.
getNewLabel <- function(prefix) {
    i <- 1
    candidateName <- paste0(prefix, ".", i)
    refdefs <- get("refDefinitions", envir = .gridSVGEnv)
    while(candidateName %in% names(refdefs)) {
        i <- i + 1
        candidateName <- paste0(prefix, getSVGoption("id.sep"), i)
    }
    candidateName
}

getLabelID <- function(label) {
    ut <- get("usageTable", envir = .gridSVGEnv)
    suffix <- ut[ut$name %in% label & ut$type == "ref", "suffix"]
    prefixName(paste0(label, getSVGoption("id.sep"), suffix))
}


# This function allows us to collect a viewport that, when pushed into
# from the ROOT viewport, should allow us to draw in the same drawing
# environment that we called this function from.
getAbsoluteVp <- function(vp = current.viewport(),
                          tm = current.transform()) {
  transloc <- c(0, 0, 1) %*% tm
  loc <- (transloc / transloc[3])[-3]

  viewport(x = unit(loc[1], "inches"),
           y = unit(loc[2], "inches"),
           width = convertWidth(unit(1, "npc"), "inches"),
           height = convertHeight(unit(1, "npc"), "inches"),
           xscale = vp$xscale,
           yscale = vp$yscale,
           gp = get.gpar(),
           just = c("left", "bottom"),
           angle = vp$angle)
}

# Used for the case where we want to work out locations in absolute units
# *but* do not want to leave the current viewport. Useful when flattening
# locations for things such as mask areas.
getAbsoluteOffset <- function(tm = current.transform()) {
    transloc <- c(0, 0, 1) %*% tm
    loc <- (transloc / transloc[3])[-3]
    unit(loc, "inches")
}
