viewportCreate <- function(vpname, newname = NULL,
                           vpPath.sep = getSVGoption("vpPath.sep")) {
    coords <- gridSVGCoords()
    if (is.null(coords))
        stop("gridSVGCoords() must be initialised")
    rootvp <- coords$ROOT
    if (is.null(rootvp))
        stop("the ROOT viewport must have coords info set")
    targetvp <- coords[[vpname]]
    if (is.null(vpname))
        stop(paste("the viewport", sQuote(vpname), "must have coords info set, see", sQuote("gridSVGCoords")))
    # Avoid having a vpPath as a viewport name
    if (is.null(newname)) {
        splitname <- strsplit(vpname, vpPath.sep)[[1]]
        vpname <- tail(splitname, 1)
    } else {
        vpname <- newname
    }

    npcx <- targetvp$x / rootvp$width
    npcy <- targetvp$y / rootvp$height
    npcw <- targetvp$width / rootvp$width
    npch <- targetvp$height / rootvp$height
    viewport(x = unit(npcx, "npc"), y = unit(npcy, "npc"),
             width = unit(npcw, "npc"), height = unit(npch, "npc"),
             angle = targetvp$angle,
             just = c("left", "bottom"), name = vpname,
             xscale = targetvp$xscale, yscale = targetvp$yscale)
}

viewportConvertX <- function(vpname, x, from, to = "svg") {
    currCoords <- validCoordsInfo(vpname)
    vpCoords <- currCoords[[vpname]]
    if (vpCoords$angle != 0) {
        stop("Viewport angle non-zero; consider using viewportConvertPos()")
    }
    if (from == "svg")
        x <- x - vpCoords$x
    width <- convertx(vpCoords, x, from, to)
    if (to == "svg")
        width <- width + vpCoords$x
    width
}

viewportConvertY <- function(vpname, x, from, to = "svg") {
    currCoords <- validCoordsInfo(vpname)
    vpCoords <- currCoords[[vpname]]
    if (vpCoords$angle != 0) {
        stop("Viewport angle non-zero; consider using viewportConvertPos()")
    }
    if (from == "svg")
        x <- x - vpCoords$y
    height <- converty(vpCoords, x, from, to)
    if (to == "svg")
        height <- height + vpCoords$y
    height
}

viewportConvertPos <- function(vpname, x, y, from, to = "svg") {
    currCoords <- validCoordsInfo(vpname)
    vpCoords <- currCoords[[vpname]]
    if (from == "svg") {
        x <- x - vpCoords$x        
        y <- y - vpCoords$y
    }
    width <- convertx(vpCoords, x, from, to)
    height <- converty(vpCoords, y, from, to)
    if (vpCoords$angle != 0) {
        theta <- -vpCoords$angle/180*pi
        w <- cos(theta)*width + sin(theta)*height
        h <- -sin(theta)*width + cos(theta)*height
        width <- w
        height <- h
    }
    if (to == "svg") {
        width <- width + vpCoords$x
        height <- height + vpCoords$y
    }
    list(x=width, y=height)
}

viewportConvertWidth <- function(vpname, x, from, to) {
    currCoords <- validCoordsInfo(vpname)
    vpCoords <- currCoords[[vpname]]
    convertx(vpCoords, x, from, to, FALSE)
}

viewportConvertHeight <- function(vpname, x, from, to) {
    currCoords <- validCoordsInfo(vpname)
    vpCoords <- currCoords[[vpname]]
    converty(vpCoords, x, from, to, FALSE)
}

convertx <- function(vpCoords, x, from, to, loc=TRUE) {
    i <- toInches(from, x,
                  vpCoords$width, vpCoords$xscale, vpCoords$inch, loc)
    u <- toUnit(to, i,
                vpCoords$width, vpCoords$xscale, vpCoords$inch, loc)
    u 
}

converty <- function(vpCoords, x, from, to, loc=TRUE) {
    i <- toInches(from, x,
                  vpCoords$height, vpCoords$yscale, vpCoords$inch, loc)
    u <- toUnit(to, i,
                vpCoords$height, vpCoords$yscale, vpCoords$inch, loc)
    round(u, 2)
}

viewportConvertDim <- function(vpname, w, h, from, to) {
    currCoords <- validCoordsInfo(vpname)
    vpCoords <- currCoords[[vpname]]
    wi <- toInches(from, w,
                   vpCoords$width, vpCoords$xscale, vpCoords$inch, FALSE)
    hi <- toInches(from, h,
                   vpCoords$height, vpCoords$yscale, vpCoords$inch, FALSE)
    if (vpCoords$angle != 0) {
        theta <- -vpCoords$angle/180*pi
        w <- cos(theta)*wi + sin(theta)*hi
        h <- -sin(theta)*wi + cos(theta)*hi
        wi <- w
        hi <- h
    }
    wu <- toUnit(to, wi,
                 vpCoords$width, vpCoords$xscale, vpCoords$inch, FALSE)
    hu <- toUnit(to, hi,
                 vpCoords$height, vpCoords$yscale, vpCoords$inch, FALSE)
    list(w=wu, h=hu)
}

toInches <- function(from, unitValue,
                     vpDimSize, nativeScale, dimInchSize, loc) {
    if (from == "inches")
        return(unitValue)
    
    nativeToInches <- function(nativeValue, nativeScale, vpDimSize,
                               dimInchSize, loc) {
        if (loc) {
            dist <- nativeValue - nativeScale[1]
        } else {
            dist <- nativeValue
        }
        nativeUnitSize <- vpDimSize / abs(nativeScale[2] - nativeScale[1])
        dist * nativeUnitSize / dimInchSize
    }
  
    npcToInches <- function(npcValue, vpDimSize, dimInchSize) {
        (npcValue * vpDimSize) / dimInchSize
    }

    if (from == "native") {
        result <- nativeToInches(unitValue, nativeScale, vpDimSize, dimInchSize,
                                 loc)
    } else if (from == "npc") {
        result <- npcToInches(unitValue, vpDimSize, dimInchSize)
    } else if (from == "svg") {
        result <- unitValue / dimInchSize
    } else {
        result <- convertUnit(unit(unitValue, from), "inches", valueOnly = TRUE)
    }

    result
}

toUnit <- function(to, unitValue, vpDimSize, nativeScale, dimInchSize, loc) {
    if (to == "inches")
        return(unitValue)

    inchesToNative <- function(inchesValue, nativeScale, vpDimSize,
                               dimInchSize, loc) {
        npc <- (inchesValue * dimInchSize) / vpDimSize
        vpRange <- nativeScale[2] - nativeScale[1]
        if (loc) {
            (npc * vpRange) + nativeScale[1]
        } else {
            (npc * vpRange)
        }
    }
  
    inchesToNpc <- function(inchesValue, vpDimSize, dimInchSize) {
        (inchesValue * dimInchSize) / vpDimSize
    }

    if (to == "native") {
        result <- inchesToNative(unitValue, nativeScale, vpDimSize, dimInchSize,
                                 loc)
    } else if (to == "npc") {
        result <- inchesToNpc(unitValue, vpDimSize, dimInchSize)
    } else if (to == "svg") {
        result <- unitValue * dimInchSize
    } else {
        result <- convertUnit(unit(unitValue, "inches"), to, valueOnly = TRUE)
    }

    result
}

