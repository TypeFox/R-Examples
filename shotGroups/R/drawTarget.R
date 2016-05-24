## return target definition including ring radii
getTarget <-
function(x, unit="cm", dstTarget=100, conversion="m2cm") {
    UseMethod("getTarget")
}

getTarget.character <-
function(x, unit="cm", dstTarget=100, conversion="m2cm") {
    if(!(x %in% names(shotGroups::targets))) {
        stop("Target unknown, see help(targets) for a list")
    }

    x <- shotGroups::targets[[x]]
    NextMethod("getTarget")
}

getTarget.default <-
function(x, unit="cm", dstTarget=100, conversion="m2cm") {
    unit <- match.arg(unit, choices=c("cm", "mm", "m", "in", "ft", "yd",
                                      "deg", "MOA", "SMOA", "rad", "mrad", "mil"))

    ## add ring radii
    if(!is.null(x$convert$ringW)) {
        x$convert$ringR <- with(x$convert,
                                c(ringD10i/2, ringD10/2,  # ring 10 with inner sub-division
                                  (ringD10/2) + seq(ringW, ringW*(x$nRings-1), by=ringW)))
    }

    ## add ring center vertical offset if present
    if(!is.null(x$convert$ringWV)) {
        x$convert$ringRV <- with(x$convert,
                                 c(ringD10Vi/2, ringD10V/2,  # ring 10 with inner sub-division
                                   (ringD10V/2) + seq(ringWV, ringWV*(x$nRings-1), by=ringWV)))
    }

    ## infer distance unit from conversion
    unitDst <- getUnits(conversion, first=TRUE)  # unit for distance

    ## convert all available measurements (ring width and radii, etc.)
    x$unitConv <- unit                           # add unit converted to
    x$inUnit <- if(unit %in% c("deg", "MOA", "SMOA", "rad", "mrad", "mil"))  {  # angular size
        ringConv <- with(x, paste0(unitDst, "2", unitTarget))
        with(x, Map(getMOA, convert, conversion=ringConv,
                    dst=dstTarget, type=unit))
    } else {                                     # absolute size
        ringConv <- with(x, getConvFac(paste0(unitTarget, "2", unit)))
        with(x, Map(function(y, rc) { rc*y }, y=convert, rc=ringConv))
    }

    x                                    # return selected target
}

## draw a target
drawTarget <-
function(x, unit="cm", dstTarget=100, conversion="m2cm", add=FALSE, cex=par("cex")) {
    unit <- match.arg(unit, choices=c("cm", "mm", "m", "in", "ft", "yd",
                                      "deg", "MOA", "SMOA", "rad", "mrad", "mil"))

    ## get chosen target including measures converted to unit
    target <- getTarget(x, unit=unit, dstTarget=dstTarget, conversion=conversion)

    ## open plot or add to existing plot?
    if(!add) {
        if(!is.null(target$inUnit$ringR)) {
            xLims <- with(target$inUnit, range(c(-ringR, ringR)))
            yLims <- if(!is.null(target$inUnit$ringRV)) {
                with(target$inUnit, range(c(-ringRV, ringRV)))
            } else {
                with(target$inUnit, range(c(-ringR, ringR)))
            }
        } else if(!is.null(target$inUnit$ringD8)) {
            xLims <- with(target$inUnit, range(c(-ringD8/2,  ringD8/2)))
            yLims <- with(target$inUnit, range(c(-ringD8V/2, ringD8V/2)))
        } else {
            warning("Cannot determine plot limits")
            xLims <- c(-1, 1)
            yLims <- c(-1, 1)
        }
        plot(0, 0, type="n", xlim=xLims, ylim=yLims, xlab=NA, ylab=NA, asp=1)
    }                                    # if(add)

    ## do we need a special drawing function?
    if(!is.null(target$draw)) {
        fun <- eval(target$draw)
        fun(target, cex=cex)
    } else {
        drawTarget_default(target, cex=cex)
    }

    return(invisible(target))
}

## draw a target
drawTarget_default <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x$inUnit, symbols(rep(0, length(ringR)), rep(0, length(ringR)), add=TRUE,
                           circles=rev(ringR), bg=rev(x$cols), fg=rev(x$colsTxt),
                           inches=FALSE))
     # abline(v=0, h=0, col="lightgray")    # add point of aim

    ## add ring numbers except for bullseye (ring number maxVal)
    ## bullseye has inner sub-division -> start numbers on ring 3
    rings1 <- with(x, seq(from=maxVal-nRings+1, to=maxVal-1, length.out=nRings-1)) # left side of center
    rings2 <- c(rings1, rev(rings1))                         # both sides
    pos1   <- with(x$inUnit, ringR[3:length(ringR)] - (ringW/2)) # right side of center
    pos2   <- c(-rev(pos1), pos1)                            # both sides
    cols1  <- with(x$inUnit, x$colsTxt[3:length(ringR)])     # right side of center
    cols2  <- c(rev(cols1), cols1)                           # both sides

    text(pos2, 0, cex=cex, label=rings2, col=cols2)
    text(0, pos2, cex=cex, label=rings2, col=cols2)
}

#####-----------------------------------------------------------------------
## ISSF special targets
drawTarget_ISSF25mRF <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x$inUnit, symbols(rep(0, length(ringR)), rep(0, length(ringR)), add=TRUE,
                           circles=rev(ringR), bg=rev(x$cols), fg=rev(x$colsTxt),
                           inches=FALSE))
    with(x$inUnit, rect(-ringR[length(ringR)],  -extra/2,
                        -ringR[4]+(0.15*ringW),  extra/2, col="white", border=NA))
    with(x$inUnit, rect( ringR[4]-(0.15*ringW), -extra/2,
                         ringR[length(ringR)],   extra/2, col="white", border=NA))
    # abline(v=0, h=0, col="lightgray")    # add point of aim

    ## add ring numbers except for bullseye (ring number maxVal)
    ## bullseye has inner sub-division -> start numbers on ring 3
    rings1 <- with(x, seq(from=maxVal-nRings+1, to=maxVal-1, length.out=nRings-1)) # left side of center
    rings2 <- c(rings1, rev(rings1))                         # both sides
    pos1   <- with(x$inUnit, ringR[3:length(ringR)] - (ringW/2)) # right side of center
    pos2   <- c(-rev(pos1), pos1)                            # both sides
    cols1  <- with(x$inUnit, x$colsTxt[3:length(ringR)])     # right side of center
    cols2  <- c(rev(cols1), cols1)                           # both sides

    text(0, pos2, cex=cex, label=rings2, col=cols2)
}

drawTarget_ISSF50m <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x$inUnit, symbols(rep(0, length(ringR)), rep(0, length(ringR)), add=TRUE,
                           circles=rev(ringR), bg=rev(x$cols), fg=rev(x$colsTxt),
                           inches=FALSE))
    with(x$inUnit, symbols(0, 0, add=TRUE, circles=extra, bg=x$cols[1], fg=NA,
                           inches=FALSE))
    with(x$inUnit, symbols(rep(0, length(ringR)-3), rep(0, length(ringR)-3), add=TRUE,
                           circles=rev(ringR)[-(1:3)], bg=rev(x$cols)[-(1:3)],
                           fg=rev(x$colsTxt)[-(1:3)], inches=FALSE))

    # abline(v=0, h=0, col="lightgray")    # add point of aim

    ## add ring numbers 1-8
    ## bullseye has inner sub-division -> start numbers on ring 4
    rings1 <- with(x, seq(from=maxVal-nRings+1, to=maxVal-2, length.out=nRings-2)) # left side of center
    rings2 <- c(rings1, rev(rings1))                         # both sides
    pos1   <- with(x$inUnit, ringR[4:length(ringR)] - (ringW/2)) # right side of center
    pos2   <- c(-rev(pos1), pos1)                            # both sides
    cols1  <- with(x$inUnit, x$colsTxt[4:length(ringR)])     # right side of center
    cols2  <- c(rev(cols1), cols1)                           # both sides

    text(pos2, 0, cex=cex, label=rings2, col=cols2)
    text(0, pos2, cex=cex, label=rings2, col=cols2)
}

## draw ISSF 300m target
drawTarget_ISSF300m <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x$inUnit, symbols(rep(0, length(ringR)), rep(0, length(ringR)), add=TRUE,
                           circles=rev(ringR), bg=rev(x$cols), fg=rev(x$colsTxt),
                           inches=FALSE))
    # abline(v=0, h=0, col="lightgray")    # add point of aim

    ## add ring numbers except for bullseye (ring number maxVal)
    ## bullseye has inner sub-division -> start numbers on ring 3
    rings <- with(x, seq(from=maxVal-nRings+1, to=maxVal-1, length.out=nRings-1)) # left side of center
    cols  <- with(x$inUnit, x$colsTxt[length(ringR):3])          # right side of center
    pos   <- with(x$inUnit, -ringR[length(ringR):3] + (ringW/2)) # right side of center
    ang   <- with(x, txtRot)
    angX  <- with(x, txtRot*pi/180)
    posLU <- cbind(pos, 0) %*% with(x, cbind(c(cos(  angX), sin(  angX)), c(-sin(  angX), cos(  angX))))
    posRU <- cbind(pos, 0) %*% with(x, cbind(c(cos(3*angX), sin(3*angX)), c(-sin(3*angX), cos(3*angX))))
    posRL <- cbind(pos, 0) %*% with(x, cbind(c(cos(5*angX), sin(5*angX)), c(-sin(5*angX), cos(5*angX))))
    posLL <- cbind(pos, 0) %*% with(x, cbind(c(cos(7*angX), sin(7*angX)), c(-sin(7*angX), cos(7*angX))))

    text(posLU[ , 1], posLU[ , 2], cex=cex, label=rings, col=cols, srt=  ang)
    text(posRU[ , 1], posRU[ , 2], cex=cex, label=rings, col=cols, srt=7*ang)
    text(posRL[ , 1], posRL[ , 2], cex=cex, label=rings, col=cols, srt=5*ang)
    text(posLL[ , 1], posLL[ , 2], cex=cex, label=rings, col=cols, srt=3*ang)
}

#####-----------------------------------------------------------------------
## NRA HPR targets
drawTarget_NRA_HPR <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x$inUnit, symbols(rep(0, length(ringR)), rep(0, length(ringR)), add=TRUE,
                           circles=rev(ringR), bg=rev(x$cols), fg=rev(x$colsTxt),
                           inches=FALSE))
    # abline(v=0, h=0, col="lightgray")    # add point of aim

    ## add ring numbers
    ## bullseye has inner sub-division -> start numbers on ring 2
    rings1 <- with(x, seq(from=maxVal-nRings+1, to=maxVal, length.out=nRings)) # left side of center
    outRng <- 1:(length(rings1)-2)
    outPos <- 3:length(rings1)
    rings2 <- c(rings1[outRng], rev(rings1[outRng]))         # both sides
    pos1   <- with(x$inUnit, ringR[2:length(ringR)] - (ringW/2)) # right side of center
    pos2   <- c(-rev(pos1[outPos]), pos1[outPos])            # both sides
    cols1  <- with(x$inUnit, x$colsTxt[2:length(ringR)])     # right side of center
    cols2  <- c(rev(cols1[outPos]), cols1[outPos])           # both sides

    text(pos2, 0,      cex=cex, label=rings2,      col=cols2)
    text(0, pos1[1:2], cex=cex, label=rev(rings1)[1:2], col=cols1[1:2])
    text(0, 0, cex=cex, label="X", col=cols1[1])
}

#####-----------------------------------------------------------------------
## BDMP special targets
drawTarget_BDMP3 <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x$inUnit, symbols(rep(0, length(ringR)), rep(0, length(ringR)), add=TRUE,
                           circles=rev(ringR), bg=rev(x$cols), fg=rev(x$colsTxt),
                           inches=FALSE))
    with(x$inUnit, rect(-extra/2, -extra/2, extra/2, extra/2, col="black", border=NA))
    with(x$inUnit, symbols(rep(0, length(ringR)), rep(0, length(ringR)), add=TRUE,
                           circles=rev(ringR), bg=NA, fg=rev(x$colsTxt), inches=FALSE))
    with(x$inUnit, symbols(0, 0, add=TRUE, circles=ringR[1], bg=x$cols[1], fg=NA, inches=FALSE))
    # abline(v=0, h=0, col="lightgray")    # add point of aim
}

drawTarget_BDMP4 <-
function(x, cex=par("cex")) {
    ## background black square
    with(x$inUnit, rect(-extra/2, -extra/2, extra/2, extra/2, col="black", border=NA))
    ## draw circles first, start from the outside
    with(x$inUnit, symbols(rep(0, length(ringR)), rep(0, length(ringR)), add=TRUE,
                           circles=rev(ringR), bg=rev(x$cols), fg=rev(x$colsTxt),
                           inches=FALSE))
    # abline(v=0, h=0, col="lightgray")    # add point of aim
}

#####-----------------------------------------------------------------------
## DSU special targets
drawTarget_DSUa <-
function(x, cex=par("cex")) {
    ## central rectangles and circle
    with(x$inUnit, rect(-ringD8/2,  -ringD8V/2,  ringD8/2,  ringD8V/2,
                        col=x$cols[2], border=x$cols[2]))
    with(x$inUnit, rect(-ringD8/2+ringD8thick, -ringD8V/2+ringD8thick,
                         ringD8/2-ringD8thick,  ringD8V/2-ringD8thick,
                        col=x$cols[1], border=x$cols[1]))
    with(x$inUnit, rect(-ringD9/2,  -ringD9V/2,  ringD9/2,  ringD9V/2,
                        col=x$cols[2], border=x$cols[2]))
    with(x$inUnit, rect(-ringD9/2+ringD9thick, -ringD9V/2+ringD9thick,
                         ringD9/2-ringD9thick,  ringD9V/2-ringD9thick,
                         col=x$cols[1], border=x$cols[1]))
    with(x$inUnit, rect(-ringD10/2, -ringD10V/2, ringD10/2, ringD10V/2,
                        col=x$cols[2], border=x$cols[1]))
    with(x$inUnit, drawCircle(c(0, 0), ringD10i/2, bg=x$cols[1], fg=x$cols[2]))

    ## offset circles
    with(x$inUnit, drawCircle(c(0, ringD9V/2 + ringOff + ringOffD/2),
                              radius=ringOffD/2, bg=x$cols[2], fg=x$cols[2]))
    with(x$inUnit, drawCircle(c(0, ringD9V/2 + ringOff + ringOffD/2),
                              radius=ringOffDi/2, bg=x$cols[1], fg=x$cols[1]))
    with(x$inUnit, rect(-crossW, ringD9V/2 + ringOff + ringOffD/2 - crossThick/2,
                         crossW, ringD9V/2 + ringOff + ringOffD/2 + crossThick/2,
                         col=x$cols[2], border=x$cols[2]))
    with(x$inUnit, rect(-crossThick/2, ringD9V/2 + ringOff + ringOffD/2 - crossW,
                         crossThick/2, ringD9V/2 + ringOff + ringOffD/2 + crossW,
                         col=x$cols[2], border=x$cols[2]))

    ## angles
    ## LU
    with(x$inUnit, rect(-ringD9/2 - angleOff,          ringD9V/2 + angleOff,
                        -ringD9/2 - angleOff + angleL, ringD9V/2 + angleOff - angleThick,
                         col=x$cols[2], border=x$cols[2]))
    with(x$inUnit, rect(-ringD9/2 - angleOff,              ringD9V/2 + angleOff,
                        -ringD9/2 - angleOff + angleThick, ringD9V/2 + angleOff - angleL,
                         col=x$cols[2], border=x$cols[2]))
    ## LD
    with(x$inUnit, rect(-ringD9/2 - angleOff,          -ringD9V/2 - angleOff,
                        -ringD9/2 - angleOff + angleL, -ringD9V/2 - angleOff + angleThick,
                         col=x$cols[2], border=x$cols[2]))
    with(x$inUnit, rect(-ringD9/2 - angleOff,              -ringD9V/2 - angleOff,
                        -ringD9/2 - angleOff + angleThick, -ringD9V/2 - angleOff + angleL,
                         col=x$cols[2], border=x$cols[2]))
    ## RD
    with(x$inUnit, rect(ringD9/2 + angleOff,          -ringD9V/2 - angleOff,
                        ringD9/2 + angleOff - angleL, -ringD9V/2 - angleOff + angleThick,
                        col=x$cols[2], border=x$cols[2]))
    with(x$inUnit, rect(ringD9/2 + angleOff,              -ringD9V/2 - angleOff,
                        ringD9/2 + angleOff - angleThick, -ringD9V/2 - angleOff + angleL,
                        col=x$cols[2], border=x$cols[2]))
    ## RU
    with(x$inUnit, rect(ringD9/2 + angleOff,          ringD9V/2 + angleOff,
                        ringD9/2 + angleOff - angleL, ringD9V/2 + angleOff - angleThick,
                        col=x$cols[2], border=x$cols[2]))
    with(x$inUnit, rect(ringD9/2 + angleOff,              ringD9V/2 + angleOff,
                        ringD9/2 + angleOff - angleThick, ringD9V/2 + angleOff - angleL,
                        col=x$cols[2], border=x$cols[2]))
}

drawTarget_DSUb <-
function(x, cex=par("cex")) {
    ## draw circle sectors first, start from the outside
    sapply(rev(seq_along(x$inUnit$ringR)), function(i) {
        ang <- (180/pi)*with(x$inUnit, atan((ringRV-ringR) / ringR))
        with(x$inUnit, drawDSUOval(c(0, 0),
                            h=ringRV[i]-ringR[i],
                            radius=ringR[i],
                            angle=ang[i],
                            fg=x$colsTxt[i],
                            bg=x$cols[i]))
    })

    rings1 <- with(x, seq(from=maxVal-nRings+1, to=maxVal-1, length.out=nRings-1)) # left side of center
    rings2 <- c(rings1, rev(rings1))                         # both sides
    posX1  <- with(x$inUnit, ringR[3:length(ringR)] - (ringW/2)) # x: right side of center
    posX2  <- c(-rev(posX1), posX1)                          # both sides
    posY1  <- with(x$inUnit, ringRV[3:length(ringRV)] - (ringWV/2)) # y: right side of center
    posY2  <- c(-rev(posY1), posY1)                          # y: both sides
    cols1  <- with(x$inUnit, x$colsTxt[3:length(ringR)])     # right side of center
    cols2  <- c(rev(cols1), cols1)                           # both sides

    text(posX2, 0, cex=cex, label=rings2, col=cols2)
    text(0, posY2, cex=cex, label=rings2, col=cols2)
}
