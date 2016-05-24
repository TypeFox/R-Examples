
.gridGraphicsEnv <- new.env()

init <- function(dl, width=NULL, height=NULL) {
    if (dev.cur() == 1) {
        dev.new()
    }
    # The graphics device we will draw onto
    assign("pd", dev.cur(), .gridGraphicsEnv)
    din <- par("din")
    if (is.null(width))
        width <- din[1]
    if (is.null(height))
        height <- din[2]
    # An off-screen graphics device 
    pdf(NULL, width=width, height=height)
    assign("rd", dev.cur(), .gridGraphicsEnv)
    # NULL out the saved display list then replay it in order
    # to restore basic graphics settings (like device background)
    dlStub <- dl
    dlStub[1] <- list(NULL)
    replayPlot(dlStub)
    # Go back to device we will draw onto
    dev.set(playDev())
    initInnerAlpha()
    initFigureAlpha()
    initPlotIndex()
    initPlotAlpha()
    initWindowIndex()
    initWindowAlpha()
    initWindowPlotAlpha()
    initClip()
    # Remove any grobname indices
    rm(list=ls(envir=.gridGraphicsEnv, pattern=".gridGraphicsIndex$"),
       envir=.gridGraphicsEnv)
}

shutdown <- function() {
    # Close the off-screen graphics device
    dev.set(recordDev())
    dev.off()
    # Make sure we go back to the device we are replaying onto
    dev.set(playDev())
    invisible()
}

initClip <- function() {
    assign("currentClip", NULL, .gridGraphicsEnv)
}

setClip <- function(x, y, w, h) {
    assign("currentClip", c(x, y, w, h), .gridGraphicsEnv)
}

getClip <- function() {
    get("currentClip", .gridGraphicsEnv)
}

playDev <- function() {
    get("pd", .gridGraphicsEnv)
}
    
recordDev <- function() {
    get("rd", .gridGraphicsEnv)
}

indexFuns <- function() {
    index <- 0
    init <- function() {
        index <<- 0
    }
    increment <- function() {
        index <<- index + 1
    }
    get <- function() {
        index
    }
    list(init=init, increment=increment, get=get)
}

pif <- indexFuns()
initPlotIndex <- pif$init
incrementPlotIndex <- pif$increment
plotIndex <- pif$get

wif <- indexFuns()
initWindowIndex <- wif$init
incrementWindowIndex <- wif$increment
windowIndex <- wif$get

alphaIndexFuns <- function() {
    index <- 0
    init <- function() {
        index <<- 0
    }
    increment <- function() {
        index <<- index + 1
    }
    get <- function() {
        if (index == 0) {
            ""
        } else {
            paste(rep(LETTERS[(index - 1) %% 26 + 1], (index - 1) %/% 26 + 1),
                  collapse="")
        }
    }
    set <- function(x) {
        if (x == "") {
            index <<- 0
        } else {
            n <- nchar(x)
            chars <- rev(strsplit(x, "")[[1]])
            index <<- sum(sapply(chars, function(y) which(LETTERS == y)))
        }
    }
    list(init=init, increment=increment, get=get, set=set)
}

iiaf <- alphaIndexFuns()
initInnerAlpha <- iiaf$init
incrementInnerAlpha <- iiaf$increment
innerAlpha <- iiaf$get

fiaf <- alphaIndexFuns()
initFigureAlpha <- fiaf$init
incrementFigureAlpha <- fiaf$increment
figureAlpha <- fiaf$get

piaf <- alphaIndexFuns()
initPlotAlpha <- piaf$init
incrementPlotAlpha <- piaf$increment
plotAlpha <- piaf$get

wiaf <- alphaIndexFuns()
initWindowAlpha <- wiaf$init
incrementWindowAlpha <- wiaf$increment
windowAlpha <- wiaf$get

wpiaf <- alphaIndexFuns()
initWindowPlotAlpha <- wpiaf$init
windowPlotAlpha <- wpiaf$get
setWindowPlotAlpha <- wpiaf$set

prefixFuns <- function() {
    prefix <- "graphics"
    get <- function() {
        prefix
    }
    set <- function(x) {
        prefix <<- as.character(x)
    }
    list(get=get, set=set)
}
pf <- prefixFuns()

prefix <- pf$get
setPrefix <- pf$set

grobname <- function(label, unique=FALSE) {
    if (unique) {
        paste(prefix(), label, sep="-")
    } else {
        stub <- paste(prefix(), "plot", plotIndex(), label, sep="-")
        indexName <- paste0(stub, ".gridGraphicsIndex")
        if (exists(indexName, .gridGraphicsEnv)) {
            index <- get(indexName, .gridGraphicsEnv)
        } else {
            index <- 1
        }
        assign(indexName, index + 1, .gridGraphicsEnv)
        paste(stub, index, sep="-")
    }
}

vpname <- function(type, clip=FALSE) {
    switch(type,
           root=paste(prefix(), type, sep="-"),
           inner=paste(prefix(), paste0(type, innerAlpha()), sep="-"),
           figure={
               if (clip) {
                   paste(prefix(), type, paste0(plotIndex(), figureAlpha()),
                         "clip", sep="-")
               } else {
                   paste(prefix(), type, paste0(plotIndex(), figureAlpha()),
                         sep="-")
               }
           },
           plot={
               if (clip) {
                   paste(prefix(), type, paste0(plotIndex(), plotAlpha()),
                         "clip", sep="-")
               } else {
                   paste(prefix(), type, paste0(plotIndex(), plotAlpha()),
                         sep="-")
               }
           },
           window=paste(prefix(), type, plotIndex(),
                         paste0(windowIndex(), windowAlpha()), sep="-"),
           # NOTE that the "window" "plot" vp uses a potentially different
           #      alpha than the "plot" vp (see clip.R for why)
           windowplot={
               if (clip) {
                   paste(prefix(), "plot",
                         paste0(plotIndex(), windowPlotAlpha()),
                         "clip", sep="-")
               } else {
                   paste(prefix(), "plot",
                         paste0(plotIndex(), windowPlotAlpha()),
                         sep="-")
               }
           })
}
