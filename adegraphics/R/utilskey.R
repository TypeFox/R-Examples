setMethod(
    f = "addkey",
    signature = "ADEg",
    definition = function(object) {
        object@adeg.par$plegend$drawKey <- TRUE
        object@g.args$key <- createkey(object)
        return(object)
    })

setMethod(
    f = "createkey",
    signature = "ADEg",
    definition = function(object) {
        if(object@adeg.par$plegend$drawKey){
            res <- object@g.args$key
        }
        else
            res <- NULL
        return(res)
    })

setMethod(
    f = "createkey",
    signature = "S2.value",
    definition = function(object) {
         return(.createkeyvalue(object, type = "S2"))
    })

setMethod(
    f = "createkey",
    signature = "T.value",
    definition = function(object) {
        return(.createkeyvalue(object, type = "T"))
    })

setMethod(
    f = "createkey",
    signature = "S2.class",
    definition = function(object) {
         return(.createkeyclass(object, type = "S2"))
    })

setMethod(
    f = "createkey",
    signature = "S1.class",
    definition = function(object) {
         return(.createkeyclass(object, type = "S1"))
    })

setMethod(
    f = "createkey",
    signature = "Tr.class",
    definition = function(object) {
         return(.createkeyclass(object, type = "Tr"))
    })

setMethod(
    f = "createcolorkey",
    signature = "ADEg",
    definition = function(object) {
        if(object@adeg.par$plegend$drawColorKey){
            res <- object@g.args$legend
        }
        else
            res <- NULL
        return(res)
    })

setMethod(
    f = "createcolorkey",
    signature = "T.image",
    definition = function(object) {
        ## add a small space before the colorkey
        trellis.par.set(layout.widths = list(axis.key.padding = 1))
        return(.createcolorkeyimage(object, type = "T"))
    })

setMethod(
    f = "createcolorkey",
    signature = "S2.image",
    definition = function(object) {
        ## add a small space before the colorkey
        trellis.par.set(layout.widths = list(axis.key.padding = 1))
        return(.createcolorkeyimage(object, type = "S2"))
    })


.createkeyvalue <- function(object, type = c("T", "S2")) {
    type <- match.arg(type)
    res <- NULL
    if(object@adeg.par$plegend$drawKey){
        res <- list()
        res$points$pch <- .symbol2pch(object@g.args$symbol)
        cstnormal <- 5 ## same value in adeg.panel.value
        if(object@g.args$method == "size"){
            center <- object@g.args$center
            breaks <- unique(c(object@s.misc$breaks.update, signif(center, 5)))
            maxsize <- max(abs(breaks))
            breaks <- breaks[order(breaks, decreasing = FALSE)]
            l0 <- length(breaks)
            breaks <- (breaks[1:(l0 - 1)] + breaks[2:l0]) / 2
            res$text$lab <- as.character(breaks)
            size <- breaks - center
            res$points$cex <- .proportional_map(size, maxsize) * object@adeg.par$ppoints$cex[1] * cstnormal
            res$points$fill <- object@adeg.par$ppoints$col[ifelse(breaks < center, 1, 2)]
            res$points$col <- object@adeg.par$ppoints$col[ifelse(breaks < center, 2, 1)]
            
        } else if(object@g.args$method == "color"){
            breaks <- object@s.misc$breaks.update
            l0 <- length(breaks)
            res$points$cex <- object@adeg.par$ppoints$cex[1] * cstnormal / 2 * object@adeg.par$plegend$size
            res$text$lab <- paste("[", breaks[l0], ";", breaks[l0 - 1], "]", sep = "")
            for(i in (l0 - 1):2)
                res$text$lab <- c(res$text$lab, paste("]", breaks[i], ";", breaks[i - 1], "]", sep = ""))
            res$points$fill <- object@adeg.par$ppoints$col[1:length(res$text$lab)]
            res$points$col <- "black"
        }
        res$columns <- length(res$text$lab)
        res$border <- TRUE
        res$between <- 0.1 * object@adeg.par$plegend$size
        res$between.columns <- 0 * object@adeg.par$plegend$size
        res$padding.text <- 1.2 * max(res$points$cex)
        res$text$cex <- object@adeg.par$plegend$size
        if(is.null(object@g.args$key$space)){
            if(type == "T"){
                res$x <- 0
                res$y <- 0
            } else {
                res$corner <- c(0,0)
                res$x <- 0.01
                res$y <- 0.01
            }
        }
        res$background <- object@adeg.par$pbackground$col
        res <- modifyList(res, as.list(object@g.args$key), keep.null = TRUE)
    }
    return(res)
}


.createkeyclass <- function(object, type = c("S1", "S2", "Tr")) {
    type <- match.arg(type)
    res <- NULL
  
    if(object@adeg.par$plegend$drawKey){
        res <- list()
        if(object@data$storeData) 
            fac <- as.factor(object@data$fac)
        else 
            fac <- as.factor(eval(object@data$fac, envir = sys.frame(object@data$frame)))

        res$text$lab <- levels(fac)
        res$text$col <- object@adeg.par$plabels$col

        if(object@adeg.par$ppoints$cex > 0){
            res$points$pch <- object@adeg.par$ppoints$pch
            res$points$col <- object@adeg.par$ppoints$col
            res$points$fill <- object@adeg.par$ppoints$fill
            
        } else if(!is.null(object@g.args$chullSize)){
            if(object@g.args$chullSize > 0){
                res$rectangles$border <- object@adeg.par$ppolygons$border
                res$rectangles$col <- object@adeg.par$ppolygons$col
                res$rectangles$alpha <- object@adeg.par$ppolygons$alpha
            }
            
        } else if(object@g.args$ellipseSize > 0){
            res$rectangles$border <- object@adeg.par$pellipses$border
            res$rectangles$col <- object@adeg.par$pellipses$col
            res$rectangles$alpha <- object@adeg.par$pellipses$alpha
            
        } else if(object@g.args$starSize > 0){
            res$lines$col <- object@adeg.par$plines$col
            res$lines$lty <- object@adeg.par$plines$lty
            res$lines$lwd <- object@adeg.par$plines$lwd
        }
        
        res$between <- 0.1 * object@adeg.par$plegend$size
        res$between.columns <- 0 * object@adeg.par$plegend$size
        res$text$cex <- object@adeg.par$plegend$size
        
        if(is.null(object@g.args$key$space)){
            if(type == "S2"){
                res$corner <- c(0,0)
                res$x <- 0.01
                res$y <- 0.01
            }
            
            if(type == "S1"){
                res$corner <- c(0,1)
                res$x <- 0.01
                res$y <- 0.99
            }
            
            if(type == "Tr"){
                res$corner <- c(1,1)
                res$x <- 0.99
                res$y <- 0.99
            }
        }
        
        res$background <- object@adeg.par$pbackground$col
    
        res <- modifyList(res, as.list(object@g.args$key), keep.null = TRUE)
    }   
    return(res)
}


.createcolorkeyimage <- function(object, type = c("T", "S2")) {
    type <- match.arg(type)
    res <- NULL
    if(object@adeg.par$plegend$drawColorKey){
        res <- list(right = list(fun = draw.colorkey, args = list(key = list(col = object@adeg.par$ppoints$col, at = object@s.misc$breaks.update))))
    }
    return(res)
}
