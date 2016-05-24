setClass(
  Class = "T.image",
  contains = "ADEg.T"
)


setMethod(
  f = "prepare",
  signature = "T.image",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)

    if(object@data$storeData) {
      coordsx <- object@data$coordsx
      coordsy <- object@data$coordsy
      z <- as.vector(as.matrix(object@data$dftab))
      dftab <- object@data$dftab
      labelsx <- object@data$labelsx
      labelsy <- object@data$labelsy
     } else {
      coordsx <- eval(object@data$coordsx, envir = sys.frame(object@data$frame))
      coordsy <- eval(object@data$coordsy, envir = sys.frame(object@data$frame))
      z <- as.vector(as.matrix(eval(object@data$dftab, envir = sys.frame(object@data$frame))))
      dftab <- eval(object@data$dftab, envir = sys.frame(object@data$frame))
      labelsx <- eval(object@data$labelsx, envir = sys.frame(object@data$frame))
      labelsy <- eval(object@data$labelsy, envir = sys.frame(object@data$frame))
    }
    
    if(is.null(object@g.args$breaks))
      object@s.misc$breaks.update <- pretty(z, object@g.args$nclass)
    else
      object@s.misc$breaks.update <- object@g.args$breaks
    
    object@s.misc$breaks.update <- breakstest(object@s.misc$breaks.update, z, n = length(object@s.misc$breaks.update))
    n <- length(object@s.misc$breaks.update)
    
    if(!is.null(object@g.args$col)) {
      if(length(object@g.args$col) < (n - 1))
        stop(paste("not enough colors defined, at least ", (n - 1), " colors expected", sep = ""), call. = FALSE)
      adegtot$ppoints$col <- object@g.args$col[1:(n - 1)]  ## color given by the user
  	} else {
      if(is.null(object@adeg.par$ppoints$col))
        adegtot$ppoints$col <- adegtot$ppalette$quanti(n - 1)
    }
  
    ## inspired by level.colors from lattice

    if(adegtot$plegend$drawColorKey)
      adegtot$ptable$y$pos <- "left"
    if(is.null(object@adeg.par$pgrid$col))
      adegtot$pgrid$col <- "black"
    if(is.null(object@adeg.par$pgrid$lwd))
      adegtot$pgrid$lwd <- 0.6
    if(is.null(object@adeg.par$pgrid$draw))
      adegtot$pgrid$draw <- FALSE ## no cells border by default
    
    if(is.null(labelsx))
        adegtot$ptable$x$tck <- 0
    if(is.null(labelsy))
        adegtot$ptable$y$tck <- 0
        
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    ## TODO:  improve the code to avoid some repetition with the parent function
    wx <- range(coordsx)
    dx <- (diff(wx) + 1) / length(coordsx)
    wy <- range(coordsy)
    dy <- (diff(wy) + 1) / length(coordsy)
    ## add an half cell at both sides
    object@g.args$xlim <- wx + c(-0.5, 0.5) * dx
    object@g.args$ylim <- wy + c(-0.5, 0.5) * dy 
    assign(name_obj, object, envir = parent.frame())
  })

  
setMethod(
  f = "panel",
  signature = "T.image",
  definition = function(object, x, y) {
    ## x is data$coordsx and y is data$coordsy
    if(object@data$storeData)
      dftab <- as.matrix(object@data$dftab)
    else
      dftab <- as.matrix(eval(object@data$dftab, envir = sys.frame(object@data$frame)))
    
    xx <- x[!is.na(x)]
    yy <- y[!is.na(y)]
    zz <- as.vector(dftab)
    dx <- diff(sort(xx)) / 2
    dy <- diff(sort(yy)) / 2
    dx <- c(dx[1], dx)
    dy <- c(dy[1], dy)

    ## draw values
    panel.levelplot.raster(x = xx[col(dftab)], y = yy[row(dftab)], z = zz, subscripts = TRUE, col.regions = object@adeg.par$ppoints$col,
                           at = object@s.misc$breaks.update, contour = FALSE, region = TRUE)

    ## draw grid (cells border)
    if(object@adeg.par$pgrid$draw) {
      xbis <- c(min(xx) - dx[1], xx + dx, max(xx) + dx[length(dx)])
      ybis <- c(min(yy) - dy[1], yy + dy, max(yy) + dy[length(dy)])
      panel.abline(h = ybis, v = xbis, col = object@adeg.par$pgrid$col, lwd = object@adeg.par$pgrid$lwd, lty = object@adeg.par$pgrid$lty)
    }
  })


## TODO: decider quelle classe on prend en compte
## a faire: verifier espacement correct de coordsx et coordsy
## que faire de la sous grille?
## attention, coordsx et coordsy ne serve qu'a donner l'ordre de trace, ils seront considere comme egalement espace, sinon fonction a revoir
table.image <- function(dftab, coordsx = 1:ncol(as.matrix(dftab)), coordsy = nrow(as.matrix(dftab)):1, labelsx = NULL, labelsy = NULL, nclass = 3, breaks = NULL, col = NULL, plot = TRUE, 
                        storeData = TRUE, add = FALSE, pos = -1, ...) {
    
    ## 4 different types can be used as tab :
    ## distance matrix (dist), contingency table (table), data.frame or matrix
    thecall <- .expand.call(match.call())
    dftab <- eval(thecall$dftab, envir = sys.frame(sys.nframe() + pos))
    
    ## modify coordsx/coordsy positions (we use only the order not the values)
    thecall$coordsx <- call("rank", thecall$coordsx, ties.method = "first")
    thecall$coordsy <- call("rank", thecall$coordsy, ties.method = "first")
    
    if(inherits(dftab, "dist")) {
        if(missing(labelsx)){
            thecall$labelsx <- labelsx <- NULL
            if(!is.null(attr(dftab, "Labels")))
                if(storeData)
                    labelsx <- attr(dftab, "Labels")
                else
                    thecall$labelsx <- call("attr", thecall$dftab, "Labels")
        }
        if(missing(labelsy)){
            thecall$labelsy <- labelsy <- NULL
            if(!is.null(attr(dftab, "Labels")))
                if(storeData)
                    labelsy <- attr(dftab, "Labels")
                else
                    thecall$labelsy <- call("attr", thecall$dftab, "Labels")
        }
        ## coordsx and coordsy should be identical for dist objects (symmetric)
        thecall$coordsx <- call(":", 1, call("attr", thecall$dftab, "Size"))
        thecall$coordsy <- call(":", call("attr", thecall$dftab, "Size"), 1)
        
    } else { ## data.frame, matrix, table
        if(missing(labelsy)){
            thecall$labelsy <- labelsy <- NULL
            if(!is.null(rownames(dftab)))
                if(storeData)
                    labelsy <- rownames(dftab)
                else
                    thecall$labelsy <- call("rownames", thecall$dftab)
        }
        
        if(missing(labelsx)){
            thecall$labelsx <- labelsx <- NULL
            if(!is.null(colnames(dftab)))
                if(storeData)
                    labelsx <- colnames(dftab)
                else
                    thecall$labelsx <- call("colnames", thecall$dftab)
        }
    }
  
    ## parameters sorted
    sortparameters <- sortparamADEg(...)
    
    ## creation of the ADEg object
    if(length(sortparameters$rest))
        warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    g.args <- c(sortparameters$g.args, list(breaks = breaks, nclass = nclass, col = col))
    if(storeData)
  	tmp_data <- list(dftab = dftab, coordsx = coordsx, coordsy = coordsy, labelsx = labelsx, labelsy = labelsy, frame = sys.nframe() + pos, storeData = storeData)
    else
        tmp_data <- list(dftab = thecall$dftab, coordsx = thecall$coordsx, coordsy = thecall$coordsy, labelsx = thecall$labelsx, labelsy = thecall$labelsy, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "T.image", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())
    
    ## preparation of the graph
    prepare(object)
    setlatticecall(object)
    if(add)
        object <- add.ADEg(object)
    else
        if(plot)
            print(object)
    invisible(object) 
}
