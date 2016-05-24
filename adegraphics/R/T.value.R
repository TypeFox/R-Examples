setClass(
    Class = "T.value",
    contains = "ADEg.T"
    )


setMethod(
    f = "prepare",
    signature = "T.value",
    definition = function(object) {
        name_obj <- deparse(substitute(object))
        
        ## pre-management of graphics parameters      
        oldparamadeg <- adegpar()
        on.exit(adegpar(oldparamadeg))
        adegtot <- adegpar(object@adeg.par)
        
        if(object@data$storeData) {
            z <- as.vector(as.matrix(object@data$dftab))
            dftab <- object@data$dftab
            labelsx <- object@data$labelsx
            labelsy <- object@data$labelsy
        } else {
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
        
        if(is.null(object@adeg.par$ppoints$cex))
            adegtot$ppoints$cex <- 1
        
        if(is.null(object@adeg.par$ppoints$alpha))
            adegtot$ppoints$alpha <- 1
        
        if(is.null(labelsx))
            adegtot$ptable$x$tck <- 0
        if(is.null(labelsy))
            adegtot$ptable$y$tck <- 0
        
        ## symbols for z = center
        if(!is.null(object@g.args$centerpar)) {
            default <- list(pch = 4, cex = 1, col = "black")
            if(is.list(object@g.args$centerpar))
                object@g.args$centerpar <- modifyList(default, object@g.args$centerpar, keep.null = TRUE)
            else
                object@g.args$centerpar <- default
        }
        
        if(!is.null(object@g.args$col)) {
            switch(object@g.args$method,
                   size = {
                       if(length(object@g.args$col) != 2 & !inherits(dftab, "table") & !inherits(dftab, "dist"))
                           stop("if method size choosen, col vector should be size 2", call. = FALSE)
                       adegtot$ppoints$col <- object@g.args$col ## color given by the user
                   },
                   color = {
                       if(length(object@g.args$col) < (n - 1))
                           stop(paste("not enough colors defined for method color, at least ", (n - 1), " colors expected", sep = "") , call. = FALSE)
                       adegtot$ppoints$col <- object@g.args$col[1:(n - 1)] ## color given by the user
                   })
        } else {
            if(object@g.args$method == "color")
                adegtot$ppoints$col <- adegtot$ppalette$quanti(n - 1)
            else if(inherits(dftab, "table") | inherits(dftab, "dist")) {
                adegtot$ppoints$col <- adegtot$ppalette$quanti(2)
            }
            else
                adegtot$ppoints$col <- adegtot$ppalette$quanti(2)
        }
        
        ## object modification before calling inherited method
        object@adeg.par <- adegtot
        callNextMethod() ## prepare graph
        
        assign(name_obj, object, envir = parent.frame())
    })


setMethod(
    f = "panel",
    signature = "T.value",
    definition = function(object, x, y) {
        if(object@data$storeData)
            dftab <- as.matrix(object@data$dftab)
        else
            dftab <- as.matrix(eval(object@data$dftab, envir = sys.frame(object@data$frame)))
        adeg.panel.values(x = x[col(dftab)], y = y[row(dftab)], z = as.vector(dftab), center = object@g.args$center, method = object@g.args$method,
                          symbol = object@g.args$symbol, ppoints = object@adeg.par$ppoints, breaks = object@s.misc$breaks.update, centerpar = object@g.args$centerpar)
    })


table.value <- function(dftab, coordsx = 1:ncol(as.matrix(dftab)), coordsy = nrow(as.matrix(dftab)):1, labelsx, labelsy, breaks = NULL, method = c("size", "color"),
                        symbol = c("square", "circle", "diamond", "uptriangle", "downtriangle"), col = NULL, nclass = 3, center = 0, centerpar = NULL, plot = TRUE,
                        storeData = TRUE, add = FALSE, pos = -1, ...) {
    
    ## 4 different types can be used as tab :
    ## distance matrix (dist), contingency table (table), data.frame or matrix
    
    ## evaluation of some parameters
    thecall <- .expand.call(match.call())
    thecall$method <- match.arg(method)
    thecall$symbol <- match.arg(symbol)
    dftab <- eval(thecall$dftab, envir = sys.frame(sys.nframe() + pos))
    if(any(is.na(dftab)))
        stop("NA entries not accepted")
    
    if(inherits(dftab, "dist")) {
        if(missing(labelsx)){
            thecall$labelsx <- labelsx <- NULL
            if(!is.null(attr(dftab, "Labels")))
                if(storeData)
                    labelsx <- attr(dftab, "Labels")
    		else
                    thecall$labelsx <- call("attr", thecall$dftab, "Labels")
        } 
        
        if(missing(labelsy)) {
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
        if(missing(labelsy)) {
            thecall$labelsy <- labelsy <- NULL
            if(!is.null(rownames(dftab)))
                if(storeData)
                    labelsy <- rownames(dftab)
                else
                    thecall$labelsy <- call("rownames", thecall$dftab)
        }
        
        if(missing(labelsx)) {
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
    g.args <- c(sortparameters$g.args, list(breaks = breaks, method = thecall$method, symbol = thecall$symbol, center = thecall$center, col = col, nclass = nclass, centerpar = centerpar))
    if(storeData)
        tmp_data <- list(dftab = dftab, coordsx = coordsx, coordsy = coordsy, labelsx = labelsx, labelsy = labelsy, frame = sys.nframe() + pos, storeData = storeData)
    else
        tmp_data <- list(dftab = thecall$dftab, coordsx = thecall$coordsx, coordsy = thecall$coordsy, labelsx = thecall$labelsx, labelsy = thecall$labelsy, frame = sys.nframe() + pos, storeData = storeData)
    
    if(inherits(dftab, "table")) {
        condres <- pmatch(c("ablineX", "ablineY", "meanX", "meanY"), names(sortparameters$rest))
        if(any(!is.na(condres))) {
            tmplist <- sortparameters$rest[condres[!is.na(condres)]]
            names(tmplist) <- c("ablineX", "ablineY", "meanX", "meanY")[which(!is.na(condres))]
            sortparameters$rest <- sortparameters$rest[-condres[(!is.na(condres))]]
            g.args <- c(g.args, tmplist)
        }
        g.args[c("ablineX", "ablineY", "meanX", "meanY")[which(is.na(condres))]] <- FALSE
        object <- new(Class = "T.cont", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())
    } else
        object <- new(Class = "T.value", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())

    if(length(sortparameters$rest))
        warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
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

