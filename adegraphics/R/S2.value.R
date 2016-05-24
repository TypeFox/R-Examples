#########################################################
###                     s.value                        ##
#########################################################

## TO DO: calcul place legend, taille des points
## Remarque ==> pour size, si couleur selon <0 ou >0 il faut s'assurer que 0 ne sera donc pas dans un intervalle? (inclus ex [-1, 1])
setClass(
    Class = "S2.value",
    contains = "ADEg.S2",
    )


setMethod(
    f = "initialize",
    signature = "S2.value",
    definition = function(.Object, data = list(dfxy = NULL, z = NULL, xax = 1, yax = 2, frame = 0, storeData = TRUE), ...) {
        .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.S2 initialize
        .Object@data$z <- data$z
        return(.Object)
    })


setMethod(
    f = "prepare",
    signature = "S2.value",
    definition = function(object) {
        name_obj <- deparse(substitute(object))
        
        ## pre-management of graphics parameters      
        oldparamadeg <- adegpar()
        on.exit(adegpar(oldparamadeg))
        adegtot <- adegpar(object@adeg.par)
        
        if(object@data$storeData)
            z <- object@data$z
        else
            z <- eval(object@data$z, envir = sys.frame(object@data$frame))
        
        if(is.null(object@adeg.par$ppoints$alpha))
            adegtot$ppoints$alpha <- 0.9
        if(is.null(object@adeg.par$ppoints$cex))
            adegtot$ppoints$cex <- 1
        if(is.null(object@adeg.par$porigin$include) & (any(names(object@g.args) %in% c("Sp", "nbobject"))))
            adegtot$porigin$include <- FALSE
        
        if(is.null(object@g.args$breaks))
            object@s.misc$breaks.update <- pretty(z, object@g.args$nclass)
        else
            object@s.misc$breaks.update <- object@g.args$breaks 
        
        object@s.misc$breaks.update <- breakstest(object@s.misc$breaks.update, z, n = length(object@s.misc$breaks.update))
        n <- length(object@s.misc$breaks.update)
        
        ## symbols for z = center
        if(!is.null(object@g.args$centerpar)) {
            default <- list(pch = 4, cex = 1, col = "black")
            if(is.list(object@g.args$centerpar))
                object@g.args$centerpar <- modifyList(default, object@g.args$centerpar, keep.null = TRUE)
            else
                object@g.args$centerpar <- default
        }
        
        if(is.null(object@adeg.par$psub$position))
            adegtot$psub$position <- "topleft"
        
        if(!is.null(object@g.args$col)) {
            switch(object@g.args$method,
                   size = {
                       if(length(object@g.args$col) != 2)
                           stop("if method size choosen, col vector should be size 2", call. = FALSE)
                       adegtot$ppoints$col <- object@g.args$col ## color given by the user
                   },
                   color = {
                       if(length(object@g.args$col) < (n - 1))
                           stop(paste("not enough colors defined for method color, at least ", (n - 1), " colors expected", sep = ""), call. = FALSE)
                       adegtot$ppoints$col <- object@g.args$col[1:(n - 1)]  ## color given by the user
                   })
        } else {
            if(object@g.args$method == "color")
                adegtot$ppoints$col <- adegtot$ppalette$quanti(n - 1)
            else
                adegtot$ppoints$col <- adegtot$ppalette$quanti(2)
        }
        
        ## object modification before calling inherited method
        object@adeg.par <- adegtot
        callNextMethod() ## prepare graph
         
        assign(name_obj, object, envir = parent.frame())
    })


## Draw symbols according to the different methods
setMethod(
    f = "panel",
    signature = "S2.value",
    definition = function(object, x, y) {
        if(object@data$storeData)
            zorig <- object@data$z
        else
            zorig <- eval(object@data$z, envir = sys.frame(object@data$frame))
        
        adeg.panel.values(x = x, y = y, z = zorig, method = object@g.args$method, symbol = object@g.args$symbol, ppoints = object@adeg.par$ppoints, 
                          breaks = object@s.misc$breaks.update, centerpar = object@g.args$centerpar, center = object@g.args$center)
    })


s.value <- function(dfxy, z, breaks = NULL, xax = 1, yax = 2, method = c("size", "color"), symbol = c("square", "circle", "diamond", "uptriangle", "downtriangle"),
                    col = NULL, nclass = 4, center = 0, centerpar = NULL, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
    
    ## evaluation of some parameters
    thecall <- .expand.call(match.call())
    thecall$method <- match.arg(method)
    if(thecall$method == "color") {
        if(center != 0 | !is.null(centerpar))
            warning("'center' and 'centerpar' are not used with 'color' method", call. = FALSE)      
        center <- 0
        centerpar <- NULL
    }
    thecall$center <- center
    thecall$centerpar <- centerpar

    thecall$symbol <- match.arg(symbol)
    df <- try(as.data.frame(eval(thecall$dfxy, envir = sys.frame(sys.nframe() + pos))), silent = TRUE)
    z <- eval(thecall$z, envir = sys.frame(sys.nframe() + pos))
    if((class(df) == "try-error") | is.null(thecall$dfxy)) ## non convenient dfxy argument
        stop("non convenient selection for dfxy (can not be converted to dataframe)", call. = FALSE)
    if(NROW(df) != NROW(z))
        stop("dfxy and z should have the same number of rows", call. = FALSE)
    
    ## parameters sorted
    sortparameters <- sortparamADEg(...)
    
    ## facets
    if(!is.null(facets)) {
        if((length(xax) == 1 & length(yax) == 1) & NCOL(z) == 1)
            object <- multi.facets.S2(thecall, sortparameters$adepar, samelimits = sortparameters$g.args$samelimits)
        else 
            stop("Facets are not allowed with multiple xax/yax or multiple z", call. = FALSE)
    }
    
    ## multiple axes
    else if((length(xax) > 1 | length(yax) > 1)) {
        if(NCOL(z) == 1)
            object <- multi.ax.S2(thecall)
        else 
            stop("Multiple xax/yax are not allowed with multiple z", call. = FALSE)
    }
    
    ## multiple z
    else if(NCOL(z) > 1) {
        object <- multi.variables.S2(thecall, "z")
    }
    
    ## simple ADEg graphic
    else {
        if(length(sortparameters$rest))
            warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
        
        ## creation of the ADEg object
        g.args <- c(sortparameters$g.args, list(method = thecall$method, symbol = thecall$symbol, center = center, breaks = breaks, col = col,
                                                nclass = nclass, centerpar = centerpar))
        if(storeData)
            tmp_data <- list(dfxy = dfxy, xax = xax, yax = yax, z = z, frame = sys.nframe() + pos, storeData = storeData)
        else
            tmp_data <- list(dfxy = thecall$dfxy, xax = xax, yax = yax, z = thecall$z, frame = sys.nframe() + pos, storeData = storeData)
        object <- new(Class = "S2.value", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = as.call(thecall))
        
        ## preparation of the graph
        prepare(object)
        setlatticecall(object)
        if(add)
            object <- add.ADEg(object)
    }
    
    if(! add & plot)
        print(object)
    invisible(object)
}
