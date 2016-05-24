##*****************************************************************************
## AUTHOR: Yves Deville
##
## CONTAINS
##
## o The 'plot' and 'lines' methods for the "Renouv" class.
##
## o Misc utility functions to build legends on Return levels plot
##   in a semi-automatic fashion.
##
## DO NOT ROXYGENISE!!!
##
##*****************************************************************************

##' Make translucient colors.
##'
##' @title Make translucient colors
##'
##' @param colors A vector of colors in a format that can be understood by
##' \code{\link{col2rgb}}.
##'
##' @param alpha Level of opacity ("0" means fully transparent and
##' "max" means opaque).  After recycling to reach the required
##' length, this value or vector is used as \code{alpha} in
##' \code{\link{rgb}}.
##'
##' @return A vector of translucient (or semi-transparent) colors.
##'
translude <- function (colors, alpha = 0.6) {
    L <- pmax(length(colors), length(alpha))
    colors <- rep(colors, length.out = L)
    alpha <- rep(alpha, length.out = L)
    rgb <- as.matrix(col2rgb(colors)/255)
    colors2 <- rgb(red = rgb["red", ], green = rgb["green", ], 
                   blue = rgb["blue", ], alpha = alpha)
}


`RLpar0` <- function(mono = TRUE) {
    ##=======================================================================
    ## Set default graphical parameters for Return Level plots
    ##========================================================================
    if (mono) {
        l.cols <- c("black", "black")
        p.cols <- "black"
        l.typs <- c("solid", "dashed", "dotted")
    } else {
        l.cols <- c("SteelBlue4", "orangered", "SpringGreen3", "purple",
                    "firebrick")
        p.cols <- c("black", "orange", "SpringGreen1", "mistyrose", "firebrick1")
        l.typs <- c("solid", "solid", "solid")
    }
    block.cols <- c("orangered", "SpringGreen3", "purple", "DarkCyan", "SteelBlue3",
                    "GoldenRod", "ForestGreen", "firebrick", "darkOliveGreen3",
                    "SlateBlue")
    block.bg <- c("yellow", "SpringGreen1", "mistyrose", "cyan", "SteelBlue1",
                  "gold", "ForestGreen", "tomato", "darkOliveGreen1", "lavender")
    block.pch <- c(21:25, 22) 
    
    l.cols <- rep(l.cols, length.out = 10L)
    p.cols <- rep(p.cols, length.out = 10L)
    l.typs <- rep(l.typs, length.out = 10L)
    block.cols <- rep(block.cols, length.out = 10L)
    block.bg <- rep(block.bg, length.out = 10L)
    block.pch <- rep(block.pch, length.out = 10L)
    .RLpar <- list()
    ## quantile (return level) curve
    .RLpar[["quant"]] <- list(type = "l", col = "black", lwd = 2, lty = "solid")
    ## sample data
    .RLpar[["OT"]] <- list(col = "black", pch = 16, cex = 0.8, bg = "black")
    ## confidence levels
    prov <- list()
    nconf <- 6L
    for (i in 1L:nconf) {
        prov[[paste("conf", i, sep = "")]] <-
            list(lty = i+1L, col = l.cols[i], lwd = 2)
    }
    .RLpar[["conf"]] <- prov
    ## historical MAX data
    prov <- list()
    nMAX <- 10L
    for (i in 1L:nMAX) {
        prov[[paste("block", i, sep = "")]] <-
            list(col = block.cols[i], pch = block.pch[i], cex = 1.1, lwd = 2,
                 bg = block.bg[i])
    }
    .RLpar[["MAX"]] <- prov
    ## OTS data 
    prov <- list()
    nOTS <- 10L
    for (i in 1L:nOTS) {
        prov[[paste("block", i, sep = "")]] <-
            list(col = block.cols[i], pch = block.pch[i], cex = 1.1, lwd = 2,
                 bg = block.bg[i])
    }
    .RLpar[["OTS"]] <- prov
    ## Modify the default values ???
    .RLpar
}

##'  New version
##'  

RLpar <- function(mono = TRUE,
                  trace = 0L, ...) {
    
    ##=======================================================================
    ## Set the graphical parameters for Return Level plots OR CHANGE THEM
    ##========================================================================
    
    mc <- match.call(expand.dots = TRUE)
    
    newPar <- as.list(mc)
    newPar[[1]] <- NULL
    newPar[["mono"]] <- NULL
    newPar[["trace"]] <- NULL
    
    oldPar <- RLpar0(mono = mono)
    if (length(newPar) == 0L) return(oldPar)
    
    dp <- duplicated(names(newPar))
    if (any(dp)) {
        warning("dupplicated par names: ",
                paste("'", names(newPar)[dp], "'", sep = "", collapse = ", "))
    }
    
    uOldPar <- unlist(oldPar)
    
    for (i in 1L:length(newPar)) {
        nm0 <- grep(names(newPar)[i], names(uOldPar))
        if (length(nm0) == 0) {
            warning("no par names matching \"", names(newPar)[i], "\"")
        }
        ## no longer possible
        if (any(is.na(nm0))) {
            warning("'from' contains unused names:\n",
                    sprintf("\"%s\"", names(newPar[[i]])[is.na(nm0)]))
        }
        
        ## Replacement of the values
        Names0 <- names(uOldPar)[nm0]
        ## cat("Names0\n"); print(Names0)
        
        ## This can not be vectorized since nmVec is a vector with elements for
        ## the hierarchical levels 1, 2, ...
        for (Name0 in Names0){
            nmVec <- unlist(strsplit(Name0, split = "\\."))
            if (trace) cat("par.name = ", nmVec, ", value = ", newPar[[i]], "\n")
            
            ## tried 'eval' on 2014-12-01
            oldPar[[nmVec]] <- newPar[[i]]
            ## cat("new value  ",  oldPar[[nmVec]], "\n")
        }
    }
    
    oldPar
}

##*****************************************************************************
`OldRLpar` <- function(mono = TRUE, ...) {
    
    ##=======================================================================
    ## Set the graphical parameters for Return Level plots OR CHANGE THEM
    ##========================================================================
    
    mc <- match.call(expand.dots = TRUE)
    
    newPar <- as.list(mc)
    newPar[[1]] <- NULL
    newPar[["mono"]] <- NULL
    
    oldPar <- RLpar0(mono = mono)
    
    dp <- duplicated(names(newPar))
    
    if (any(dp)) {
        warning("dupplicated par names: ",
                paste("'", names(newPar)[dp], "'", sep = "", collapse = ", "))
    }
    
    uOldPar <- unlist(oldPar)
    nm0 <- match(names(newPar), names(uOldPar))
    if (any(is.na(nm0))) {
        warning("'from' contains unused names:\n",
                sprintf("\"%s\"", names(newPar)[is.na(nm0)]))
    }
    
    ## Replacement of the values
    Names0 <- names(newPar)[!is.na(nm0)]
    
    ## This can not be vectorized since nmVec is a vector with elements for
    ## the hierarchical levels 1, 2, ...
    for (Name0 in Names0){
        nmVec <- unlist(strsplit(Name0, split = "\\."))
        oldPar[[nmVec]] <- newPar[[Name0]]
    }
    
    oldPar
}


##*****************************************************************************
`RLlegend.ini` <- function(x = "topleft", bty = "n", ...) {
    .RLlegend <- list(x = x, bty = bty, ...)
    assign(".RLlegend", .RLlegend, envir = RenextEnvir)
    invisible(.RLlegend)
}

##*****************************************************************************
`RLlegend.show` <- function() {
    .RLlegend <- get(".RLlegend", envir = RenextEnvir, mode = "list")
    do.call("legend", .RLlegend)
}

##*****************************************************************************
## Not exported, but useful for debuging through Renext:::RLlegend.print()
##
`RLlegend.print` <- function() {
    .RLlegend <- get(".RLlegend", envir = RenextEnvir, mode = "list")
    print(.RLlegend)
}

##*****************************************************************************
plot.Renouv <- function(x,
                        pct.conf = x$pct.conf,
                        show = list(OT = TRUE, quant = TRUE, conf = TRUE, MAX = TRUE,
                            OTS = TRUE),
                        mono = TRUE,
                        predict = FALSE, 
                        par = NULL,
                        legend = TRUE,
                        label = NULL,
                        problim = NULL,
                        Tlim = NULL,
                        main = NULL,
                        xlab = "periods",
                        ylab = "level",
                        posOptions = NULL,
                        ## maxBlocks = 10L,
                        byBlockStyle = NULL,
                        ...) {
    
    ##==========================================================================
    ## The plot method
    ##==========================================================================
    
    mc <- match.call(expand.dots = TRUE)
        
    yLim <- rangeLev.Renouv(x, Tlim = Tlim)
    labs <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
    
    if (is.null(main)) main <- ""

    ## compute suitable xlim 
    if (!is.null(problim)) {
      
      if ( !is.numeric(problim) || length(problim) != 2L ||
          any(is.na(problim)) || any(problim <= 0) || any(problim >= 1) ) {
        stop("invalid limits in 'problim'.")
      }
      
      if (!is.null(Tlim)) {
        stop("only one of 'problim' and 'Tlim' arguments can be provided")
      }
      xLim <-  -log(x$estimate["lambda"] * c(1 - problim))
      
    } else {
        if (!is.null(Tlim)) { 
        if ( !is.numeric(Tlim) || length(Tlim) != 2 || any(is.na(Tlim)) ||
            any(Tlim < 0) ) {
          stop("invalid limits in 'Tlim'.")
        }
        if (Tlim[1] < 0.01) Tlim[1] <- 0.01
        xLim <-  log(Tlim)
      } else {
        xMin <- max(c(0,  log(min(x$ret.lev$period))), na.rm = TRUE)
        xLim <- c(xMin, log(max(x$pred$period)))  
      }
    }
    
    ## prepare (empty) plot
    plot(x = xLim,
         y = yLim,
         type = "n",
         main = main,
         xlab = xlab,
         ylab = ylab,
         xaxt = "n",
         ...)
    ## add x-axis and grid
    axis(side = 1, at = log(labs), labels = labs)
    abline(v = log(labs), col = "gray", lty = "dotted")
    abline(h = pretty(par()$usr[3:4]), col = "gray", lty = "dotted")
    ## prepare legend

    if (is.null(label)) label <- deparse(substitute(x))
    if (legend) RLlegend.ini()
  
    ## Call 'lines' method
    lines.Renouv(x = x,
                 pct.conf = pct.conf,
                 show = show,
                 mono = mono,
                 predict = predict,
                 par = par,
                 legend = FALSE,
                 label = label,
                 posOptions = posOptions,
                 ## maxBlocks = maxBlocks,
                 byBlockStyle = byBlockStyle,
                 ...)
    
    ## draw legend
    if (legend) {
        RLlegend.show()
    }
  }

##*****************************************************************************
`lines.Renouv` <- function(x,
                           pct.conf = x$pct.conf,
                           show = NULL,
                           mono = TRUE,
                           predict = FALSE,
                           par = NULL,
                           legend = FALSE,
                           label = NULL,
                           posOptions = NULL,
                           ## maxBlocks = 10L,
                           byBlockStyle = NULL,
                           ...) {
   
    DEBUG <- FALSE
    RLpar0 <- par

    ## After this, 'byBlockStyle' is a list with both elements "MAX" and "OTS"
    ## existing
    if (is.null(byBlockStyle)){
        byBlockStyle <- list()
        if (is.null(x$x.OT)) {
            ## maybe block maxima with OTS data ?
            byBlockStyle$MAX <- FALSE
            byBlockStyle$OTS <- TRUE
        } 
    } else {
        byBlockStyle <- as.list(byBlockStyle)
        if (!(all(names(byBlockStyle) %in% c("MAX", "OTS")))) 
            warning("'byBlockStyle' names differing from \"MAX\" or \"OTS\": ",
                    "ignored")
    }
    ## now fill poissibly missing elements
    for (type in c("MAX", "OTS")) {
        if (is.null(byBlockStyle[[type]])) {
            ht <- x[[paste("history", type, sep = ".")]]
            if (ht$flag) byBlockStyle[[type]] <- (nlevels(ht$block) <= 10L)
            else byBlockStyle[[type]] <- FALSE 
        }
    }
        
    ##*************************************************************************
    ## lines method for objects with class "Renouv" 
    ##*************************************************************************

    nx.OT <- length(x$x.OT)   

    if (is.null(posOptions)) {
        ST <- SandT(x)
    } else {
        parList <- c(list(object = x), posOptions)
        ST <- do.call(SandT, parList)
    }
    
    ##=========================================================================
    ## shown objects
    ##=========================================================================
    .show <- list(quant = TRUE, conf = FALSE, OT = FALSE, MAX = FALSE,
                  OTS = FALSE)
    
    if (!is.null(show)) {
        if (!is.list(show)) stop("'show' must be a list")
        for (type in c("OT", "quant", "conf")) {
            if (!is.null(show[[type]])) .show[[type]] <- show[[type]]
        }
    }
        
    for (type in c("MAX", "OTS")) {
        .show[[type]] <- transShow(show[[type]], x = x, type = type)
        if (byBlockStyle[[type]] && any(.show[[type]]) && !all(.show[[type]]))
            warnings("some but not all blocks of type \"", type, "\" are shown",
                     " with 'byBlockStyle' FALSE, this may be misleading. ",
                     "Consider using 'byBlockStyle' TRUE for this type.")
    }

    ## cat("show after\n"); print(.show)    
    .par <- RLpar(mono = mono)
    if (!is.null(RLpar0)) {
        if (!is.list(RLpar0)) stop("'par' must be a list")
        .par <- rReplace(from = RLpar0, to = .par)
    }
    
    ##=========================================================================
    ## analyse 'label' and prepare legend
    ##=========================================================================
    if (is.null(label)) label <- deparse(substitute(x))
   
    if (is(label, "character")) {

        Label <- list(quant = paste(label, "quant"))
        if ( .show[["conf"]] && length(pct.conf) ) {  ## pct
            Label$conf <- list()
            for (il in 1L:length(pct.conf)) {
                Label$conf[[paste("conf", il, sep = "")]] <-
                    paste(label, " ", pct.conf[il], "%", sep = "")
            }
        }
        
        if (.show[["OT"]] && (nx.OT > 0L) ) Label$OT <- paste(label, "sample")

        ## fill Label[["MAX"]] and Label[["OTS"]] with suitable character vectors
        for (type in c("MAX", "OTS")) {
            ht <- x[[paste("history", type, sep = ".")]]
            if (ht$flag && any(.show[[type]])) {
                nB <- nlevels(ht$block)
                if (nB == 0L) {
                    stop("nlevels(x[[history.", type, "]]$block) is zero")
                }
                if (byBlockStyle[[type]]) {
                    if (all(nchar(ht$blockNames) > 0L)) {
                        Label[[type]] <- as.list(paste(label, ht$blockNames))
                        names(Label[[type]]) <- paste("block", 1L:nB, sep ="")
                    } else {
                        for (ib in 1L:nB) {
                            Label[[type]][[paste("block", ib, sep = "")]] <-
                                paste(label, " ", type, " block", ib,  sep = "")
                        } 
                    } 
                } else {
                    Label[[type]] <- as.list(rep(paste(label, type), nB))
                    names(Label[[type]]) <- paste("block", 1L:nB, sep ="")
                }
            }
        }
        
    } else if (is(label, "list")) Label <- label

    
    ##=========================================================================
    ## predict if necessary
    ##=========================================================================
    ## limits in years with min >= 1 year
    xLim <- par()$usr[c(1L, 2L)]
    if (xLim[1L] < 0)  xLim[1L] <- 0

    .RLlegend <- try(get(".RLlegend", RenextEnvir))
    oc <- class(.RLlegend)
    if ( is.null(.RLlegend) || oc == "try-error" ) {
        warning("'.RLlegend' not found: re-initialising")
        .RLlegend <- RLlegend.ini()
    }
    
    if (.show[["quant"]] || (.show[["conf"]] && length(pct.conf)) ) {
    
        if (predict) {
            x.g <- seq(from = xLim[1L], to = xLim[2L], length.out = 100L)
            if (length(pct.conf) > 0) lev.conf0 <- pct.conf / 100
            else lev.conf0 <- 0.95
            Data <- predict.Renouv(object = x, newdata = exp(x.g),
                                   level = lev.conf0)
        } else {
            Data <- x$ret.lev
        }
        x.g <- log(Data$period) ## x.g is replaced
       
    }
    
    ##=========================================================================
    ## plot the quantile curve
    ##=========================================================================
    if ( .show[["quant"]]) {
        parList <- c(list(x = x.g, y = Data$quant), .par[["quant"]])
        do.call(lines, parList)
        ## XXX
        ## lines(x = x.g, y = Data$quant, col = "green", type = "o", lty = 2)
        ##-----------------------------------------------------------------------
        ## update legend
        ##-----------------------------------------------------------------------
        if (is.null(Label[["quant"]])) {
            stop("'label' must be a character or a list with a 'quant' element")
        }
        .RLlegend$legend <- c(.RLlegend$legend, eval(Label[["quant"]]))
        lty.prov <- eval(.par[["quant"]][["lty"]])
        par(lty = lty.prov)
        lty.prov <- par()$lty
        .RLlegend[["lty"]] <- c(.RLlegend[["lty"]], lty.prov)
        .RLlegend[["pch"]] <- c(.RLlegend[["pch"]], NA)
        .RLlegend[["pt.bg"]] <- c(.RLlegend[["pt.bg"]], NA)
        for (nm in c("col", "lwd")) {
            .RLlegend[[nm]] <- c(.RLlegend[[nm]], eval(.par[["quant"]][[nm]]))
        }
    }
    ##=========================================================================
    ## plot confidence limits
    ##=========================================================================
    if ( .show[["conf"]] && length(pct.conf) ) {
        if (DEBUG) cat("showing 'conf' lines ... \n")
        cnames <- colnames(Data)
        candLnames <- match(paste("L", pct.conf, sep = "."), cnames)
        candUnames <- match(paste("U", pct.conf, sep = "."), cnames)
        ind <- !is.na(candLnames) & !is.na(candUnames)
        for (i in 1L:length(pct.conf))  {
            if (ind[i]) {
                iL <- candLnames[i]
                iU <- candUnames[i]
                parList <- c(list(x = x.g, y = Data[, iL]), .par[["conf"]][[i]])
                do.call(lines, parList)
                parList$y <- Data[ , iU]
                do.call(lines, parList)
                ## update legend
                nm <- paste("conf", i, sep = "")
                if (is.null(Label$conf[[nm]])) {
                    stop("'label' must be a character or a list with a '",
                         nm, "' element")
                }
                .RLlegend$legend <- c(.RLlegend$legend, eval(Label$conf[[nm]]))
                ## translate 'lty' if needed using 'par' function
                lty.prov <- .par[["conf"]][[nm]][["lty"]]
                par(lty = lty.prov)
                lty.prov <- par( )$lty
                .RLlegend[["lty"]] <- c(.RLlegend[["lty"]], eval(lty.prov))
                .RLlegend[["pch"]] <- c(.RLlegend[["pch"]], NA)
                .RLlegend[["pt.bg"]] <- c(.RLlegend[["pt.bg"]], NA)
                for (nm2 in c("col", "lwd")) {
                    .RLlegend[[nm2]] <- c(.RLlegend[[nm2]],
                                          eval(.par[["conf"]][[nm]][[nm2]]))
                }
                if (DEBUG) {
                    cat("YYY\n")
                    print(.RLlegend)
                }
            } else {
                warning("confidence limits for level ",
                        pct.conf[i], "% not found in data. Use 'predict = TRUE'")   
            }
        }
    }
    ##=========================================================================
    ## show sample points if wanted, and if there are some
    ##=========================================================================
    if (.show[["OT"]] && (nx.OT > 0L) ) {
        if (DEBUG) cat("showing 'OT' data ... \n")
        ind <- ST$groupNames[ST$group] == "OT"

        ## if ( (length(selName) > 0L) && !is.null(selName$OT) ) {
        ##     if (is.null(names(ST$x)))
        ##         warning("OT vector with no names. 'selName' ignored")
        ##     else ind <- ind & grepl(selName$OT, names(ST$x))
        ## }
        
        parList <- c(list(x = log(ST$T[ind]), y = ST$x[ind]), .par[["OT"]])
        do.call(points, parList)
        .RLlegend$legend <- c(.RLlegend$legend, eval(Label$OT))
        .RLlegend[["lty"]] <- c(.RLlegend[["lty"]], NA)
        .RLlegend[["col"]] <- c(.RLlegend[["col"]], eval(.par[["OT"]][["col"]]))
        .RLlegend[["lwd"]] <- c(.RLlegend[["lwd"]], NA)
        .RLlegend[["pch"]] <- c(.RLlegend[["pch"]], eval(.par[["OT"]][["pch"]]))
        .RLlegend[["pt.bg"]] <- c(.RLlegend[["pt.bg"]], eval(.par[["OT"]][["bg"]]))
        if (DEBUG) {
            print(.RLlegend)
            cat("... 'OT' data shown\n")
        }
    }
    
    ##=======================================================================
    ## changed on 2014-11-28:  one loop for MAX and OTS data
    ## shorter  and easier to read / maintain
    ##=======================================================================

    typesBlock <- sapply(.show, any)[c("MAX", "OTS")]
    typesBlock <- names(typesBlock[typesBlock])
    
    for (type in typesBlock) {
            
        if (DEBUG) cat("showing '", type, "' data ... \n==================\n")
        maxB <- length(RLpar0()[[type]])
        ht <- x[[paste("history", type, sep = ".")]]
        nB <- nlevels(ht$block)
        ibShown <- 0
        
        for (ib in 1L:nB) {
            
            if (.show[[type]][ib]) {
                
                ibShown <- ibShown + 1L

                if (DEBUG) cat("showing", type, "block", ib, " ... \n")
                nm <- paste(type, ".block", ib, sep = "")
                ind.ib <- ST$groupNames[ST$group] == nm
                
                if (byBlockStyle[[type]]) ibGraph <- 1 + (ib -1L) %% (maxB -1L)
                else ibGraph <- 1L
                nm <- paste("block", ib, sep = "")
                nmGraph <- paste("block", ibGraph, sep = "")
                
                if (ht$r[ib]) {
                    
                    parList <- c(list(x = log(ST$T[ind.ib]), y = ST$x[ind.ib]),
                                 .par[[type]][[ibGraph]])
                    do.call(points, parList)

                    ## update legend
                    if (byBlockStyle[[type]] || (ibShown == 1L)) {

                        if (is.null(Label[[type]][[nm]])) {
                            stop("'Label' must be a character or a list with ",
                                 "a \", type, \" element containing a '", nm,
                                 "' element")
                        }
                        
                        .RLlegend$legend <- c(.RLlegend$legend, Label[[type]][[nm]])
                        
                        ## lines properties.     
                        .RLlegend[["lty"]] <- c(.RLlegend[["lty"]], NA)
                        .RLlegend[["lwd"]] <- c(.RLlegend[["lwd"]], NA)
                
                        ## points properties 'lwd' is for empty symbols pch = 21 to 26 
                        .RLlegend[["col"]] <- c(.RLlegend[["col"]],
                                                eval(.par[[type]][[nmGraph]][["col"]]))
                        .RLlegend[["pch"]] <- c(.RLlegend[["pch"]],
                                                eval(.par[[type]][[nmGraph]][["pch"]]))
                        .RLlegend[["pt.bg"]] <- c(.RLlegend[["pt.bg"]],
                                                  eval(.par[[type]][[nmGraph]][["bg"]])) 
                        .RLlegend[["pt.lwd"]] <-
                            c(.RLlegend[["pt.lwd"]],
                              eval(.par[[type]][[nmGraph]][["lwd"]]))
                    }
                    
                } else {  ## EMPTY block this can happen only for OTS blocks
         
                    ## draw an horizontal segment...
                    segments(x0 = par()$usr[1], x1 = log(ht$effDuration[ib]), 
                             y0 = ht$threshold[ib], y1 = ht$threshold[ib],
                             col = .par[[type]][[ibGraph]]$col,
                             lwd = .par[[type]][[ibGraph]]$lwd)

                    ## and materialise its upper end-point
                    points(x = log(ht$effDuration[ib]), y = ht$threshold[ib],
                           col = .par[[type]][[ibGraph]]$col,
                           pch = .par[[type]][[ibGraph]]$pch,
                           bg = .par[[type]][[ibGraph]]$bg,
                           lwd = .par[[type]][[ibGraph]]$lwd)
                    
                    ## update legend
                    if (byBlockStyle[[type]] || (ibShown == 1L)) {

                        .RLlegend$legend <- c(.RLlegend$legend, eval(Label[[type]][[nm]]))

                        ## lines properties
                        .RLlegend[["lty"]] <- c(.RLlegend[["lty"]], "solid")
                        .RLlegend[["lwd"]] <- c(.RLlegend[["lwd"]], .par[[type]][[ibGraph]]$lwd)
                        
                        ## points properties 'lwd' is for empty symbols pch = 21 to 26 
                        .RLlegend[["col"]] <- c(.RLlegend[["col"]],
                                                eval(.par[[type]][[nmGraph]][["col"]]))
                        .RLlegend[["pch"]] <- c(.RLlegend[["pch"]],
                                                eval(.par[[type]][[nmGraph]][["pch"]]))
                        .RLlegend[["pt.bg"]] <- c(.RLlegend[["pt.bg"]],
                                                  eval(.par[[type]][[nmGraph]][["bg"]])) 
                        .RLlegend[["pt.lwd"]] <- c(.RLlegend[["pt.lwd"]],
                                                   eval(.par[[type]][[nmGraph]][["lwd"]]))
                                                   
                    }
                }
            }
                
        }             
        if (DEBUG) {
            print(.RLlegend)
            cat("... '", type,"'data shown\n")
        }
    }
    
    assign(".RLlegend", .RLlegend, envir = RenextEnvir)
    
    if (DEBUG) {
        print(.RLlegend)
        cat("exiting lines.Renouv\n")
    }
    
    
    ## par(opar)
    
}

