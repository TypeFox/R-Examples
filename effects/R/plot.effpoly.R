# Plot method for effpoly objects

# modified by Michael Friendly: added ci.style="bands" & alpha.band= arg
# modified by Michael Friendly: added lwd= argument for llines (was lwd=2)
# 2013-11-06: fixed drop dimension when only one focal predictor. John
# 2014-10-10: namespace fixes. John
# 2014-12-05: made key.args more flexible. John
# 2014-03-22: use wide columns by default only when x for legend not set. J. Fox

plot.effpoly <- function(x,
    type=c("probability", "logit"),
    x.var=which.max(levels),
    rug=TRUE,
    xlab,
    ylab=paste(x$response, " (", type, ")", sep=""), 
    main=paste(effect, "effect plot"),
    colors, symbols, lines, cex=1.5, lwd=2,
    factor.names=TRUE, ci.style, band.colors, band.transparency=0.3,
    style=c("lines", "stacked"), 
    confint=(style == "lines" && !is.null(x$confidence.level)), 
    transform.x=NULL, ticks.x=NULL, xlim=NULL,
    ylim, rotx=0, alternating=TRUE, roty=0, grid=FALSE,
    layout, key.args=NULL,
    row=1, col=1, nrow=1, ncol=1, more=FALSE, use.splines=TRUE, ...){     
    ci.style <- if(missing(ci.style)) NULL else 
        match.arg(ci.style, c("bars", "lines", "bands", "none"))
    type <- match.arg(type)
    style <- match.arg(style)
    effect.llines <- llines
    has.se <- !is.null(x$confidence.level) 
    if (confint && !has.se) stop("there are no confidence limits to plot")
    if (style == "stacked"){
        if (type != "probability"){
            type <- "probability"
            warning('type set to "probability" for stacked plot')
        }
        if (confint){
            confint <- FALSE
            warning('confint set to FALSE for stacked plot')
        }
    }
    if (missing(colors)){
        if (style == "stacked"){
            colors <- if (x$model == "multinom") rainbow_hcl(length(x$y.levels))
            else sequential_hcl(length(x$y.levels))
        }
        else colors <- palette()
    }
    if (missing(band.colors)) band.colors <- colors
    if (missing(symbols)) symbols <- 1:length(colors)
    if (missing(lines)) lines <- 1:length(colors)
    .mod <- function(a, b) ifelse( (d <- a %% b) == 0, b, d)
    .modc <- function(a) .mod(a, length(colors))
    .mods <- function(a) .mod(a, length(symbols))
    .modl <- function(a) .mod(a, length(lines))
    effect <- paste(sapply(x$variables, "[[", "name"), collapse="*")
    split <- c(col, row, ncol, nrow)
    n.predictors <- length(names(x$x))
    y.lev <- x$y.lev
    n.y.lev <- length(y.lev)
    ylevel.names <- make.names(paste("prob",y.lev))
    colnames(x$prob) <- colnames(x$logit) <- 
        colnames(x$lower.logit) <- colnames(x$upper.logit) <- 
        colnames(x$lower.prob) <- colnames(x$upper.prob)<- ylevel.names
    x.frame <-as.data.frame(x)
    predictors <- names(x.frame)[1:n.predictors]
    levels <- if (n.predictors==1) length (x.frame[,predictors])
    else sapply(apply(x.frame[, predictors, drop=FALSE], 2, unique), length)
    if (is.character(x.var)) {
        which.x <- which(x.var == predictors)
        if (length(which.x) == 0) stop(paste("x.var = '", x.var, "' is not in the effect.", sep=""))
        x.var <- which.x
    }
    x.vals <- x.frame[, names(x.frame)[x.var]]    
    response <-matrix(0, nrow=nrow(x.frame), ncol=n.y.lev)
    for (i in 1:length(x$y.lev)){
        level <- which(colnames(x$prob)[i] == ylevel.names)
        response[,i] <- rep(x$y.lev[level], length(response[,i]))
    }
    prob <- as.vector(x$prob)
    logit <- as.vector(x$logit)
    response <- as.vector(response)
    if (has.se){
        lower.prob <- as.vector(x$lower.prob)
        upper.prob <- as.vector(x$upper.prob)
        lower.logit <- as.vector(x$lower.logit)
        upper.logit <- as.vector(x$upper.logit)
    }
    response <- factor(response, levels=y.lev)
    data <- data.frame(prob, logit)
    if (has.se) data <- cbind(data, data.frame(lower.prob, upper.prob, lower.logit, upper.logit))
    data[[x$response]] <- response
    for (i in 1:length(predictors)){
        data <-cbind(data, x.frame[predictors[i]])
    }
    levs <- levels(x$data[[predictors[x.var]]])
    n.predictor.cats <- sapply(data[, predictors[-c(x.var)], drop=FALSE], 
        function(x) length(unique(x)))
    if (length(n.predictor.cats) == 0) n.predictor.cats <- 1
    ci.style <- if(is.null(ci.style)) {
        if(is.factor(x$data[[predictors[x.var]]])) "bars" else "bands"} else ci.style
    if( ci.style=="none" ) confint <- FALSE
    ### no confidence intervals if confint == FALSE or ci.style=="none"
    if (!confint){ # plot without confidence bands
        layout <- if (missing(layout)){
            lay <- c(prod(n.predictor.cats[-(n.predictors - 1)]), 
                prod(n.predictor.cats[(n.predictors - 1)]), 1)
            if (lay[1] > 1) lay else lay[c(2, 1, 3)]
        }
        else layout
        if (style == "lines"){ # line plot
            if (n.y.lev > min(c(length(colors), length(lines), length(symbols))))
                warning('Colors, lines and symbols may have been recycled')
            if (is.factor(x$data[[predictors[x.var]]])){ # x-variable a factor
                key <- list(title=x$response, cex.title=1, border=TRUE,
                    text=list(as.character(unique(response))),
                    lines=list(col=colors[.modc(1:n.y.lev)], lty=lines[.modl(1:n.y.lev)], lwd=lwd),
                    points=list(pch=symbols[.mods(1:n.y.lev)], col=colors[.modc(1:n.y.lev)]),
                    columns = if ("x" %in% names(key.args)) 1 else find.legend.columns(n.y.lev))
                for (k in names(key.args)) key[k] <- key.args[k]
                result <- xyplot(eval(if (type=="probability") 
                    parse(text=if (n.predictors==1) 
                        paste("prob ~ as.numeric(", predictors[x.var], ")")
                        else paste("prob ~ as.numeric(", predictors[x.var],") | ", 
                            paste(predictors[-x.var], collapse="*")))
                    else parse(text=if (n.predictors==1) 
                        paste("logit ~ as.numeric(", predictors[x.var], ")")
                        else paste("logit ~ as.numeric(", predictors[x.var],") | ", 
                            paste(predictors[-x.var], collapse="*")))), 
                    strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                    panel=function(x, y, subscripts, rug, z, x.vals, ...){
                        if (grid) panel.grid()
                        for (i in 1:n.y.lev){
                            sub <- z[subscripts] == y.lev[i]
                            good <- !is.na(y[sub])
                            effect.llines(x[sub][good], y[sub][good], lwd=lwd, type="b", col=colors[.modc(i)], lty=lines[.modl(i)],
                                pch=symbols[i], cex=cex, ...)
                        }
                    },
                    ylab=ylab,
                    ylim= if (missing(ylim))
                        if (type == "probability") range(prob) else range(logit)
                    else ylim,
                    xlab=if (missing(xlab)) predictors[x.var] else xlab,
                    x.vals=x$data[[predictors[x.var]]], 
                    rug=rug,
                    z=response,
                    scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx),
                        y=list(rot=roty),  
                        alternating=alternating),
                    main=main,
#                    key=c(key, key.args),
                    key=key,
                    layout=layout,
                    data=data, ...)
                result$split <- split
                result$more <- more
                class(result) <- c("plot.eff", class(result))    		
            }
            else { # x-variable numeric
                if(use.splines) effect.llines <- spline.llines # added 10/17/13
                nm <- predictors[x.var]
                x.vals <- x$data[[nm]]   
                if (nm %in% names(ticks.x)){
                    at <- ticks.x[[nm]]$at
                    n <- ticks.x[[nm]]$n
                }
                else{
                    at <- NULL
                    n <- 5
                }
                xlm <- if (nm %in% names(xlim)){
                    xlim[[nm]]
                }
                else range.adj(data[nm]) # range(x.vals)
                tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
                    trans <- transform.x[[nm]]$trans
                    make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
                }
                else {
                    trans <- I
                    make.ticks(xlm, link=I, inverse=I, at=at, n=n)
                }
                key <- list(title=x$response, cex.title=1, border=TRUE,
                    text=list(as.character(unique(response))), 
                    lines=list(col=colors[.modc(1:n.y.lev)], lty=lines[.modl(1:n.y.lev)], lwd=lwd),
                    columns = if ("x" %in% names(key.args)) 1 else find.legend.columns(n.y.lev))
                for (k in names(key.args)) key[k] <- key.args[k]
                result <- xyplot(eval(if (type=="probability") 
                    parse(text=if (n.predictors==1) paste("prob ~ trans(", predictors[x.var], ")")
                        else paste("prob ~ trans(", predictors[x.var],") |", 
                            paste(predictors[-x.var], collapse="*")))
                    else parse(text=if (n.predictors==1) paste("logit ~ trans(", predictors[x.var], ")")
                        else paste("logit ~ trans(", predictors[x.var],") | ", 
                            paste(predictors[-x.var], collapse="*")))), 
                    strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                    panel=function(x, y, subscripts, rug, z, x.vals, ...){
                        if (grid) panel.grid()
                        if (rug) lrug(trans(x.vals))
                        for (i in 1:n.y.lev){
                            sub <- z[subscripts] == y.lev[i]
                            good <- !is.na(y[sub])
                            effect.llines(x[sub][good], y[sub][good], lwd=lwd, type="l", col=colors[.modc(i)], lty=lines[.modl(i)], ...)
                        }
                    },
                    ylab=ylab,
                    xlim=suppressWarnings(trans(xlm)),
                    ylim= if (missing(ylim))
                        if (type == "probability") range(prob) else range(logit)
                    else ylim,
                    xlab=if (missing(xlab)) predictors[x.var] else xlab,
                    x.vals=x$data[[predictors[x.var]]], 
                    rug=rug,
                    z=response,
                    scales=list(x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx), y=list(rot=roty),
                        alternating=alternating),
                    main=main,
#                    key=c(key, key.args),
                    key=key,
                    layout=layout,
                    data=data, ...)
                result$split <- split
                result$more <- more
                class(result) <- c("plot.eff", class(result))			
            }
        }
        else { # stacked plot
            if (n.y.lev > length(colors))
                stop(paste('Not enough colors to plot', n.y.lev, 'regions'))
            key <- list(text=list(lab=rev(y.lev)), rectangle=list(col=rev(colors[1:n.y.lev])))
            for (k in names(key.args)) key[k] <- key.args[k]
            if (is.factor(x$data[[predictors[x.var]]])){ # x-variable a factor
                result <- barchart(eval(parse(text=if (n.predictors == 1) 
                    paste("prob ~ ", predictors[x.var], sep="")
                    else paste("prob ~ ", predictors[x.var]," | ", 
                        paste(predictors[-x.var], collapse="*")))), 
                    strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                    groups = response,
                    col=colors,
                    horizontal=FALSE, 
                    stack=TRUE, 
                    data=data, 
                    ylim=if (missing(ylim)) 0:1 else ylim,
                    ylab=ylab, 
                    xlab=if (missing(xlab)) predictors[x.var] else xlab,
                    scales=list(x=list(rot=rotx), y=list(rot=roty), 
                        alternating=alternating),
                    main=main,
#                    key=c(key, key.args),
                    key=key,
                    layout=layout)
                result$split <- split
                result$more <- more
                class(result) <- c("plot.eff", class(result))			
            }
            else { # x-variable numeric
                if(use.splines) effect.llines <- spline.llines # added 10/17/13
                nm <- predictors[x.var]
                x.vals <- x$data[[nm]]   
                if (nm %in% names(ticks.x)){
                    at <- ticks.x[[nm]]$at
                    n <- ticks.x[[nm]]$n
                }
                else{
                    at <- NULL
                    n <- 5
                }
                xlm <- if (nm %in% names(xlim)){
                    xlim[[nm]]
                }
                else range.adj(data[nm]) # range(x.vals)
                tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
                    trans <- transform.x[[nm]]$trans
                    make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
                }
                else {
                    trans <- I
                    make.ticks(xlm, link=I, inverse=I, at=at, n=n)
                }
                result <- densityplot(eval(parse(text=if (n.predictors == 1)
                    paste("~ trans(", predictors[x.var], ")", sep="")
                    else paste("~ trans(", predictors[x.var], ") | ",
                        paste(predictors[-x.var], collapse="*")))),
                    probs=x$prob,
                    strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                    panel =  function(x, subscripts, rug, x.vals, probs=probs, col=colors, ...){
                        fill <- function(x, y1, y2, col){
                            if (length(y2) == 1) y2 <- rep(y2, length(y1))
                            if (length(y1) == 1) y1 <- rep(y1, length(y2))
                            panel.polygon(c(x, rev(x)), c(y1, rev(y2)), col=col)
                        }
                        n <- ncol(probs)
                        Y <- t(apply(probs[subscripts,], 1, cumsum))
                        fill(x, 0, Y[,1], col=col[1])
                        for (i in 2:n){
                            fill(x, Y[,i-1], Y[,i], col=col[i])
                        }
                        if (rug) lrug(trans(x.vals))
                    },
                    rug=rug,
                    x.vals=x$data[[predictors[x.var]]],
                    data=x$x,
                    xlim=suppressWarnings(trans(xlm)),
                    ylim=if (missing(ylim)) 0:1 else ylim,
                    ylab=ylab,
                    xlab=if (missing(xlab)) predictors[x.var] else xlab,
                    scales=list(x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx), y=list(rot=roty),
                        alternating=alternating),
                    main=main,
#                    key=c(key, key.args),
                    key=key,
                    layout=layout, ...)
                result$split <- split
                result$more <- more
                class(result) <- c("plot.eff", class(result))
            }
        }
    }
    ### with confidence bands
    else{ # plot with confidence bands
        layout <- if(missing(layout)) c(prod(n.predictor.cats), length(levels(response)), 1) 
        else layout
        if (type == "probability"){
            lower <- lower.prob
            upper <- upper.prob
        }
        else {
            lower <- lower.logit
            upper <- upper.logit
        }
        ### factor
        if (is.factor(x$data[[predictors[x.var]]])){ # x-variable a factor
            levs <- levels(x$data[[predictors[x.var]]])
            result <- xyplot(eval(if (type=="probability") 
                parse(text=if (n.predictors==1) 
                    paste("prob ~ as.numeric(", predictors[x.var],") |", x$response)
                    else paste("prob ~ as.numeric(", predictors[x.var],") |", 
                        paste(predictors[-x.var], collapse="*"), 
                        paste("*", x$response)))
                else parse(text=if (n.predictors==1) 
                    paste("logit ~ as.numeric(", predictors[x.var],") |", x$response)
                    else paste("logit ~ as.numeric(", predictors[x.var],")|", 
                        paste(predictors[-x.var], collapse="*"), 
                        paste("*", x$response)))),
                par.strip.text=list(cex=0.8),							
                strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                panel=function(x, y, subscripts, x.vals, rug, lower, upper, ... ){
                    if (grid) panel.grid()
                    good <- !is.na(y)
                    effect.llines(x[good], y[good], lwd=lwd, type="b", pch=19, col=colors[1], cex=cex, ...)
                    subs <- subscripts+as.numeric(rownames(data)[1])-1		
                    if (ci.style == "bars"){
                        larrows(x0=x[good], y0=lower[subs][good], 
                            x1=x[good], y1=upper[subs][good], 
                            angle=90, code=3, col=colors[.modc(2)], length=0.125*cex/1.5)
                    }
                    else  if(ci.style == "lines"){
                        effect.llines(x[good], lower[subs][good], lty=2, col=colors[.modc(2)])
                        effect.llines(x[good], upper[subs][good], lty=2, col=colors[.modc(2)])
                    }
                    else { if(ci.style == "bands") {		
                        panel.bands(x[good], y[good],
                            lower[subs][good], upper[subs][good],
                            fill=band.colors[1], alpha=band.transparency)
                    }}
                },
                
                
                ylab=ylab,
                ylim= if (missing(ylim)) c(min(lower), max(upper)) else ylim,
                xlab=if (missing(xlab)) predictors[x.var] else xlab,
                main=main,
                x.vals=x$data[[predictors[x.var]]],
                rug=rug,
                lower=lower,
                upper=upper, 
                scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx), 
                    y=list(rot=roty), alternating=alternating),
                layout=layout,
                data=data, ...)
            result$split <- split
            result$more <- more
            class(result) <- c("plot.eff", class(result))
        }
        else { # x-variable numeric
            if(use.splines) effect.llines <- spline.llines # added 10/17/13
            nm <- predictors[x.var]
            x.vals <- x$data[[nm]]   
            if (nm %in% names(ticks.x)){
                at <- ticks.x[[nm]]$at
                n <- ticks.x[[nm]]$n
            }
            else{
                at <- NULL
                n <- 5
            }
            xlm <- if (nm %in% names(xlim)){
                xlim[[nm]]
            }
            else range.adj(data[nm]) # range(x.vals)
            tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
                trans <- transform.x[[nm]]$trans
                make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
            }
            else {
                trans <- I
                make.ticks(xlm, link=I, inverse=I, at=at, n=n)
            }
            result <- xyplot(eval(if (type=="probability") 
                parse(text=if (n.predictors==1) 
                    paste("prob ~ trans(", predictors[x.var],") |", x$response)
                    else paste("prob ~ trans(", predictors[x.var],") |", 
                        paste(predictors[-x.var], collapse="*"), 
                        paste("*", x$response)))
                else parse(text=if (n.predictors==1) 
                    paste("logit ~ trans(", predictors[x.var],") |", x$response)
                    else paste("logit ~ trans(", predictors[x.var],") |", 
                        paste(predictors[-x.var], collapse="*"), 
                        paste("*", x$response)))
            ),
                par.strip.text=list(cex=0.8),							
                strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                panel=function(x, y, subscripts, x.vals, rug, lower, upper, ... ){
                    if (grid) panel.grid()
                    if (rug) lrug(trans(x.vals))
                    good <- !is.na(y)
                    effect.llines(x[good], y[good], lwd=lwd, col=colors[1], ...)
                    subs <- subscripts+as.numeric(rownames(data)[1])-1	
                    if (ci.style == "bars"){
                        larrows(x0=x[good], y0=lower[subs][good], 
                            x1=x[good], y1=upper[subs][good], 
                            angle=90, code=3, col=colors[.modc(2)], length=0.125*cex/1.5)
                    }
                    else  if(ci.style == "lines"){
                        effect.llines(x[good], lower[subs][good], lty=2, col=colors[.modc(2)])
                        effect.llines(x[good], upper[subs][good], lty=2, col=colors[.modc(2)])
                    } 
                    else { if(ci.style == "bands") {	
                        panel.bands(x[good], y[good],
                            lower[subs][good], upper[subs][good],
                            fill=band.colors[1], alpha=band.transparency)
                    }}
                },
                ylab=ylab,
                xlim=suppressWarnings(trans(xlm)),
                ylim= if (missing(ylim)) c(min(lower), max(upper)) else ylim,
                xlab=if (missing(xlab)) predictors[x.var] else xlab,
                main=main,
                x.vals=x$data[[predictors[x.var]]],
                rug=rug,
                lower=lower,
                upper=upper, 
                scales=list(y=list(rot=roty), x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx),
                    alternating=alternating),
                layout=layout,
                data=data, ...)
            result$split <- split
            result$more <- more
            class(result) <- c("plot.eff", class(result))
        }
    }
    result
}

