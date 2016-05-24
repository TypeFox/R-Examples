plot.sempls <- function(x, ...){
    col <- list(...)$col
    #if(is.null(col) && require(colorspace)){
    #    for(i in 0:(length(x$model$latent)-1)){
    #        n <- sapply(x$model$blocks, length)
    #        col_tmp <- rainbow_hcl(n[i+1], start = 90+(i*777), end = -30+(i*777))
    #        col <- append(col, col_tmp)
    #    }
    #}
    if(!is.null(col)){
        old_col <- trellis.par.get("superpose.line")$col
        trellis.par.set(superpose.line=list(col=col))
    }
    MVs <- NULL
    wghtev <- x$weights_evolution[!x$weights_evolution$weights %in% c(0, NA),]
    #ymin <- min(x$weights_evolution$iteration)
    #ymax <- max(x$weights_evolution$iteration)
    print(xyplot(weights ~ iteration | LVs,
                   groups=MVs, type="a",
                   data=wghtev,
                   as.table=TRUE,
                   auto.key=list(lines=TRUE,
                                 points=FALSE, # new
                                 space="right",
                                 title="MVs",
                                 ...),
                   main="Evolution of Outer Weights",
                   xlab="Iteration",
                   ylab="Outer Weights",
                   #ylim=ymin:ymax,
                   ...))
    if(!is.null(col)){
        trellis.par.set(superpose.line=list(col=old_col))
    }
    invisible(wghtev)
}

### Alternatives:
#xyplot(weights ~ iteration|LVs, data=tmp2, groups=MVs, type = "a", auto.key =list(space = "right", points = FALSE, lines = TRUE), ylim=range(weights))
#tmp <- tmp[tmp$weights!=0,]
#tmp4$LVs <- factor(tmp4$LVs, levels=ecsi$model$latent)
#xyplot(weights ~ iteration|LVs, data=tmp4, groups=MVs, as.table=TRUE, type="b", auto.key =list(space = "right", points = FALSE, lines = TRUE))
#xyplot(weights ~ iteration, data=tmp4, groups=MVs, as.table=TRUE, type="l", col=1)


### lattice:::densityplot
densityplot.sempls <- function(x, data, use=c("fscores", "prediction", "residuals"), ...){
    use <- match.arg(use)
    if(use=="fscores")         val <- x$factor_scores
    else if(use=="prediction") val <- predict(x)
    else if(use=="residuals")  val <- residuals(x)

    Y <- data.frame(NULL)
    exogenous <- exogenous(x$model)
    for(i in x$model$latent){
        if(i %in% exogenous & use!="fscores") next
        tmp <- data.frame(value=val[,i], name=i)
        Y <- rbind(Y, tmp)
    }
    Y$name <- factor(Y$name, levels=x$model$latent)
    #dots <- list(...)
    if(is.null(list(...)$main)){
        main=paste(deparse(substitute(x)), "\n",
                  ifelse(use=="fscores", "factor scores", use))
    }
    if(is.null(list(...)$sub)){
        sub=paste("Exogenous LVs: ", paste(exogenous, collapse=", "))
    }
    densityplot(~value|name, data=Y, main=main, sub=sub, as.table=TRUE, ...)
    #densityplot(~value|name, data=Y, main=main, sub=sub, as.table=TRUE,
    #            font.sub=1, cex.font=0.5,...)

 }

densityplot.bootsempls <- function(x, data, pattern="beta", subset=NULL, ...){
    ind <- grep(pattern, colnames(x$t))
    ifelse(is.null(subset),
           params <- colnames(x$t)[ind],
           ifelse(is.character(subset), params <- subset, params <- colnames(x$t)[subset])
           )
    Y <- data.frame(NULL)
    for(i in params){
        if(round(var(x$t[,i], na.rm=TRUE), digits=4)==0) next
        tmp <- data.frame(value=x$t[,i], name=i)
        Y <- rbind(Y, tmp)
    }
    densityplot(~value|name, data=Y, as.table=TRUE, ...)
 }

## lattice::parallel (Version < 0.20-5)
## lattice::parallelplot (Version >= 0.20-5)
parallelplot.bootsempls <- function(x, data, pattern="beta", subset=NULL, reflinesAt,
                                col=c("grey", "darkred", "darkred", "black"),
                                lty=c("solid", "solid", "dashed", "dotted"), ...)
{
    ifelse(is.null(subset), ind <- grep(pattern, colnames(x$t)), ind <- subset)
    
    lower <- summary(x, ...)$table$Lower
    upper <- summary(x, ...)$table$Upper
    Y <- rbind(x$t, x$t0, lower, upper, deparse.level=0)
    
    if(!missing(reflinesAt)){
        Y <- rbind(Y, matrix(rep(reflinesAt, each=ncol(x$t)),
                             nrow=length(reflinesAt), byrow=TRUE))
        origin <- c(rep("1resample", x$nboot), "2sample", "3ci", "3ci",
                    rep("4reflines", times=length(reflinesAt)))
        Y <- data.frame(Y, origin)
    }
    else Y <- data.frame(Y, origin=c(rep("1resample", x$nboot),
                              "2sample", "3ci", "3ci"))
    
    parallelplot(~Y[ind], data=Y, groups=origin,
                   common.scale=TRUE, col=col, lty=lty, ...)
}


mvplot <- function(model, ...){
  UseMethod("mvplot", model)
}

mvplot.plsm <- function(model, data, LVs, ask=TRUE, ...){
    try(data <- data[, model$manifest], silent=TRUE)
    if(inherits(data, "try-error")) stop("The 'models' manifest variables must be contained in 'data'!")

    long <- reshape(data, v.names="value",  ids=rownames(data), idvar="ids",
                    times=names(data), timevar="MV", varying=list(names(data) ),
                    direction="long")
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))
    par(ask=ask)
    if(missing(LVs)) indx <- model$latent
    else indx <- LVs
    if(length(indx)==1) par(ask=FALSE)
    charts <- list()
    for(i in indx){
        tab <- as.data.frame(xtabs(~ value + MV ,
                                   data=long[long$MV %in% model$block[[i]],]))
        charts[[i]] <- barchart(Freq ~ value| MV, data=tab, main=i, ...)
        print(charts[[i]])
    }
    invisible(charts)
}

