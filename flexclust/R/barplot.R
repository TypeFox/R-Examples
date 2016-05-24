#
#  Copyright (C) 2006-2013 Friedrich Leisch
#  $Id: barplot.R 14 2013-07-02 09:56:24Z leisch $
#

setMethod("barplot", "kccasimple",
function (height, bycluster = TRUE, oneplot = TRUE,
          data = NULL, FUN=colMeans, 
          main = deparse(substitute(height)), 
          which = 1:height@k,
          names.arg = NULL, oma=par("oma"),
          col=NULL, mcol="darkred", srt=45, horiz=FALSE, ...)
{
    object <- height
    opar <- par(c("mfrow", "oma", "mgp", "xpd", "oma"))
    on.exit(par(opar))
    n <- length(which)

    if(is.null(col))
        col <- LightColors
    
    col <- rep(col, length=object@k)

    if(is.null(data)){
        cluster <- object@cluster
        centers <- object@centers
        size <- info(object, "size")
        if(is(object, "kcca"))
            datacent <- object@xcent
        else
            datacent <- NULL
    }
    else{
        cluster <- predict(object, data)
        centers <- matrix(NA, nrow=object@k, ncol=ncol(data))
        colnames(centers) <- colnames(data)
        for(k in 1:object@k){
            ok <- cluster==k
            if(any(ok))
                centers[k,] <- FUN(data[ok,])
        }
        size <- tabulate(cluster, nbins=object@k)
        datacent <- FUN(data)
    }
    
    ylim <- range(centers)
    
    if (is.null(names.arg))
        names.arg <- colnames(centers)
    if (bycluster) {
        par(xpd=NA)
        if (oneplot) {
            if (n <= 3) {
                par(mfrow = c(n, 1), oma=oma)
            }
            else {
                par(mfrow = c(ceiling(n/2), 2), oma=oma)
            }
        }
        for (k in which) {
            mid <- barplot(centers[k, ], col = col[k], 
                           names.arg = "", ylim = ylim, ...)
            if (!is.null(names.arg)){
                text(mid + 0.005 * min(srt, 90-srt) * (mid[2] - mid[1]),
                     ylim[1] - par("cxy")[2], adj = ifelse(srt==0, 0.5, 1),
                     srt = srt, paste(names.arg, "  "))
            }
            if (!is.null(mcol) && !is.null(datacent)) {
                points(mid, datacent, pch = 16, col = mcol)
                points(mid, datacent, type = "h", col = mcol)
            }
            title(main = paste("Cluster ", k, ": ", size[k], 
                " points (", round(100 * size[k]/sum(size)), "%)", sep = ""))
        }
    }
    else {
        a <- ceiling(sqrt(ncol(centers)))
        if (oneplot) {
            par(mfrow = c(a, ceiling(ncol(centers)/a)))
        }
        for (k in 1:ncol(centers)) {
            barplot(centers[which, k], col = col[which], 
                    ylim = ylim, xlab="", ...)
            title(main = names.arg[k])
            if (!is.null(mcol) && !is.null(datacent)) {
                abline(h = datacent[k], col = mcol)
            }
        }
    }
})

###**********************************************************

setMethod("barchart", "kccasimple",
function(x, data,  xlab="", strip.labels=NULL, strip.prefix="Cluster ",
         col=NULL, mcol="darkred", mlcol=mcol, which=NULL, legend=FALSE,
         shade=FALSE, diff=NULL, ...)
{
    if(is.null(strip.labels)){
        SIZE <- info(x, "size")
        strip.labels <-
            paste(strip.prefix, 1:x@k, ": ", SIZE, " (",
                  round(100 * SIZE/sum(SIZE)), "%)", sep="")
    }

    if(is.null(mcol)) mcol <- NA
    if(is.null(mlcol)) mlcol <- NA

    b <- Barchart(x=x@centers, m=x@xcent, strip.labels=strip.labels,
                  xlab=xlab, col=col, mcol=mcol, mlcol=mlcol, which=which,
                  shade=shade, diff=diff, ...)

    if(legend){
        plot(b, more=TRUE)

        pushViewport(viewport(x=0.2,y=0, width=0.8,height=0.05,
                              just=c("left","bottom")))

        grid.text("Population center:", 0.1, 0.5, just=1)
        grid.segments(x0=0.12, y0=0.5, x1=0.2, y1=0.5,
                      gp=gpar(col=mlcol))
        grid.points(0.2, 0.5, pch=16,
                    size=unit(0.5, "char"), gp=gpar(col=mcol))

        grid.text("Cluster centers :", 0.48, 0.5, just=1)

        if(shade){        
            grid.rect(0.5,0.75,width=0.1,height=0.25, just=c(0,0.5),
                      gp=gpar(fill=flxColors(color="light")[3]))
            grid.rect(0.5,0.25,width=0.1,height=0.25, just=c(0,0.5),
                      gp=gpar(col=flxColors(color="dark", grey=TRUE)))

            grid.text("relevant difference", 0.62, 0.75, just=0)
            grid.text("irrelevant difference", 0.62, 0.25, just=0)
        }
        else{
            grid.rect(0.5, 0.5, width=0.1, height=0.25, just=c(0,0.5),
                      gp=gpar(fill=col[1]))
        }
        popViewport(1)
        lattice:::lattice.setStatus(print.more = FALSE)
        return(NULL)
    }
    else
        return(b)
})


### Currently not exported, hence document arguments here:
### x: matrix of cluster centers
### m: vector with location of total population
### labels: text for panel header strips (default is rownames(x))
### REST: see barchart method for kccasimple objects
Barchart <- function(x, m, which=NULL, col=NULL, mcol="darkred",
                     mlcol=mcol, strip.labels=NULL, xlab="",
                     shade=FALSE, diff=NULL, ...)
{
    x <- as.matrix(x)
    m <- as.vector(m)

    if(!is.null(strip.labels))
        rownames(x) <- rep(strip.labels, length=nrow(x))
    
    if(is.null(col))
        col <- flxColors(color="light")

    col <- rep(col, length=nrow(x))

    if(is.null(which))
        which <- seq(1, ncol(x))

    x <- x[,which]
    ## sonst musz man die barplots von unten nach oben lesen
    x <- x[,ncol(x):1]
    m <- rev(m[which])

    x <- as.data.frame(as.table(x))

    panel <- createBarchartPanel(m=m, col=col, mcol=mcol, mlcol=mlcol,
                                 shade=shade, diff=diff)

    barchart(Var2~Freq|Var1, data=x, panel=panel, as.table=TRUE,
             xlab=xlab, ...)
}


createBarchartPanel <- function(m, col, mcol, mlcol, shade, diff)
{
    KKK <- 1
    KKKplus <- function() KKK <<- KKK+1

    if(is.null(diff))
        diff <- c(max(m)/4, 0.5)
    else
        diff <- rep(diff, length=2)
            
    grey <- flxColors(color="dark", grey=TRUE)


    mypanel <- function(x, y, ...)
    {
        COL <- rep("white", length(x))
        MCOL <- rep(grey, length=length(x))
        MLCOL <- rep(grey, length=length(x))
        BCOL <- rep(grey, length=length(x))
            
        if(length(shade)==1){
            if(shade){
                d1 <- abs(x-m) >= diff[1]
                d2 <- abs((x-m)/m) >= diff[2]
                shade <- d1|d2
            }
            else{
                shade <- rep(TRUE, length(x))
            }
        }
        else{
            if(is.matrix(shade)) shade <- shade[KKK,]
            ### reverse to match reversing in Barchart() above
            shade <- rev(rep(as.logical(shade), length=length(x)))
        }
        
        COL[shade] <- col[KKK]
        MCOL[shade] <- mcol
        MLCOL[shade] <- mlcol
        BCOL[shade] <- "black"
        
        MCOL[is.na(x)] <- NA
        MLCOL[is.na(x)] <- NA
        
        if(!all(is.na(MLCOL)))
            grid.segments(x0=0, y0=1:length(x), x1=m, y1=1:length(x),
                          gp=gpar(col=MLCOL),
                          default.units="native")

        panel.barchart(x, y, col=COL, border=BCOL, ...)

        if(!all(is.na(MCOL)))
            grid.points(m, 1:length(x), pch=16,
                        size=unit(0.5, "char"), gp=gpar(col=MCOL))
        
        grid.segments(1, 1, 4, 4)
        KKKplus()
    }
    return(mypanel)
}

###**********************************************************

propBarchart <- function(x, g, alpha=0.05, correct="holm",
                         test="prop.test", sort=FALSE,
                         strip.prefix="", strip.labels=NULL,
                         which=NULL, ...)
{
    call <- match.call()
    x <- as(x, "matrix")
    if(sort && is.null(which))
        which <- rev(order(colMeans(x)))

    if(!is.null(which))
        x <- x[,which,drop=FALSE]
    storage.mode(x) <- "integer"
    if(!all.equal(sort(unique(as.vector(x))), 0:1))
        stop("x must be a binary matrix")
    
    g <- as.factor(g)
    b <- 100 * as.matrix(aggregate(x, list(g), mean, na.rm=TRUE)[,-1])
    rownames(b) <- levels(g)

    ltab <- table(g)
    if(is.null(strip.labels))
        strip.labels <- paste(strip.prefix, names(ltab), ": ", ltab, sep="")
    else
        if(length(unique(strip.labels))!=nrow(b))
            stop("need as many unique strip.labels as non-empty groups in g")

    if(is.character(test)) test <- get(test, mode="function")
    
    p <- pa <- apply(x, 2, function(z) test(table(g, z))$p.value)
    if(!is.null(correct))
        pa <- p.adjust(p, method=correct)

    m <- 100*colMeans(x, na.rm=TRUE)

    cpval <- format.pval(pa, digits=3)
    cpval[pa>alpha] <- "."
    TAB <- cbind(t(round(b)), all=round(m), p.value=cpval)

    new("propBarchart", chart=Barchart(b, m, strip.labels=strip.labels,
                                       shade=pa<=alpha, col="grey", ...),
        gprop = b, tprop = m, p.value=pa, table=TAB)
}

setMethod("show", "propBarchart", function(object) plot(object@chart))

setMethod("summary", "propBarchart",
function(object, ...)
{
    print(object@table, quote=FALSE, ...)
})
