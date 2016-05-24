#' Function to visualise a data frame using advanced boxplot
#'
#' \code{visBoxplotAdv} is supposed to visualise a data frame using advanced boxplot. In addition to boxplot, a scatter plot is also drawn with various methods to avoid co-incident points so that each point is visible (with fine-controling the color and plotting character). Also, these points can be pies or thermometers, which allows an additional proportation data to be visualised as well. 
#'
#' @param formula a formula, such as 'y ~ grp', where 'y' is a numeric vector of data values to be split into groups according to the grouping variable 'grp' (usually a factor)
#' @param data a data.frame (or list) from which the variables in 'formula' should be taken.
#' @param orientation the orientation. It can be one of "vertical" for the vertical orientation, "horizontal" for the horizontal orientation
#' @param method the method for arranging the points. It can be one of "swarm" for arranging points in increasing order (if a point would overlap an existing point, it is shifted sideways (along the group axis) by a minimal amount sufficient to avoid overlap), "center" for first discretizing the values along the data axis (in order to create more efficient packing) and then using a square grid to produce a symmetric swarm, "hex" for first discretization and then arranging points in a hexagonal grid, and "square" for first discretization and then arranging points in a square grid
#' @param corral the method to adjust points that would be placed outside their own group region. It can be one of "none" for not adjusting runaway points, "gutter" for collecting runaway points along the boundary between groups, "wrap" for wrapping runaway points to produce periodic boundaries, "random" for placing runaway points randomly in the region, and "omit" for omitting runaway points
#' @param corralWidth the width of the "corral" in user coordinates
#' @param cex size of points relative to the default. This must be a single value
#' @param breaks breakpoints (optional).  If NULL, breakpoints are chosen automatically
#' @param spacing relative spacing between points
#' @param labels labels for each group. Recycled if necessary. By default, these are inferred from the data
#' @param at numeric vector giving the locations where the swarms should be drawn; defaults to '1:n' where n is the number of groups
#' @param add whether to add to an existing plot
#' @param log whether to use a logarithmic scale on the data axis
#' @param xlim limits for x-axis
#' @param ylim limits for y-axis
#' @param xlab labels for x-aixs
#' @param ylab labels for y-aixs
#' @param pch plotting characters, specified by group and recycled if necessary. In additon to the convertional pch values, it can also be "circles", "thermometers", or "pies". For "pies" (or "thermometers"), users can also specify the proportional values (see below "pwpie") to visualise another information in the pie (or themometer) chart
#' @param col plotting colors, specified by group and recycled if necessary
#' @param bg plotting background, specified by group and recycled if necessary
#' @param pwpch point-wise version of pch
#' @param pwcol point-wise version of col
#' @param pwbg point-wise version of bg
#' @param pwpie point-wise proportion used when drawing pies or themometers
#' @param do.plot whether to draw main plot
#' @param do.boxplot whether to draw boxplot. It only works when the main plot is drawn
#' @param boxplot.notch whether to draw a notch in the boxplot. If the notches of two plots do not overlap this is 'strong evidence' that the two medians differ
#' @param boxplot.border the color for the outlines of the boxplots
#' @param boxplot.col the color for the bodies of the boxplots
#' @param ... additional graphic parameters for the plot
#' @return
#' A data frame with plotting information. It has the same row names as the input data
#' @note none
#' @export
#' @importFrom Biobase pData
#' @seealso \code{\link{visBoxplotAdv}}
#' @include visBoxplotAdv.r
#' @examples
#' data(TCGA_mutations)
#' pd <- Biobase::pData(TCGA_mutations)
#' # only tumor types "LAML" or "BLCA"
#' data <- pd[pd$TCGA_tumor_type=="LAML" | pd$TCGA_tumor_type=="BLCA",]
#' labels <- levels(as.factor(data$TCGA_tumor_type))
#' # colors for gender
#' pwcol <- as.numeric((data$Gender))
#' # pie for relative age
#' pwpie <- data$Age/(max(data$Age))
#' out <- visBoxplotAdv(formula=time ~ TCGA_tumor_type, data=data, pch="pies", pwcol=pwcol, pwpie=pwpie)
#' legend("topright", legend=levels(data$Gender), box.col="transparent", pch=19, col=unique(pwcol))

visBoxplotAdv <- function(formula, data, orientation=c("vertical","horizontal"), method=c("center","hex","square","swarm"), corral=c("none","gutter","wrap","random","omit"), corralWidth, cex=1, spacing=1, breaks=NULL, labels, at=NULL, add=FALSE, log=FALSE, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, pch=c("circles","thermometers","pies")[1], col=graphics::par("col"), bg=NA, pwpch=NULL, pwcol=NULL, pwbg=NULL, pwpie=NULL, do.plot=TRUE, do.boxplot=TRUE, boxplot.notch=FALSE, boxplot.border="#888888C0", boxplot.col="transparent", ...)
{

    if(missing(formula) || (length(formula)!=3)){
        stop("The input formula is missing or incorrect!")
    }
    
    if(is.null(pwpie) || length(pwpie) != nrow(data)) pwpie <- rep(1, nrow(data))
    
    method <- match.arg(method)
    corral <- match.arg(corral)
    orientation <- match.arg(orientation)
    
    m <- match.call(expand.dots=FALSE)
    if (is.matrix(eval(m$data, parent.frame()))){
        m$data <- as.data.frame(data)
    }
    
    if(is.null(rownames(data))){
        rnames <- 1:nrow(data)
    }else{
        rnames <- rownames(data)
    }
    rownames(data) <- rnames
    
    m[[1]] <- as.name("model.frame")
    flag <- names(m) %in% c("pwpch","pwcol","pwbg","pwpie")
    m <- m[union(1:3,which(flag))]
    
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    f <- mf[-response]
    
    f <- f[names(f) %in% attr(attr(mf, "terms"), "term.labels")]
    if(!is.null(mf$'(pwpch)')){
        pwpch <- split(mf$'(pwpch)', f)
    }
    if(!is.null(mf$'(pwcol)')){
        pwcol <- split(mf$'(pwcol)', f)
    }
    if(!is.null(mf$'(pwbg)')){
        pwbg <- split(mf$'(pwbg)',f)
    }
    
    ##################################
    if(!is.null(mf$'(pwpie)')){
        pwpie <- split(mf$'(pwpie)',f)
    }
    
    pwrnames <- split(rnames,f)
    ##################################
    
    dlab <- as.character(formula)[2]
    glab <- as.character(formula)[3]
    
    ## divides the data in the vector 'x' into the groups defined by 'f'
    ## a list of vectors containing the values for the groups
    ## components of the list are named by the levels of 'f'
    x <- split(x=mf[[response]], f=f)
    n.groups <- length(x)
    
    if(length(cex) > 1) {
        stop('the parameter "cex" must have length 1')
    }
    
    # group labels
    if(missing(labels) || is.null(labels)) {
        if(is.null(names(x))) {
            if(n.groups == 1) {
                labels <- NA
            }else{
                labels <- 1:n.groups
            }
        }else{
            labels <- names(x)
        }
    }else{
        labels <- rep(labels, length.out=n.groups)
    }
    
    # at-labels
    if (is.null(at)){
        at <- 1:n.groups
    }else if (length(at) != n.groups){
        stop(gettextf("'at' must have length equal to %d, the number of groups", n.groups), domain = NA)
    }
    
    ## this function returns a "group" vector, to complement "unlist"
    unlistGroup <- function(x, nms=names(x)){
        ## If times consists of a single integer, the result consists of the whole input repeated this many times.
        ## If times is a vector of the same length as x (after replication by each), the result consists of x[1] repeated times[1] times, x[2] repeated times[2] times and so on
        rep(nms, sapply(x, length))
    }

    x.val <- unlist(x)
    x.gp <- unlistGroup(x, nms=labels)
    if((range(x.val, finite=T)[1] <= 0) && log) warning('values <= 0 omitted from logarithmic plot')

    n.obs <- length(x.val)
    n.obs.per.group <- sapply(x, length)

    ## Resolve xlim, ylim, dlim, xlab, ylab
    if(log) {
        dlim <- 10 ^ (grDevices::extendrange(log10(x.val[x.val > 0])))
    }else{
        dlim <- grDevices::extendrange(x.val, f=0.01)
    }
    
    glim <- c(min(at)-0.5, max(at)+0.5)
    
    if(orientation=="horizontal") { 
        if(is.null(ylim)){
            ylim <- glim
        }
        if(is.null(xlim)){
            xlim <- dlim
        }
        if(is.null(xlab)){
            xlab <- dlab
        }
        if(is.null(ylab)){
            ylab <- glab
        }
    }else if(orientation=="vertical"){
        if(is.null(xlim)){
            xlim <- glim
        }
        if(is.null(ylim)){
            ylim <- dlim
        }
        if(is.null(ylab)){
            ylab <- dlab
        }
        if(is.null(xlab)){
            xlab <- glab
        }
    }

    #### Resolve plotting characters and colors
    if(is.null(pwpch)) {
        pch.out <- unlistGroup(x, nms=rep(pch, length.out=n.groups))
    }else{
        if(is.list(pwpch)) {
            names(pwpch) <- names(x)
            stopifnot(all(sapply(pwpch, length)==n.obs.per.group))
            pch.out <- unlist(pwpch)
        }else{
            pch.out <- pwpch
        }
    }
    stopifnot(length(pch.out) == n.obs)

    if(is.null(pwcol)) {
        col.out <- unlistGroup(x, nms=rep(col, length.out=n.groups))
    } else {
        if(is.list(pwcol)) {
            names(pwcol) <- names(x)
            stopifnot(all(sapply(pwcol, length) == n.obs.per.group))
            col.out <- unlist(pwcol)
        } else {
            col.out <- pwcol
        }
    }
    stopifnot(length(col.out) == n.obs)

    if(is.null(pwbg)) {
        bg.out <- unlistGroup(x, nms=rep(bg, length.out=n.groups))
    } else {
        if(is.list(pwbg)) {
            names(pwbg) <- names(x)
            stopifnot(all(sapply(pwbg, length) == n.obs.per.group))
            bg.out <- unlist(pwbg)
        } else {
            bg.out <- pwbg
        }
    }
    stopifnot(length(bg.out) == n.obs)
    
    #######################################
    if(is.null(pwpie)) {
        pie.out <- unlistGroup(x, nms=rep(1, length.out=n.groups))
    } else {
        if(is.list(pwpie)) {
            names(pwpie) <- names(x)
            stopifnot(all(sapply(pwpie, length) == n.obs.per.group))
            pie.out <- unlist(pwpie)
        } else {
            pie.out <- pwpie
        }
    }
    stopifnot(length(pie.out) == n.obs)
    #######################################
    names(pwrnames) <- names(x)    
    rnames.out <- unlist(pwrnames)


    #### Set up the plot
    if(do.plot & !add) {
        plot(xlim, ylim, type='n', axes=F, log=ifelse(log, ifelse(orientation=="horizontal", 'x', 'y'), ''), xlab=xlab, ylab=ylab)
    }

    #### Calculate the size of a plotting character along group- or data-axis
    sizeMultiplier <- graphics::par('cex') * cex * spacing
    if(orientation=="horizontal") {
        size.g <- graphics::yinch(0.08, warn.log = FALSE) * sizeMultiplier
        size.d <- graphics::xinch(0.08, warn.log = FALSE) * sizeMultiplier
    } else {# vertical
        size.g <- graphics::xinch(0.08, warn.log = FALSE) * sizeMultiplier
        size.d <- graphics::yinch(0.08, warn.log = FALSE) * sizeMultiplier
    }


    ############################################
    ## A function called 'calculateSwarm'
    calculateSwarm <- function(x, dsize, gsize) {
        if(length(x) == 0) return(numeric(0))
        out <- data.frame(x = x / dsize, y = 0, i = seq(along = x))
        out <- out[order(out$x), ]
        if(nrow(out) > 1) {
            for (i in 2:nrow(out)) {
                xi <- out$x[i]
                yi <- out$y[i]
                pre <- out[1:(i - 1), ] # previous points
                wh <- xi - pre$x < 1  # which ones are potentially overlapping
                wh[is.na(wh)] <- FALSE  # missing values are not potentially overlapping
                if(any(wh)) {
                    pre <- pre[wh, ]
                    pre <- pre[order(abs(pre$y)), ]
                    poty.off <- sqrt(1 - ((xi - pre$x) ^ 2)) # potential y offset
                    poty <- c(0, pre$y + poty.off, pre$y - poty.off) # potential y values
                    poty.bad <- sapply(poty, function(y) { # check for overlaps
                        any(((xi - pre$x) ^ 2 + (y - pre$y) ^ 2) < 0.999)
                    })
                    poty[poty.bad] <- Inf
                    out$y[i] <- poty[which.min(abs(poty))]
                } else {
                    out$y[i] <- 0
                }
            }
        }
        out <- out[order(out$i), ]
        out[is.na(out$x), 'y'] <- NA  # missing x values should have missing y values
        out$y * gsize
    }


    # jitter points horizontally
    swarmX <- function(x, y, xsize = graphics::xinch(0.08, warn.log = FALSE), ysize = graphics::yinch(0.08, warn.log = FALSE), log = NULL, cex = graphics::par("cex")){ 
        if(is.null(log)) {
            log <- paste(ifelse(graphics::par('xlog'), 'x', ''), ifelse(graphics::par('ylog'), 'y', ''), sep = '')
        }
        xlog <- 'x' %in% strsplit(log, NULL)[[1L]]
        ylog <- 'y' %in% strsplit(log, NULL)[[1L]]
        xy <- grDevices::xy.coords(x = x, y = y, recycle = TRUE, log = log)
        stopifnot((length(unique(xy$x)) <= 1))
        if(xlog) xy$x <- log10(xy$x)
        if(ylog) xy$y <- log10(xy$y)
        x.new <- xy$x + calculateSwarm(xy$y, dsize = ysize * cex, gsize = xsize * cex)
        out <- data.frame(x = x.new, y = y)
        if(xlog) out$x <- 10 ^ out$x
        out
    }

    # jitter points vertically
    swarmY <- function(x, y, xsize = graphics::xinch(0.08, warn.log = FALSE), ysize = graphics::yinch(0.08, warn.log = FALSE), log = NULL, cex = graphics::par("cex")) { 
        if(is.null(log)) {
            log <- paste(ifelse(graphics::par('xlog'), 'x', ''), ifelse(graphics::par('ylog'), 'y', ''), sep = '')
        }
        xlog <- 'x' %in% strsplit(log, NULL)[[1L]]
        ylog <- 'y' %in% strsplit(log, NULL)[[1L]]
        xy <- grDevices::xy.coords(x = x, y = y, recycle = TRUE, log = log)
        stopifnot((length(unique(xy$y)) <= 1))
        if(xlog) xy$x <- log10(xy$x)
        if(ylog) xy$y <- log10(xy$y)
        y.new <- xy$y + calculateSwarm(xy$x, dsize = xsize * cex, gsize = ysize * cex)
        out <- data.frame(x = x, y = y.new)
        if(ylog) out$y <- 10 ^ out$y
        out
    }
    
    floating.pie.asp <- function (xpos, ypos, x, edges = 200, radius = 1, col = NULL, startpos = 0, ...){
        u <- graphics::par("usr")
        user.asp <- diff(u[3:4])/diff(u[1:2])
        p <- graphics::par("pin")
        inches.asp <- p[2]/p[1]
        asp <- user.asp/inches.asp
        if (!is.numeric(x) || any(is.na(x) | x < 0)) 
            stop("floating.pie: x values must be non-negative")
        x <- c(0, cumsum(x)/sum(x))
        dx <- diff(x)
        nx <- length(dx)
        col <- if (is.null(col)) 
            grDevices::rainbow(nx)
        else rep(col, length.out = nx)
        if (length(i <- which(dx == 1))) {
            graphics::symbols(xpos, ypos, circles = radius, inches = FALSE, add = TRUE, bg = col[i], fg = col[which(dx == 0)])
        }else {
            bc <- 2 * pi * (x[1:nx] + dx/2) + startpos
            for (i in seq_len(nx)) {
                n <- max(2, floor(edges * dx[i]))
                t2p <- 2 * pi * seq(x[i], x[i + 1], length = n) + startpos
                xc <- c(cos(t2p) * radius + xpos, xpos)
                yc <- c(sin(t2p) * radius * asp + ypos, ypos)
                graphics::polygon(xc, yc, col = col[i], ...)
            }
        }
    }

    
    ############################################
    
    ## Calculate point positions g.pos, d.pos 
    if(method == 'swarm') {
        if(orientation=="horizontal") {
            g.offset <- lapply(x, function(a) swarmY(x = a, y = rep(0, length(a)), cex = sizeMultiplier)$y)
        } else {
            g.offset <- lapply(x, function(a) swarmX(x = rep(0, length(a)), y = a, cex = sizeMultiplier)$x)
        }
        d.pos <- x
    } else {
        if(method == 'hex') size.d <- size.d * sqrt(3) / 2

        if(log) {
            if(is.null(breaks)){
                breaks <- 10 ^ seq(log10(dlim[1]), log10(dlim[2]) + size.d, by = size.d)
            }
            if(length(breaks) == 1 && is.na(breaks[1])) {
                d.index <- x
                d.pos <- x
            } else {
                mids <- 10 ^ ((log10(head(breaks, -1)) + log10(tail(breaks, -1))) / 2)
                d.index <- lapply(x, cut, breaks = breaks, labels = FALSE)
                d.pos <- lapply(d.index, function(a) mids[a])
            }
        } else {
            if(is.null(breaks)){
                breaks <- seq(dlim[1], dlim[2] + size.d, by = size.d)
            }
            if(length(breaks) == 1 && is.na(breaks[1])) {
                d.index <- x
                d.pos <- x
            } else {
                mids <- (head(breaks, -1) + tail(breaks, -1)) / 2
                d.index <- lapply(x, cut, breaks = breaks, labels = FALSE)
                d.pos <- lapply(d.index, function(a) mids[a])
            }
        }

        x.index <- lapply(d.index, function(v) {
            if(length(stats::na.omit(v)) == 0) return(v)
            v.s <- lapply(split(v, v), seq_along)
            if(method == 'center')
                v.s <- lapply(v.s, function(a) a - mean(a))
            else if(method == 'square')
                v.s <- lapply(v.s, function(a) a - floor(mean(a)))
            else if(method == 'hex') {
                odd.row <- (as.numeric(names(v.s)) %% 2) == 1
                v.s[odd.row] <- lapply(v.s[odd.row], function(a) a - floor(mean(a)) - 0.25)
                v.s[!odd.row] <- lapply(v.s[!odd.row], function(a) a - ceiling(mean(a)) + 0.25)
            }
            unsplit(v.s, v)
        }) 

        g.offset <- lapply(1:n.groups, function(i) x.index[[i]] * size.g)
    }

    ## now check for runaway points (if "corral" has been set)
    if(corral != 'none') {
        if(missing(corralWidth)) {
            if(n.groups > 1) {
                corralWidth <- min(at[-1] - at[-n.groups]) - (2 * size.g)
            } else {
                corralWidth <- 2 * (min(diff(c(graphics::par('usr')[1], at, graphics::par('usr')[2]))) - size.g)
            }
        } else {
            stopifnot(length(corralWidth) == 1)
            stopifnot(corralWidth > 0)
        }
        halfCorralWidth <- corralWidth / 2
        if(corral == 'gutter') {
            g.offset <- lapply(g.offset, function(zz) pmin(halfCorralWidth, pmax(-halfCorralWidth, zz)))
        }
        if(corral == 'wrap') {
            g.offset <- lapply(g.offset, function(zz) ((zz + halfCorralWidth) %% (halfCorralWidth * 2)) - halfCorralWidth)
        }
        if(corral == 'random') {
            g.offset <- lapply(g.offset, function(zz) ifelse(zz > halfCorralWidth | zz < -halfCorralWidth, stats::runif(length(zz), -halfCorralWidth, halfCorralWidth), zz))
        }
        if(corral == 'omit') {
            g.offset <- lapply(g.offset, function(zz) ifelse(zz > halfCorralWidth, NA, ifelse(zz < -halfCorralWidth, NA, zz)))
        }
    }

    g.pos <- lapply(1:n.groups, function(i) at[i] + g.offset[[i]])

    out <- data.frame(x = unlist(g.pos), 
                      y = unlist(d.pos), 
                      pch = pch.out, 
                      col = col.out, 
                      bg = bg.out,
                      pie = pie.out,
                      x.orig = x.gp, 
                      y.orig = x.val,
                      stringsAsFactors = FALSE
                      )
    rownames(out) <- rnames.out
    
    ########################################################################################
    ## sort: first by x.orig and then y, col
    tmp <- data.frame(ind=1:nrow(out), x.orig=out$x.orig, y=out$y, col=out$col, pie=out$pie)
    ordering <- tmp[with(tmp, order(x.orig,y,col,pie)),]$ind
    ## make sure: for the same y, iterms are stayed together according to col and pie (proportion)
    a <- out[ordering, ]
    b <- split(a, a$x.orig)
    for(i in 1:length(b)){
        c <- b[[i]]
        #################
        d <- split(c, c$y)
        for(j in 1:length(d)){
            e <- d[[j]]
            e$x <- sort(e$x)
            d[[j]] <- e
        }
        #################
        ttmp <- d[[1]]
        if(length(d)>=2){
            for(k in 2:length(d)){
                ttmp <- rbind(ttmp, d[[k]])
            }
        }
        b[[i]] <- ttmp
    }
    out <- b[[1]]
    if(length(out)>=2){
        for(i in 2:length(b)){
            out <- rbind(out, b[[i]])
        }
    }
    ## back to the original order as the input data
    ind <- match(rownames(data), rownames(out))
    out <- out[ind,]
    ########################################################################################

    if(do.plot) {
        if(orientation=="horizontal") {
            
            if(is.numeric(out$pch)){
                graphics::points(out$y, out$x, pch=out$pch, col=out$col, bg=out$bg, cex=cex)
            }else{
                for(i in 1:length(out$pch)){
                    pch_tmp <- out$pch[i]
                    if(pch_tmp == 'circles'){
                        graphics::symbols(out$y[i], out$x[i], circles=min(size.g, size.d)/2, inches=F, fg=out$col[i], bg=out$bg[i], add=T)
                    }else if(pch_tmp == 'thermometers'){
                        graphics::symbols(out$y[i], out$x[i], thermometers=cbind(0.4*size.d, 0.8*size.g, out$pie[i]), inches=F, fg=out$col[i], bg=out$bg[i], add=T)
                    }else if(pch_tmp == 'pies'){
                        floating.pie.asp(xpos=out$y[i], ypos=out$x[i], x=cbind(out$pie[i], 1-out$pie[i]), edges=200, radius=0.8*max(size.g, size.d)/2, col=cbind(out$col[i], "transparent"), startpos=0, border=out$col[i])
                    }else{
                        graphics::points(out$x[i], out$y[i], pch=out$pch[i], col=out$col[i], bg=out$bg[i], cex=cex)
                    }
                }
            }
            
            if(!add) {
                graphics::axis(1, ...)
                graphics::axis(2, at=at, labels=labels, tick=FALSE, ...)
                graphics::box(...)
            }
            ## add boxplot
            if(do.boxplot){
                graphics::boxplot(formula=formula, data=data, at=at, names=labels, outline=F, add=T, horizontal=T, notch=boxplot.notch, border=boxplot.border, col=boxplot.col)
            }
            
        }else if(orientation=="vertical") {
            
            if(is.numeric(out$pch)){
                graphics::points(out$x, out$y, pch=out$pch, col=out$col, bg=out$bg, cex=cex)
            }else{
                for(i in 1:length(out$pch)){
                    pch_tmp <- out$pch[i]
                    if(pch_tmp == 'circles'){
                        graphics::symbols(out$x[i], out$y[i], circles=min(size.g, size.d)/2, inches=F, fg=out$col[i], bg=out$bg[i], add=T)
                    }else if(pch_tmp == 'thermometers'){
                        graphics::symbols(out$x[i], out$y[i], thermometers=cbind(0.4*size.g, 0.8*size.d, out$pie[i]), inches=F, fg=out$col[i], bg=out$bg[i], add=T)
                    }else if(pch_tmp == 'pies'){
                        floating.pie.asp(xpos=out$x[i], ypos=out$y[i], x=cbind(out$pie[i], 1-out$pie[i]), edges=200, radius=0.8*min(size.g, size.d)/2, col=cbind(out$col[i], "transparent"), startpos=0, border=out$col[i])
                    }else{
                        graphics::points(out$x[i], out$y[i], pch=out$pch[i], col=out$col[i], bg=out$bg[i], cex=cex)
                    }
                }
            }
            
            if(!add) {
                graphics::axis(2, ...)
                graphics::axis(1, at=at, labels=labels, tick=FALSE, ...)
                graphics::box(...)
            }
            ## add boxplot
            if(do.boxplot){
                graphics::boxplot(formula=formula, data=data, at=at, names=labels, outline=F, add=T, horizontal=F, notch=boxplot.notch, border=boxplot.border, col=boxplot.col)
            }
        }
    }
    
    invisible(out)
}