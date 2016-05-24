############################################################################################
## package 'secrlinear'
## utility.R
## last changed 2014-08-28 2014-09-09 2014-10-31
## 2014-11-07 allow breaks in lines
############################################################################################

replacedefaults <- function (default, user) replace(default, names(user), user)

along.line <- function (route, spacing, break.at = 0) {
    ## Interpolate points along route at interval 'spacing'
    ## route is matrix of vertices
    np <- nrow(route)
    if (nrow(route)<2) {
        return (matrix(nrow=0, ncol=2))
    }
    dxy <- route[2:np,,drop=FALSE] - route[1:(np-1),,drop=FALSE]
    d <- apply(dxy^2,1,sum)^0.5
    d[break.at] <- 0    ## disjoint lines
    cumd <- cumsum(d)
    linelength <- cumd[np-1]
    if (linelength == 0)
        return (matrix(nrow=0, ncol=2))
    else {

                ipos <- seq(spacing/2, linelength, spacing)
        first <- sapply(ipos, function(x) match(TRUE, cumd>=x))
        cumd <- c(0,cumd)
        add <- ipos - cumd[first]
        pr <- add / d[first]
        result <- route[first,, drop = FALSE] + dxy[first,] * pr
        result <- cbind(result, rep(0, nrow(result)), ipos)
        itermini <- c(cbind(break.at, break.at+1)[-1],np)
        termini <- cbind(route[itermini,], rep(1, length(itermini)), cumd[itermini])
        result <- rbind(result, termini)
        result <- result[order(result[,4]), 1:3]
        return(result)
    }
}
#------------------------------------------------------------------------------------------
## r <- matrix(runif(10), ncol=2)
## plot(r, type = 'o')
## tmp <- alongline(r, 0.1)
## points(tmp, col='red')
#------------------------------------------------------------------------------------------

sample.line <- function(x, spacing) {
    ## Function for systematic point sample
    ## Input  -- SpatialLinesDataFrame
    ## Output -- SpatialPointsDataFrame
    ## default type=2 includes endpoints
    sampleone <- function (i) {
        lsub <- x[i,]
        ns <- trunc(lgth[i] / spacing)
        if (ns>0) {
            xyL <- coordinates(lsub)[[1]]
            nLine <- length(xyL)
            nperLine <- sapply(xyL, nrow)
            xy <- do.call(rbind, xyL)
            breaks <- c(0, cumsum(nperLine)[-nLine])
            lsamp <- along.line(xy, spacing, breaks)
            nrl <- nrow(lsamp)
            pdf <- data.frame(LineID = rep(i, nrl), Terminal = as.logical(lsamp[,3]))
            pdf <- merge(pdf, lsub)
            SpatialPointsDataFrame(lsamp[,1:2], data = pdf)
        }
        else NULL
    }
    if (!is(x, "SpatialLinesDataFrame")) stop("x must be SpatialLinesDataFrame")
    lgth <- SpatialLinesLengths(x)
    nr <- nrow(x)
    ## For each feature
    results <- sapply(1:nr, sampleone)
    results <- results[!sapply(results, is.null)]
    do.call(rbind, results)
}

getCentres <- function (xy) {
    ## redundant
    nrxy <- nrow(xy)
    if (nrxy > 1) {
        xy <- (xy[-1,] + xy[-nrxy,]) / 2
        rownames(xy) <- 1:(nrxy-1)
    }
    xy
}
#------------------------------------------------------------------------------------------

make.sldf <- function (coord, f) {
    ## Form SpatialLinesDataFrame from coordinates data
    ## coord is dataframe of coordinates - must include columns 'x','y'
    ## f is vector of values by which to split coord rows
    if (missing(f) & ('LineID' %in% names(coord)))
        f <-  coord$LineID
    if (missing(f))
        coordlist <- list(coord)
    else
        coordlist <- split(coord[,c('x','y')], f)
    S0 <- lapply(coordlist, Line)
    S1 <- Lines(S0, ID = '1')
    S2 <- SpatialLines(list(S1))
    ldf <-  data.frame(ID = 1:length(S2), rownames=1)
    SpatialLinesDataFrame(S2, data = ldf)
}
#------------------------------------------------------------------------------------------

snapPointsToLinearMask <- function (xy, mask) {
    gr <- attr(mask, 'graph')
    geometryxy <- data.frame(x = V(gr)$x, y = V(gr)$y)
    matchedxy <- nearesttrap(xy, geometryxy)
    tempxy <- geometryxy[matchedxy,,drop = FALSE]
    mostattributes(tempxy) <- attributes(xy)
    tempxy
}

#------------------------------------------------------------------------------------------

showedges <- function (mask, plt = TRUE, add = FALSE, type = c('all', 'noskips', 'skips'),
                       lengths = c(0,Inf), ...) {
    type <- match.arg(type)
    gr <- attr(mask, 'graph')
    if (!is.igraph(gr))
        stop("mask does not have valid graph attribute")
    ed <- get.edges(gr, E(gr))
    x <- V(gr)$x
    y <- V(gr)$y
    LineID <- V(gr)$LineID
    OK <- (E(gr)$weight >= lengths[1]) & (E(gr)$weight <= lengths[2])
    if (type != 'all') {
        mixed <- V(gr)$LineID[ed[,1]] != V(gr)$LineID[ed[,2]]
        if (type == 'skips')
            OK <- OK & mixed
        else
            OK <- OK & !mixed
    }
    start <- ed[OK,1]
    finish <- ed[OK,2]
    LineID1 <- LineID[start]
    LineID2 <- LineID[finish]
    x1 <- x[start]
    y1 <- y[start]
    x2 <- x[finish]
    y2 <- y[finish]
    ed <- data.frame(start = start, finish = finish,
                     line1 = LineID1, line2 = LineID2,
                     x1 = x1, y1 = y1, x2 = x2, y2 = y2,
                     weight = E(gr)$weight[OK])
    if (plt) {
        if (!add) plot(mask)
        segments(x1, y1, x2, y2, ... )
        invisible (ed)
    }
    else
        ed
}
#------------------------------------------------------------------------------------------
getLineID <- function (mask, laboffset= rep(spacing(mask)*3,2), ...) {
    if (is.null(covariates(mask)$LineID))
        stop("LineID not found in covariates(mask)")
    plot(mask, ...)
    cat ("click on line \n")
    output <- data.frame(Point=numeric(0), LineID=character(0))
    repeat {
        xy1 <- as.data.frame(locator(1))
        if (nrow(xy1) < 1)
            break
        else {
            matched <- nearesttrap(xy1, mask)
            lineID <- as.character(covariates(mask)$LineID[matched])
            cat("Point ", matched, " is on line ", lineID, "\n")
            points(mask$x[matched], mask$y[matched], pch=1)
            text (mask$x[matched] + laboffset[1], mask$y[matched] + laboffset[2],
                  lineID, cex=0.7, col = 'red')
            output<- rbind(output, data.frame(Point=matched, LineID=lineID, stringsAsFactors=FALSE))
        }
    }
    invisible (output)
}

#------------------------------------------------------------------------------------------
replot <- function(mask, xlim=NULL, ylim=NULL, ...) {
    if (is.null(xlim) | is.null(ylim)) {
        cat ("click on opposite corners \n")
        bb <- as.data.frame(locator(2))
        if (nrow(bb)<2) return(NULL)
        xlim <- sort(bb[,1])
        ylim <- sort(bb[,2])
    }
    eqscplot(0,0, xlim=xlim, ylim=ylim, xlab='', ylab='', axes=F, type='n' )
    plot(mask, add = TRUE, ...)
    NULL
}
#------------------------------------------------------------------------------------------

deleteedges <- function(mask, replot = TRUE, ...) {
    gr <- attr(mask, 'graph')
    x <- V(gr)$x
    y <- V(gr)$y
    usr <- par()$usr
    cat ("click on both ends \n")
    repeat {
        xy <- as.data.frame(locator(2))
        if (nrow(xy) < 2)
            break
        else {
            matched <- nearesttrap(xy, cbind(x,y))
            matched <- sort(matched)
            gr[matched[1], matched[2]] <- NULL
            cat("Edge between ", matched[1], " and ", matched[2], " deleted\n")
            segments(x[matched[1]], y[matched[1]], x[matched[2]], y[matched[2]], col='white', lwd=3)
        }
    }
    attr(mask, 'graph') <- gr
    if (replot) {
        replot (mask, xlim = usr[1:2], ylim = usr[3:4])
        showedges (mask, ...)
    }
    invisible(mask)
}
#------------------------------------------------------------------------------------------

addedges <- function(mask, replot = TRUE, ...) {
    gr <- attr(mask, 'graph')
    x <- V(gr)$x
    y <- V(gr)$y
    grxy <- cbind(x,y)
    usr <- par()$usr
    cat ("click on both ends \n")
    repeat {
        xy <- as.data.frame(locator(2))
        if (nrow(xy) < 2)
            break
        else {
            matched <- nearesttrap(xy, grxy)
            matched <- sort(matched)
            gr[matched[1], matched[2], attr="weight"] <- edist(grxy[matched[1],,drop=F],
                                                               grxy[matched[2],,drop=F])
            cat("Edge added between ", matched[1], " and ", matched[2], "\n")
            segments(x[matched[1]], y[matched[1]], x[matched[2]], y[matched[2]], col='green', lwd=3)
        }
    }
    attr(mask, 'graph') <- gr
    if (replot) {
        replot (mask, xlim = usr[1:2], ylim = usr[3:4])
        showedges (mask, ...)
    }
    invisible(mask)
}
#------------------------------------------------------------------------------------------

cleanskips <- function (mask) {
    skips <- showedges(mask, type = 'skips', plt = FALSE)
    if (nrow(skips) > 0) {
        skips <- split(skips, paste(skips$line1, skips$line2, sep = '.'))
        getbadedges <-  function(x) {
            notOK <- (1:nrow(x)) != which.min(x$weight)  ## not shortest
            if (any(notOK))
                as.matrix(x[notOK, c('start','finish'), drop = FALSE])
            else
                matrix(nrow=0,ncol=2)
        }
        badedge <- lapply(skips, getbadedges)
        badedges <- do.call(rbind, badedge)
        if (nrow(badedges)>0) {
            gr <- attr(mask, 'graph')
            gr[from = badedges[,1], to = badedges[,2]] <- NULL
            attr(mask, 'graph') <- gr
        }
    }
    mask
}
#------------------------------------------------------------------------------------------

cleanskipsgraph <- function (gr) {
    ed <- get.edges(gr, E(gr))
    LineID <- V(gr)$LineID
    OK <- LineID[ed[,1]] != LineID[ed[,2]]
    start <- ed[OK,1]
    finish <- ed[OK,2]
    weight <- E(gr)$weight[OK]
    skips <- data.frame(start = start, finish = finish,
                     line1 = LineID[start], line2 = LineID[finish],
                     weight = weight)
    if (nrow(skips) > 0) {
        skips <- split(skips, paste(skips$line1,skips$line2,sep='.'))
        getbadedges <-  function(x) {
            notOK <- (1:nrow(x)) != which.min(x$weight)  ## not shortest
            if (any(notOK))
                as.matrix(x[notOK, c('start','finish'), drop = FALSE])
            else
                matrix(nrow=0,ncol=2)
        }
        badedge <- lapply(skips, getbadedges)
        badedges <- do.call(rbind, badedge)
        gr[from = badedges[,1], to = badedges[,2]] <- NULL
    }
    gr
}
#------------------------------------------------------------------------------------------
## suppose we have a SpatialLines object and want a SpatialLinesDataFrame
## this should do it
## df <- data.frame(z = c(1,2), row.names= sapply(slot(Sl, "lines"), function(x) slot(x, "ID")))
## Sldf <- SpatialLinesDataFrame(Sl, data = df)
