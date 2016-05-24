splitcalls <- function(pd) {
    stacks <- pd$stacks
    counts <- pd$counts
    first <- sapply(stacks, `[`, 1)
    spstacks <- split(stacks, first)
    spcounts <- split(counts, first)
    totals <- lapply(spcounts, sum)  ## maybe sapply?
    merge <- function(total, stacks, counts) {
        keep <- sapply(stacks, function(s) length(s) > 1)
        stacks <- lapply(stacks[keep], `[`, -1)
        counts <- counts[keep]
        list(total = total, stacks = stacks, counts = counts)
    }
    mapply(merge, totals, spstacks, spcounts, SIMPLIFY=FALSE)
}

makeTree <- function(fun, pd) {
    if (length(pd$stacks) > 0) {
        spcalls <- splitcalls(pd)
        funs <- names(spcalls)
        calls <- lapply(seq_along(spcalls),
                        function(i) makeTree(funs[i], spcalls[[i]]))
    }
    else calls <- NULL
    list(fun = fun, hits = pd$total, calls = calls)
}

splitRect <- function(score, left, bottom, right, top) {
    stopifnot(all(score > 0))
    cumscore <- cumsum(score) / sum(score)
    nc <- length(cumscore)
    width <- right - left
    height <- top - bottom
    if (width < height) {
        bottoms <- bottom + height * c(0, cumscore[-nc])
        tops <- bottom + height * cumscore
        lefts <- rep(left, nc)
        rights <- rep(right, nc)
    }
    else {
        bottoms <-rep(bottom, nc)
        tops <- rep(top, nc)
        lefts <- left + width * c(0, cumscore[-nc])
        rights <- left + width * cumscore
    }
    list(left = lefts, bottom = bottoms, right = rights, top = tops)
}

squarifiedTiles <- function(v, left, bottom, right, top) {
    squarify <- function(children, row,  w) {
        child <- children[1] ## **** check for length 0?
        newRow <- c(row, child)
        if (length(row) == 0 || worst(row, w) >= worst(newRow, w)) {
            remainingChildren <- children[-1]
            if (length(remainingChildren) == 0)
                layoutrow(newRow)
            else
                squarify(remainingChildren, newRow, w)
        }
        else {
            layoutrow(row)
            squarify(children, numeric(0), width())
        }
    }

    width <- function() min(right - left, top - bottom)

    layoutrow <- function(r) {
        if (right - left >= top - bottom) {
            a <- sum(r) / (top - bottom)
            b <- r / a
            rects <<- rbind(rects,
                            cbind(left,
                                  c(bottom, bottom + cumsum(b)[-length(r)]),
                                  left + a,
                                  bottom + cumsum(b)))
            left <<- left + a
        }
        else {
            b <- sum(r) / (right - left)
            a <- r / b
            rects <<- rbind(rects,
                            cbind(c(left, left + cumsum(a)[-length(r)]),
                                  bottom,
                                  left + cumsum(a),
                                  bottom + b))
            bottom <<- bottom + b
        }
    }

    worst <- function(R, w) {
        s <- sum(R)
        max(w^2 * R / s^2, s^2 / (w^2 * R))
    }

    stopifnot(all(v > 0) && all(diff(v) <= 0))

    v <- (v / sum(v)) * (right - left) * (top - bottom)

    rects <- NULL
    squarify(v, numeric(0), min(right - left, top - bottom))
    list(left = rects[,1], bottom = rects[,2],
         right = rects[,3], top = rects[,4])
}

makeTreeMapData <- function(n, left, bottom, right, top,
                            cex = 0.75, depth = 1, tile = splitRect) {
    hits <- n$hits
    calls <- n$calls
    nc <- length(calls)
    callhits <- if (nc > 0) sapply(calls, function(x) x$hits) / hits else 0
    selffrac <- 1 - sum(callhits)
    lpad <- strwidth("M", cex = cex)* 0.6 ## was 0.3
    label <- n$fun
    wd <- strwidth(label, cex = cex) + 2 * lpad
    ht <- strheight(label, cex = cex) + 2 * lpad
    if (wd <= right - left && ht <= selffrac * (top - bottom)) {
        showLabel <- TRUE
        rotate <- FALSE
        labX <- left + lpad
        labY <- top - lpad
    }
    else if (wd <= top - bottom && ht <= selffrac * (right - left)) {
        showLabel <- TRUE
        rotate <- TRUE
        labX <- right - lpad
        labY <- top - lpad
    }
    else {
        showLabel <- FALSE
        rotate <- FALSE
        labX <- 0
        labY <- 0
    }
    if (nc > 0) {
        sRight <- right
        sTop <- top
        if (showLabel)
            if (rotate) sRight <- sRight - ht
            else sTop <- sTop - ht
        w <- (sRight - left) * sqrt(1 - selffrac)
        h <- (sTop - bottom) * sqrt(1 - selffrac)
        dw <- ((sRight - left) - w) / 2
        dh <- ((sTop - bottom) - h) / 2
        sLeft <- left + (if (rotate) min(dh, dw) else dw)
        sRight <- sLeft + w
        sBottom <- bottom + (if (rotate) dh else min(dh, dw))
        sTop <- sBottom + h
        calls <- calls[order(callhits, decreasing = TRUE)]
        callhits <- sort(callhits, decreasing = TRUE)
        s <- tile(callhits, sLeft, sBottom, sRight, sTop)
        rest <- mapply(function(n, left, bottom, right, top)
                       makeTreeMapData(n, left, bottom, right, top,
                                       cex, depth + 1, tile),
                       calls, s$left, s$bottom, s$right, s$top,
                       SIMPLIFY = FALSE)
    }
    else rest <- NULL
    d <- data.frame(left = left, bottom = bottom, right = right, top = top,
                    hits = hits, depth = depth,
                    label = label, labX = labX, labY = labY,
                    showLabel = showLabel, rotate = rotate,
                    stringsAsFactors = FALSE)
    if (is.null(rest)) d
    else rbind(d, do.call(rbind, rest))
}

calleeTreeMap <- function(pd, srclines = FALSE, cex = 0.75, colormap = NULL,
                          main = "Callee Tree Map", squarify = FALSE,
                          border = NULL) {
    pd$total <- sum(pd$counts)
    plot(c(0,1), c(0,1), type = "n", xlab = "", ylab = "",
         axes = FALSE, main = main)
    if (srclines)
        pd$stacks <- refStacks(pd)
    tile <- if (squarify) squarifiedTiles else splitRect
    v <-makeTreeMapData(makeTree("", pd), 0, 0, 1, 1, cex = cex, tile = tile)
    nc <- nrow(v)
    cmap <- if (! is.null(colormap)) colormap else default.cmap
    fun <- sub(" .*", "", v$label)
    col <- cmap(fun, v$depth, v$hits)
    rect(v$left, v$bottom, v$right, v$top, col = col, border = border)
    vlr <- v[v$showLabel & v$rotate,]
    if (nrow(vlr) > 0)
        with(vlr, text(labX, labY, label, srt = -90, adj = c(0, 1), cex = cex))
    vlnr <- v[v$showLabel & ! v$rotate,]
    if (nrow(vlnr) > 0)
        with(vlnr, text(labX, labY, label, adj = c(0, 1), cex = cex))
    invisible(structure(v, class = c("proftools_calleeTreeMap", class(v))))
}

ctmIdentify <- function(x, n = 1, print = FALSE, ...) {
    val <- list()
    while (n > 0) {
        n <- n - 1
        loc <- locator(1)
        if (! is.null(loc)) {
            idx <- which(loc$x >= x$left & loc$x <= x$right &
                         loc$y >= x$bottom & loc$y <= x$top)
            if (length(idx) > 0) {
                stack <- x$label[idx][-1]
                if (print)
                    cat(stack, "\n")
                val <- c(val, list(stack))
            }
            else break
        }
        else break
    }
    val
}

identify.proftools_calleeTreeMap <- ctmIdentify
