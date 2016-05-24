
# C_text(xy.coords, labels, adj, pos, offset, vfont, cex, col, font, ...)

C_text <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:10)])
    dev.set(playDev())
    # TODO:  handle 'pos', 'offset', 'vfont'
    depth <- gotovp(par$xpd)
    xx <- tx(x[[2]]$x, par)
    yy <- ty(x[[2]]$y, par)
    labels <- x[[3]]
    adj <- x[[4]]
    just <- just(adj, par)
    adjx <- just[1]
    adjy <- just[2]
    pos <- x[[5]]
    offset <- unit(x[[6]]*par$cin[2]*par$cex, "in")
    if (!is.null(pos)) {
        n <- length(labels)
        pos <- rep(pos, length.out=n)
        adjx <- rep(0.5, length.out=n)
        # 0.3333 = yCharOffset
        adjy <- rep(0.3333, length.out=n)
        xx <- rep(xx, length.out=n)
        yy <- rep(yy, length.out=n)
        xx[pos == 2] <- xx - convertWidth(offset, "native", valueOnly=TRUE)
        xx[pos == 4] <- xx + convertWidth(offset, "native", valueOnly=TRUE)
        yy[pos == 1] <- yy - convertHeight(offset, "native", valueOnly=TRUE)
        yy[pos == 3] <- yy + convertHeight(offset, "native", valueOnly=TRUE)
        adjx[pos == 2] <- 1
        adjx[pos == 4] <- 0   
        adjy[pos == 1] <- 1 - (0.5 - 0.3333) 
        adjy[pos == 3] <- 0
    }
    vfont <- x[[7]]
    cex <- FixupCex(x[[8]]*par$cex, 1)
    cex <- ifelse(is.na(cex), par$cex, cex)
    col <- FixupCol(x[[9]], NA, par$bg)
    col <- ifelse(is.na(col), par$col, col)
    font <- FixupFont(x[[10]], NA)
    font <- ifelse(is.na(font), par$font, font)
    family <- par$family
    if (!is.null(vfont) && !is.language(labels)) {
        # Override 'font' and 'family'
        font <- vfont[2]
        family <- mapVfont(vfont[1])
    }
    grid.text(labels, xx, yy, default.units="native",
              hjust=adjx, vjust=adjy, rot=par$srt,
              gp=gpar(cex=cex, col=col, fontface=font, fontfamily=family,
                  lineheight=par$lheight),
              name=grobname("text"))
    upViewport(depth)
}

just <- function(adj, par) {
    if (is.null(adj) || length(adj) == 0) {
        adjx <- par$adj
        adjy <- NA
    } else {
        if (length(adj) == 1) {
            adjx <- adj
            adjy <- NA
        } else {
            adjx <- adj[1]
            adjy <- adj[2]
        }
    }
    c(adjx, adjy)
}

mapVfont <- function(vfont) {
    switch(vfont,
           "serif"="HersheySerif",
           "sans serif"="HersheySans",
           "script"="HersheyScript",
           "gothic english"="HersheyGothicEnglish",
           "gothic german"="HersheyGothicGerman",
           "gothic italian"="HersheyGothicItalian",
           "serif symbol"="HersheySymbol",
           "sans serif symbol"="HersheySansSymbol")
}
