graph_cred_int <-
function(hex,title='Credibility \n level',conf.int,count)
{
  # We change the function to obtain the good title
  hex$legend$right$fun <- function (legend = 1.2, inner = legend/5, cex.labels = 1, cex.title = 1.2,
    style = "colorscale", minarea = 0.05, maxarea = 0.8, mincnt = 1,
    maxcnt, trans = NULL, inv = NULL, colorcut = seq(0, 1, length = 17),
    density = NULL, border = NULL, pen = NULL, colramp = function(n) {
        LinGray(n, beg = 90, end = 15)
    }, ..., vp = NULL, draw = FALSE)
{
    style <- match.arg(style, eval(formals(grid.hexagons)[["style"]]))
    if (style %in% c("centroids", "lattice", "colorscale")) {
        if (is.null(trans)) {
            sc <- maxcnt - mincnt
            bnds <- round(mincnt + sc * colorcut)
        }
        else {
            if (!is.function(trans) && !is.function(inv))
                stop("'trans' and 'inv' must both be functions if 'trans' is not NULL")
            con <- trans(mincnt)
            sc <- trans(maxcnt) - con
            bnds <- round(inv(con + sc * colorcut))
        }
    }
    ans <- switch(style, colorscale = {
        n <- length(bnds)
        pen <- colramp(n - 1)
        hexxy <- hexcoords(dx = 1, n = 1)[c("x", "y")]
        maxxy <- max(abs(unlist(hexxy)))
        hexxy <- lapply(hexxy, function(x) 0.5 * x/maxxy)
        pol <- polygonGrob(x = 0.5 + rep(hexxy$x, n - 1), y = (rep(1:(n -
            1), each = 6) + hexxy$y)/n, id.lengths = rep(6, n -
            1), gp = gpar(fill = pen, col = border), default.units = "npc")
        txt <- textGrob(as.character(bnds), x = 0.5, y = (0:(n -
            1) + 0.5)/n, gp = gpar(cex = cex.labels), default.units = "npc")
        ttl <- textGrob(title, gp = gpar(cex = cex.title))
        key.layout <- grid.layout(nrow = 2, ncol = 2, heights = unit(c(1.5,
            1), c("grobheight", "grobheight"), data = list(ttl,
            txt)), widths = unit(c(1/n, 1), c("grobheight", "grobwidth"),
            data = list(pol, txt)), respect = TRUE)
        key.gf <- frameGrob(layout = key.layout, vp = vp)
        key.gf <- placeGrob(key.gf, ttl, row = 1, col = 1:2)
        key.gf <- placeGrob(key.gf, pol, row = 2, col = 1)
        key.gf <- placeGrob(key.gf, txt, row = 2, col = 2)
        key.gf
    }, centroids = , lattice = {
        warning("legend shows relative sizes")
        radius <- sqrt(minarea + (maxarea - minarea) * colorcut)
        n <- length(radius)
        if (is.null(pen)) pen <- 1
        if (is.null(border)) border <- pen
        hexxy <- hexcoords(dx = 1, n = 1)[c("x", "y")]
        maxxy <- max(abs(unlist(hexxy)))
        hexxy <- lapply(hexxy, function(x) 0.5 * x/maxxy)
        pol <- polygonGrob(x = 0.5 + rep(radius, each = 6) *
            rep(hexxy$x, n), y = (rep(0.5 + 1:n, each = 6) +
            rep(radius, each = 6) * hexxy$y - 1)/n, id.lengths = rep(6,
            n), gp = gpar(fill = pen, col = border), default.units = "npc")
        txt <- textGrob(as.character(bnds), x = 0.5, y = (1:n -
            0.5)/n, gp = gpar(cex = cex.labels), default.units = "npc")
        ttl <- textGrob(title, gp = gpar(cex = cex.title))
        key.layout <- grid.layout(nrow = 2, ncol = 2, heights = unit(c(1.5,
            1), c("grobheight", "grobheight"), data = list(ttl,
            txt)), widths = unit(c(1/n, 1), c("grobheight", "grobwidth"),
            data = list(pol, txt)), respect = TRUE)
        key.gf <- frameGrob(layout = key.layout, vp = vp)
        key.gf <- placeGrob(key.gf, ttl, row = 1, col = 1:2)
        key.gf <- placeGrob(key.gf, pol, row = 2, col = 1)
        key.gf <- placeGrob(key.gf, txt, row = 2, col = 2)
        key.gf
    }, nested.lattice = , nested.centroids = {
        dx <- inner/2
        dy <- dx/sqrt(3)
        hexC <- hexcoords(dx, dy, n = 1, sep = NULL)
        numb <- cut(floor(legend/inner), breaks = c(-1, 0, 2,
            4))
        if (is.na(numb)) numb <- 4
        switch(numb, {
            warning("not enough space for legend")
            return(textGrob(""))
        }, size <- 5, size <- c(1, 5, 9), size <- c(1, 3, 5,
            7, 9))
        xmax <- length(size)
        radius <- sqrt(minarea + (maxarea - minarea) * (size -
            1)/9)
        txt <- as.character(size)
        lab <- c("Ones", "Tens", "Hundreds", "Thousands", "10 Thousands",
            "100 Thousands", "Millions", "10 Millions", "100 Millions",
            "Billions")
        power <- floor(log10(maxcnt)) + 1
        yinc <- 16 * dy
        ysize <- yinc * power
        xmid <- 0
        x <- inner * (1:xmax - (1 + xmax)/2) + xmid
        n <- length(x)
        tx <- rep.int(hexC$x, n)
        ty <- rep.int(hexC$y, n)
        six <- rep.int(6:6, n)
        y <- rep.int(3 * dy - 0.75 * yinc, xmax)
        if (is.null(pen)) {
            pen <- 1:power + 1
            pen <- cbind(pen, pen + 10)
        }
        if (is.null(border)) border <- TRUE
        key.layout <- grid.layout(nrow = 1, ncol = 1, heights = unit(ysize,
            "inches"), widths = unit(legend, "inches"), respect = TRUE)
        key.gf <- frameGrob(layout = key.layout, vp = vp)
        n6 <- rep.int(6, n)
        for (i in 1:power) {
            y <- y + yinc
            key.gf <- placeGrob(key.gf, polygonGrob(x = unit(legend/2 +
                rep.int(hexC$x, n) + rep.int(x, n6), "inches"),
                y = unit(rep.int(hexC$y, n) + rep.int(y, n6),
                  "inches"), id.lengths = n6, gp = gpar(col = pen[i,
                  1], fill = if (border) 1 else pen[i, 1])),
                row = 1, col = 1)
            key.gf <- placeGrob(key.gf, polygonGrob(x = legend/2 +
                tx * rep.int(radius, six) + rep.int(x, six),
                y = ty * rep.int(radius, six) + rep.int(y, six),
                default.units = "inches", id = NULL, id.lengths = rep(6,
                  n), gp = gpar(fill = pen[i, 2], col = border)),
                row = 1, col = 1)
            key.gf <- placeGrob(key.gf, textGrob(txt, x = legend/2 +
                x, y = y - 4.5 * dy, default.units = "inches",
                gp = gpar(cex = cex.labels)), row = 1, col = 1)
            key.gf <- placeGrob(key.gf, textGrob(lab[i], x = legend/2 +
                xmid, y = y[1] + 4.5 * dy, default.units = "inches",
                gp = gpar(cex = 1.3 * cex.title)), row = 1, col = 1)
        }
        key.gf
    })
    if (draw) {
        grid.draw(ans)
        invisible(ans)
    }
    else ans
  }
  
  colorcut <- hex$panel.args.common$colorcut*max(count)
  
  perc <- rep(NA,length(colorcut))
  for ( i in 1:length(colorcut) )
  {
    perc[i] <- sum( count[ count >= colorcut[i] ] )/( sum( count ) )
  }
  hex$legend$right$args$colorcut <- unique(perc)
  
  hex$legend$right$args$maxcnt <- 100
  hex$legend$right$args$mincnt <- 0
  
  return(hex)
}