# This function is a slightly modified version of plot.ca() from package ca.
# Released under the GPL (no version specified), Copyright Michael Greenacre
# and Oleg Nenadic <onenadi at uni-goettingen.de>.
# http://cran.r-project.org/web/packages/ca/index.html
plotCorpusCa <- function (x, dim = c(1, 2), map = "symmetric",
                          what = c("all", "all"), mass = c(FALSE, FALSE),
                          contrib = c("none", "none"), col = c("blue", "red"),
                          col.text = c("black", "blue", "black", "red"),
                          font = c(3, 4, 1, 2),
                          pch = c(16, 1, 17, 24), labels = c(2, 2),
                          arrows = c(FALSE, FALSE), cex = 0.75,
                          xlab = paste("Dimension", dim[1]),
                          ylab = paste("Dimension", dim[2]), ...)
{
    obj <- x
    if (length(what) != 2) 
        what <- rep(what, length = 2)
    if (length(mass) != 2) 
        mass <- rep(mass, length = 2)
    if (length(contrib) != 2) 
        contrib <- rep(contrib, length = 2)
    if (length(col) != 2) 
        col <- rep(col, length = 2)
    if (length(col.text) != 2) 
        col.text <- rep(col.text, length = 4)
    if (length(font) != 2) 
        font <- rep(font, length = 4)
    if (length(labels) != 2) 
        labels <- rep(labels, length = 2)
    if (length(pch) != 4) 
        pch <- rep(pch, length = 4)
    if (!is.numeric(x$suprow)) {
        if (map == "colgab" | map == "colgreen") {
            if (what[1] != "none") 
                what[1] <- "active"
        }
    }
    if (!is.numeric(x$supcol)) {
        if (map == "rowgab" | map == "rowgreen") {
            if (what[2] != "none") 
                what[2] <- "active"
        }
    }
    if (min(dim) < 0) {
        swisign <- ifelse(dim < 0, -1, 1)
        dim.c <- dim(obj$rowcoord)[2]
        signmat <- diag(rep(swisign, length = dim.c))
        obj$rowcoord <- obj$rowcoord %*% signmat
        obj$colcoord <- obj$colcoord %*% signmat
        dim <- abs(dim)
    }
    K <- dim(obj$rowcoord)[2]
    I <- dim(obj$rowcoord)[1]
    J <- dim(obj$colcoord)[1]
    svF <- matrix(rep(obj$sv[1:K], I), I, K, byrow = TRUE)
    svG <- matrix(rep(obj$sv[1:K], J), J, K, byrow = TRUE)
    rpc <- obj$rowcoord * svF
    cpc <- obj$colcoord * svG
    symrpc <- obj$rowcoord * sqrt(svF)
    symcpc <- obj$colcoord * sqrt(svG)
    mt <- c("symmetric", "rowprincipal", "colprincipal", "symbiplot", 
        "rowgab", "colgab", "rowgreen", "colgreen")
    mti <- 1:length(mt)
    mtlut <- list(symmetric = list(x = rpc, y = cpc), rowprincipal = list(x = rpc, 
        y = obj$colcoord), colprincipal = list(x = obj$rowcoord, 
        y = cpc), symbiplot = list(x = symrpc, y = symcpc), rowgab = list(x = rpc, 
        y = obj$colcoord * obj$colmass), colgab = list(x = obj$rowcoord * 
        obj$rowmass, y = cpc), rowgreen = list(x = rpc, y = obj$colcoord * 
        sqrt(obj$colmass)), rowgreen = list(x = obj$rowcoord * 
        sqrt(obj$rowmass), y = cpc))
    x <- mtlut[[mti[mt == map]]][[1]]
    y <- mtlut[[mti[mt == map]]][[2]]
    x.names <- obj$rownames
    y.names <- obj$colnames
    rm(mt, mti, mtlut)
    indx <- dim(x)[1]
    indy <- dim(y)[1]
    pch.x <- rep(pch[1], dim(x)[1])
    pch.y <- rep(pch[3], dim(y)[1])
    coltext.x <- rep(col.text[1], dim(x)[1])
    coltext.y <- rep(col.text[3], dim(y)[1])
    font.x <- rep(font[1], dim(x)[1])
    font.y <- rep(font[3], dim(y)[1])
    pr <- c("none", "active", "passive", "all")
    pri <- 1:4
    if (is.na(obj$rowsup[1])) {
        sup.x <- NA
        act.x <- x
        xn.sup <- NA
        xn.act <- x.names
    }
    else {
        sup.x <- x[obj$rowsup, , drop=FALSE]
        act.x <- x[-obj$rowsup, , drop=FALSE]
        pch.x[obj$rowsup] <- pch[2]
        coltext.x[obj$rowsup] <- col.text[2]
        font.x[obj$rowsup] <- font[2]
        xn.sup <- x.names[obj$rowsup]
        xn.act <- x.names[-obj$rowsup]
    }
    if (is.na(obj$colsup[1])) {
        sup.y <- NA
        act.y <- y
        yn.sup <- NA
        yn.act <- y.names
    }
    else {
        sup.y <- y[obj$colsup, , drop=FALSE]
        act.y <- y[-obj$colsup, , drop=FALSE]
        pch.y[obj$colsup] <- pch[4]
        coltext.y[obj$colsup] <- col.text[4]
        font.y[obj$colsup] <- font[4]
        yn.sup <- y.names[obj$colsup]
        yn.act <- y.names[-obj$colsup]
    }
    prlut <- list(none = list(x = NA, y = NA), active = list(x = act.x, 
        y = act.y), supplementary = list(x = sup.x, y = sup.y), 
        all = list(x = x, y = y))
    nameslut <- list(none = list(x.names = NA, y.names = NA), 
        active = list(x.names = xn.act, y.names = yn.act), supplementary = list(x.names = xn.sup, 
            y.names = yn.sup), all = list(x.names = x.names, 
            y.names = y.names))
    pchlut <- list(none = list(x.pch = NA, y.pch = NA), active = list(x.pch = rep(pch[1], 
        dim(x)[1]), y.pch = rep(pch[3], dim(y)[1])), supplementary = list(x.pch = rep(pch[2], 
        dim(x)[1]), y.pch = rep(pch[4], dim(y)[1])), all = list(x.pch = pch.x, 
        y.pch = pch.y))
    coltextlut <- list(none = list(x.coltext = NA, y.coltext = NA), active = list(x.coltext = rep(col.text[1], 
        dim(x)[1]), y.coltext = rep(col.text[3], dim(y)[1])), supplementary = list(x.coltext = rep(col.text[2], 
        dim(x)[1]), y.coltext = rep(col.text[4], dim(y)[1])), all = list(x.coltext = coltext.x, 
        y.coltext = coltext.y))
    fontlut <- list(none = list(x.font = NA, y.font = NA), active = list(x.font = rep(font[1], 
        dim(x)[1]), y.font = rep(font[3], dim(y)[1])), supplementary = list(x.font = rep(font[2], 
        dim(x)[1]), y.font = rep(font[4], dim(y)[1])), all = list(x.font = font.x, 
        y.font = font.y))
    x <- prlut[[pri[pr == what[1]]]][[1]]
    y <- prlut[[pri[pr == what[2]]]][[2]]
    x.names <- nameslut[[pri[pr == what[1]]]][[1]]
    y.names <- nameslut[[pri[pr == what[2]]]][[2]]
    x.pch <- pchlut[[pri[pr == what[1]]]][[1]]
    y.pch <- pchlut[[pri[pr == what[2]]]][[2]]
    x.coltext <- coltextlut[[pri[pr == what[1]]]][[1]]
    y.coltext <- coltextlut[[pri[pr == what[2]]]][[2]]
    x.font <- fontlut[[pri[pr == what[1]]]][[1]]
    y.font <- fontlut[[pri[pr == what[2]]]][[2]]
    if (is.matrix(x)) {
        x <- x[, dim, drop=FALSE]
    }
    else {
        x <- matrix(x[dim], ncol = length(dim), nrow = 1)
    }
    if (is.matrix(y)) {
        y <- y[, dim, drop=FALSE]
    }
    else {
        y <- matrix(y[dim], ncol = length(dim), nrow = 1)
    }
    if (mass[1]) 
        cex.x <- 0.5 + obj$rowmass^(1/3)/max(obj$rowmass^(1/3))
    else cex.x <- 1
    if (mass[2]) 
        cex.y <- 0.5 + obj$colmass^(1/3)/max(obj$colmass^(1/3))
    else cex.y <- 1

    # For supplementary points with mass NA
    cex.x[is.na(cex.x)] <- 1
    cex.y[is.na(cex.y)] <- 1

    nc0 <- 50
    cst <- 230
    col.x <- col[1]
    col.y <- col[2]
    if (contrib[1] == "relative") {
        cind <- obj$rowmass * (rpc[, dim[1]]^2 + rpc[, dim[2]]^2)/obj$rowinertia
        cb.x <- col2rgb(col[1])
        collut.x <- rgb(seq(cst, cb.x[1, 1], length = nc0), seq(cst, 
            cb.x[2, 1], length = nc0), seq(cst, cb.x[3, 1], length = nc0), 
            maxColorValue = 255)
        xtemp <- nc0 * (cind)
        col.x <- collut.x[xtemp]
    }
    else if (contrib[1] == "absolute") {
        cind <- obj$rowmass * (rpc[, dim[1]]^2 + rpc[, dim[2]]^2)/(obj$sv[dim[1]]^2 + 
            obj$sv[dim[2]]^2)
        cb.x <- col2rgb(col[1])
        p.x <- cb.x[, 1] + (cst - cb.x[, 1])/indx
        collut.x1 <- rgb(seq(cst, p.x[1], length = nc0/2), seq(cst, 
            p.x[2], length = nc0/2), seq(cst, p.x[3], length = nc0/2), 
            maxColorValue = 255)
        collut.x2 <- rgb(seq(p.x[1], cb.x[1, 1], length = nc0/2), 
            seq(p.x[2], cb.x[2, 1], length = nc0/2), seq(p.x[3], 
                cb.x[3, 1], length = nc0/2), maxColorValue = 255)
        collut.x <- c(collut.x1, collut.x2)
        xtemp <- nc0 * (cind)
        col.x <- collut.x[xtemp]
    }
    if (contrib[2] == "relative") {
        cind <- obj$colmass * (cpc[, dim[1]]^2 + cpc[, dim[2]]^2)/obj$colinertia
        cb.y <- col2rgb(col[2])
        collut.y <- rgb(seq(cst, cb.y[1, 1], length = nc0), seq(cst, 
            cb.y[2, 1], length = nc0), seq(cst, cb.y[3, 1], length = nc0), 
            maxColorValue = 255)
        ytemp <- nc0 * cind
        col.y <- collut.y[ytemp]
    }
    if (contrib[2] == "absolute") {
        cind <- obj$colmass * (cpc[, dim[1]]^2 + cpc[, dim[2]]^2)/(obj$sv[dim[1]]^2 + 
            obj$sv[dim[2]]^2)
        cb.y <- col2rgb(col[2])
        p.y <- cb.y[, 1] + (cst - cb.y[, 1])/indy
        collut.y1 <- rgb(seq(cst, p.y[1], length = nc0/2), seq(cst, 
            p.y[2], length = nc0/2), seq(cst, p.y[3], length = nc0/2), 
            maxColorValue = 255)
        collut.y2 <- rgb(seq(p.y[1], cb.y[1, 1], length = nc0/2), 
            seq(p.y[2], cb.y[2, 1], length = nc0/2), seq(p.y[3], 
                cb.y[3, 1], length = nc0/2), maxColorValue = 255)
        collut.y <- c(collut.y1, collut.y2)
        ytemp <- nc0 * cind
        col.y <- collut.y[ytemp]
    }
    q1 <- (1:dim(x)[1])
    q2 <- (1:dim(y)[1])
    l1 <- c(x[q1, 1], y[q2, 1])
    l1 <- l1[!is.na(l1)]
    l2 <- c(x[q1, 2], y[q2, 2])
    l2 <- l2[!is.na(l2)]
    if (length(l1) == 0) 
        l1 <- c(-0.1, 0.1)
    if (length(l2) == 0) 
        l2 <- c(-0.1, 0.1)
    lim1 <- range(l1) + c(-0.05, 0.05) * diff(range(l1))
    lim2 <- range(l2) + c(-0.05, 0.05) * diff(range(l2))
    pty.backup <- par()$pty
    plot(c(x[, 1], y[, 1]), c(x[, 2], y[, 2]), xlab = xlab, ylab = ylab, 
        type = "n", axes = FALSE, asp = 1, ...)
    box()
    abline(h = 0, v = 0, lty = 3)
    axis(1)
    axis(2)
    if (!is.na(x[1]) & labels[1] != 1) {
        if (arrows[1]) {
            arrows(rep(0, length(x[, 1])), rep(0, length(x[, 
                1])), x[, 1], x[, 2], col = col.x, length = 0.1)
        }
        else {
            points(x[, 1], x[, 2], cex = cex.x, col = col.x, 
                pch = x.pch)
        }
    }
    if (!is.na(y[1]) & labels[2] != 1) {
        if (arrows[2]) {
            arrows(rep(0, length(y[, 1])), rep(0, length(y[, 
                1])), y[, 1], y[, 2], col = col.y, length = 0.1)
        }
        else {
            points(y[, 1], y[, 2], cex = cex.y, col = col.y, 
                pch = y.pch)
        }
    }
    if (labels[1] > 0 && !is.na(x[1]) &&
        labels[2] > 0 && !is.na(y[1]))
        .pointLabel(rbind(x, y), c(x.names, y.names), cex = cex * par("cex"), xpd = TRUE,
                   col=c(x.coltext, y.coltext), font=c(x.font, y.font))
    else if (labels[1] > 0 && !is.na(x[1]))
        .pointLabel(x, x.names, cex = cex * par("cex"), xpd = TRUE, col=x.coltext, font=x.font)
    else if (labels[2] > 0 && !is.na(y[1]))
        .pointLabel(y, y.names, cex = cex * par("cex"), xpd = TRUE, col=y.coltext, font=y.font)

    par(pty = pty.backup)
}

# Function taken from the directlabels package, but it is in the public domain
.pointLabel <- function(x, y = NULL, labels = seq(along = x), cex = 1,
                       method = c("SANN", "GA"),
                       allowSmallOverlap = FALSE,
                       trace = FALSE,
                       doPlot = TRUE,
                       ...)
{
  # http://en.wikipedia.org/wiki/Automatic_label_placement
  # http://www.szoraster.com/Cartography/PracticalExperience.htm
  # http://www.eecs.harvard.edu/~shieber/Projects/Carto/carto.html
  # http://i11www.iti.uni-karlsruhe.de/map-labeling/bibliography/

  if (!missing(y) && (is.character(y) || is.expression(y))) {
    labels <- y
    y <- NULL
  }
  if (is.factor(labels)) 
    labels <- as.character(labels)
  z = xy.coords(x, y, recycle = TRUE)
  x = z$x
  y = z$y
  if (length(labels) < length(x))
    labels = rep(labels, length(x))

  method <- match.arg(method)

  boundary = par()$usr
  image_width = boundary[2] - boundary[1]
  image_height = boundary[4] - boundary[3]
  if (allowSmallOverlap) # default to 2% of the image size
    nudgeFactor = .02*(abs(boundary[1] + 1i*boundary[2] - boundary[3] - 1i*boundary[4]))

  n_labels = length(x)
                        
  # There are eight possible alignment codes, corresponding to the 
  # corners and side mid-points of the rectangle
  # Codes are 1:8
  # Code 7 is the most preferred
  xBoundary = image_width * 0.01 # add a small boundary around the rectangle
  yBoundary = image_height * 0.01
  width = strwidth(labels, units = "user", cex = cex) + xBoundary
  height = strheight(labels, units = "user", cex = cex) + yBoundary
  gen_offset <- function(code)
          c(-1, -1, -1,  0,  0,  1,  1,  1)[code] * (width/2) +
    1i * c(-1,  0,  1, -1,  1, -1,  0,  1)[code] * (height/2)


  # Finds intersection area of two rectangles
  rect_intersect <- function(xy1, offset1, xy2, offset2) {
    w = pmin(Re(xy1+offset1/2), Re(xy2+offset2/2)) - pmax(Re(xy1-offset1/2), Re(xy2-offset2/2))   
    h = pmin(Im(xy1+offset1/2), Im(xy2+offset2/2)) - pmax(Im(xy1-offset1/2), Im(xy2-offset2/2))   
    w[w <= 0] = 0
    h[h <= 0] = 0
    w*h
  }

  nudge <- function(offset) {
    # Nudge the labels slightly if they overlap:
    doesIntersect = rect_intersect(xy[rectidx1] + offset[rectidx1], rectv[rectidx1],
                                    xy[rectidx2] + offset[rectidx2], rectv[rectidx2]) > 0

    pyth = abs(xy[rectidx1] + offset[rectidx1] - xy[rectidx2] - offset[rectidx2]) / nudgeFactor
    eps = 1.0e-10

    for (i in which(doesIntersect & pyth > eps)) {
      idx1 = rectidx1[i]
      idx2 = rectidx2[i]
      vect = (xy[idx1] + offset[idx1] - xy[idx2] - offset[idx2]) / pyth[idx1]
      offset[idx1] = offset[idx1] + vect
      offset[idx2] = offset[idx2] - vect
    }
    offset
  }

  objective <- function(gene) {
    offset = gen_offset(gene)

    # Allow for "bending" the labels a bit
    if (allowSmallOverlap) offset = nudge(offset)

    if (!is.null(rectidx1))
      area = sum(rect_intersect(xy[rectidx1] + offset[rectidx1], rectv[rectidx1],
                                xy[rectidx2] + offset[rectidx2], rectv[rectidx2]))
    else
      area = 0
    
    # Penalize labels which go outside the image area
    # Count points outside of the image
    n_outside = sum(Re(xy + offset - rectv/2) < boundary[1] | Re(xy + offset + rectv/2) > boundary[2] |
                    Im(xy + offset - rectv/2) < boundary[3] | Im(xy + offset + rectv/2) > boundary[4]) 
    area + n_outside * image_width * image_height
  }

  
  # Make a list of label rectangles in their reference positions,
  # centered over the map feature; the real labels are displaced
  # from these positions so as not to overlap
  # Note that some labels can be bigger than others
  xy = x + 1i * y
  rectv = width + 1i * height

  rectidx1 = rectidx2 = array(0, (length(x)^2 - length(x)) / 2)
  k=0
  for (i in 1:length(x))
    for (j in seq(len=(i-1))) {
      k = k + 1
      rectidx1[k] = i
      rectidx2[k] = j
    }
  canIntersect = rect_intersect(xy[rectidx1], 2 * rectv[rectidx1],
                                xy[rectidx2], 2 * rectv[rectidx2]) > 0
  rectidx1 = rectidx1[canIntersect]
  rectidx2 = rectidx2[canIntersect]
  if (trace) cat("possible intersects =", length(rectidx1), "\n")

  if (trace) cat("portion covered =", sum(rect_intersect(xy, rectv,xy,rectv))/(image_width*image_height),"\n")

  SANN <- function() {
    # Make some starting "genes"
    #gene = sample(1:8, n_labels, repl = TRUE)
    gene = rep(8, n_labels)
    score = objective(gene)
    bestgene = gene
    bestscore = score
    T = 2.5
    for (i in 1:50) {
      k = 1
      for (j in 1:50) {
        newgene = gene
        newgene[sample(1:n_labels, 1)] = sample(1:8,1)
        newscore = objective(newgene)
        if (newscore < score || runif(1) < 1 - exp((newscore - score) / T)) {
          k = k + 1
          score = newscore
          gene = newgene
        }
        if (score <= bestscore) {
          bestscore = score
          bestgene = gene
        }
        if (bestscore == 0 || k == 10) break
      }
      if (bestscore == 0) break
      if (trace) cat("overlap area =", bestscore, "\n")
      T = 0.9 * T
    }
  
    if (trace) cat("overlap area =", bestscore, "\n")
    nx = Re(xy + gen_offset(bestgene))
    ny = Im(xy + gen_offset(bestgene))
    list(x = nx, y = ny)
  }

  xy = SANN()

  # Taken from http://article.gmane.org/gmane.comp.lang.r.general/147787
  shadowtext <- function(xy, labels, col='black', bg='white',
                          theta=seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
      
          xy <- xy.coords(xy)
          xo <- r*strwidth('A')
          yo <- r*strheight('A')

          for (i in theta)
              text(xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ...)

          text(xy$x, xy$y, labels, col=col, ... )
  }

  if (doPlot)
    shadowtext(xy, labels, cex = cex, ...)

  invisible(xy)
}
