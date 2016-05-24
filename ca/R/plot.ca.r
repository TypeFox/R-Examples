################################################################################
# plot.ca(): Plotting ca objects (ca package 0.70)
################################################################################
plot.ca <- function(x, 
                    dim     = c(1,2), 
                    map     = "symmetric", 
                    what    = c("all", "all"), 
                    mass    = c(FALSE, FALSE), 
                    contrib = c("none", "none"), 
                    col     = c("blue", "red"), 
                    pch     = c(16, 21, 17, 24), 
                    labels  = c(2,2), 
                    arrows  = c(FALSE, FALSE), 
                    lines   = c(FALSE, FALSE),
                    lwd     = 1,
                    xlab    = "_auto_",
                    ylab    = "_auto_",
                    col.lab = c("blue", "red"),
                    ...){
  obj <- x
 # recycling input:
  if (length(what) != 2){
    what <- rep(what, length = 2)
    }
  if (length(mass) != 2){
    mass <- rep(mass, length = 2)
    }
  if (length(contrib) != 2){
    contrib <- rep(contrib, length = 2)
    }
  if (length(col) != 2){
    col <- rep(col, length = 2)
    }
  if (length(labels) != 2){
    labels <- rep(labels, length = 2)
    }
  if (length(pch) != 4){
    pch <- rep(pch, length = 4)
    }
  if (length(lines) != 2){
    lines <- rep(lines, length = 2)
    }
  if (length(col.lab) != 2){
    col.lab <- rep(col.lab, length = 2)
    }
 # check for suprow/-col and 'row-/col.gab/-green'
  if (!is.numeric(x$suprow)) {
    if (map == "colgab" | map == "colgreen") {
      if (what[1] != "none") what[1] <- "active"
      }
    }
  if (!is.numeric(x$supcol)) {
    if (map == "rowgab" | map == "rowgreen") {
      if (what[2] != "none") what[2] <- "active"
      }
    }
 # sign switching:
  if (min(dim) < 0){
    swisign      <- ifelse(dim < 0, -1, 1)
    dim.c        <- dim(obj$rowcoord)[2]
    signmat      <- diag(rep(swisign, length = dim.c))
    obj$rowcoord <- obj$rowcoord%*%signmat
    obj$colcoord <- obj$colcoord%*%signmat
    dim          <- abs(dim)
    }
 # principal coordinates:
  K      <- dim(obj$rowcoord)[2]
  I      <- dim(obj$rowcoord)[1] ; J <- dim(obj$colcoord)[1]
  svF    <- matrix(rep(obj$sv[1:K], I), I, K, byrow = TRUE)
  svG    <- matrix(rep(obj$sv[1:K], J), J, K, byrow = TRUE)
  rpc    <- obj$rowcoord * svF
  cpc    <- obj$colcoord * svG
  symrpc <- obj$rowcoord * sqrt(svF)
  symcpc <- obj$colcoord * sqrt(svG)
 # maptype
  mt <- c("symmetric", "rowprincipal", "colprincipal", "symbiplot", "rowgab", "colgab", "rowgreen", "colgreen")
  mti <- 1:length(mt)
  mtlut <- list(symmetric    = list(x = rpc, y = cpc),
                rowprincipal = list(x = rpc, y = obj$colcoord),
                colprincipal = list(x = obj$rowcoord, y = cpc),
                symbiplot    = list(x = symrpc, y = symcpc), 
                rowgab       = list(x = rpc, y = obj$colcoord * obj$colmass),
                colgab       = list(x = obj$rowcoord * obj$rowmass, y = cpc), 
                rowgreen     = list(x = rpc, y = obj$colcoord * sqrt(obj$colmass)), 
                rowgreen     = list(x = obj$rowcoord * sqrt(obj$rowmass), y = cpc)
                )
  x       <- mtlut[[mti[mt==map]]][[1]]
  y       <- mtlut[[mti[mt==map]]][[2]]
  x.names <- obj$rownames
  y.names <- obj$colnames
 # profiles to plot
  indx  <- dim(x)[1]
  indy  <- dim(y)[1]
  pch.x <- rep(pch[1],dim(x)[1])
  pch.y <- rep(pch[3],dim(y)[1])
  pr    <- c("none", "active", "passive", "all")
  pri   <- 1:4
  if (is.na(obj$rowsup[1])) {
    sup.x  <- NA
    act.x  <- x
    xn.sup <- NA
    xn.act <- x.names
    } else {
      sup.x  <- x[obj$rowsup,]
      act.x  <- x[-obj$rowsup,]
      pch.x[obj$rowsup] <- pch[2]
      xn.sup <- x.names[obj$rowsup]
      xn.act <- x.names[-obj$rowsup]
      }
  if (is.na(obj$colsup[1])) {
    sup.y  <- NA
    act.y  <- y
    yn.sup <- NA
    yn.act <- y.names
    } else {
      sup.y  <- y[obj$colsup,]
      act.y  <- y[-obj$colsup,]
      pch.y[obj$colsup] <- pch[4]
      yn.sup <- y.names[obj$colsup]
      yn.act <- y.names[-obj$colsup]
      }
  prlut <- list(none          = list(x = NA, y = NA), 
                active        = list(x = act.x, y = act.y),
                supplementary = list(x = sup.x, y = sup.y),
                all           = list(x = x, y = y))
  nameslut <- list(none          = list(x.names = NA, y.names = NA),
                   active        = list(x.names = xn.act, y.names = yn.act),
                   supplementary = list (x.names = xn.sup, y.names = yn.sup),
                   all           = list(x.names = x.names, y.names = y.names) )
  pchlut <- list(none          = list(x.pch = NA, y.pch = NA),
                 active        = list(x.pch = rep(pch[1],dim(x)[1]), y.pch = rep(pch[3],dim(y)[1])),
                 supplementary = list(x.pch = rep(pch[2],dim(x)[1]), y.pch = rep(pch[4],dim(y)[1])),
                 all           = list(x.pch = pch.x, y.pch = pch.y) )
  x       <- prlut[[pri[pr==what[1]]]][[1]]
  y       <- prlut[[pri[pr==what[2]]]][[2]]
  x.names <- nameslut[[pri[pr==what[1]]]][[1]]
  y.names <- nameslut[[pri[pr==what[2]]]][[2]]
  x.pch   <- pchlut[[pri[pr==what[1]]]][[1]]
  y.pch   <- pchlut[[pri[pr==what[2]]]][[2]]
  
 # dimensions to plot
  if(is.matrix(x)){
    x <- x[,dim]
    } else {
    x <- matrix(x[dim],ncol = length(dim), nrow = 1)
    }
  if(is.matrix(y)){
    y <- y[,dim]
    } else {
    y <- matrix(y[dim], ncol = length(dim), nrow = 1) }

 ## plot setup
 # radius/mass
  if (mass[1]){
    cex.x <- 0.5 + obj$rowmass^(1/3) / max(obj$rowmass^(1/3))
    } else {
    cex.x <- 1
    }
  if (mass[2]){
    cex.y <- 0.5 + obj$colmass^(1/3) / max(obj$colmass^(1/3))
    } else {
    cex.y <- 1
    }
 # contributions/colour intensities
  nc0 <- 50
  cst <- 230
  col.x <- col[1]
  col.y <- col[2]
  if (contrib[1] == "relative") {
    cind     <- obj$rowmass * (rpc[,dim[1]]^2 + rpc[,dim[2]]^2) / obj$rowinertia
    cb.x     <- col2rgb(col[1])
    collut.x <- rgb(seq(cst, cb.x[1, 1], length = nc0),
                    seq(cst, cb.x[2, 1], length = nc0),
                    seq(cst, cb.x[3, 1], length = nc0), maxColorValue = 255 )
    xtemp <- nc0*(cind)
    col.x <- collut.x[xtemp]
    }  else {
    if (contrib[1] == "absolute") {
      cind <- obj$rowmass*(rpc[,dim[1]]^2 + rpc[,dim[2]]^2) / (obj$sv[dim[1]]^2 + 
                obj$sv[dim[2]]^2)
      cb.x <- col2rgb(col[1])
      p.x <- cb.x[,1] + (cst - cb.x[,1])/indx
      collut.x1 <- rgb(seq(cst, p.x[1], length = nc0/2),
                       seq(cst, p.x[2], length = nc0/2),
                       seq(cst, p.x[3], length = nc0/2), maxColorValue = 255)
      collut.x2 <- rgb(seq(p.x[1], cb.x[1, 1], length = nc0/2),
                       seq(p.x[2], cb.x[2, 1], length = nc0/2),
                       seq(p.x[3], cb.x[3, 1], length = nc0/2), maxColorValue = 255)
      collut.x <- c(collut.x1, collut.x2)
      xtemp <- nc0*(cind)
      col.x <- collut.x[xtemp]
      }
    }
  if (contrib[2] == "relative") {
    cind <- obj$colmass*(cpc[,dim[1]]^2 + cpc[,dim[2]]^2) / obj$colinertia
    cb.y <- col2rgb(col[2])
    collut.y <- rgb(seq(cst, cb.y[1, 1], length = nc0),
                    seq(cst, cb.y[2, 1], length = nc0),
                    seq(cst, cb.y[3, 1], length = nc0), maxColorValue = 255 )
    ytemp <- nc0 * cind
    col.y <- collut.y[ytemp]
    } 
  if (contrib[2] == "absolute") {
    cind <- obj$colmass*(cpc[,dim[1]]^2 + cpc[,dim[2]]^2) / (obj$sv[dim[1]]^2 + 
              obj$sv[dim[2]]^2)
    cb.y <- col2rgb(col[2])
    p.y <- cb.y[,1] + (cst - cb.y[,1])/indy
    collut.y1 <- rgb(seq(cst, p.y[1], length = nc0/2),
                     seq(cst, p.y[2], length = nc0/2),
                     seq(cst, p.y[3], length = nc0/2), maxColorValue = 255 )
    collut.y2 <- rgb(seq(p.y[1], cb.y[1, 1], length = nc0/2),
                     seq(p.y[2], cb.y[2, 1], length = nc0/2),
                     seq(p.y[3], cb.y[3, 1], length = nc0/2), maxColorValue = 255 )
    collut.y <- c(collut.y1, collut.y2)
    ytemp <- nc0 * cind
    col.y <- collut.y[ytemp]
    }

## plotting:
 # determine margins
  q1 <- (1:dim(x)[1])
  q2 <- (1:dim(y)[1])
  l1 <- c(x[q1,1], y[q2,1]) ; l1 <- l1[!is.na(l1)]
  l2 <- c(x[q1,2], y[q2,2]) ; l2 <- l2[!is.na(l2)]
  if (length(l1) == 0) l1 <- c(-.1, .1)
  if (length(l2) == 0) l2 <- c(-.1, .1)
  lim1 <- range(l1) + c(-.05, .05) * diff(range(l1))
  lim2 <- range(l2) + c(-.05, .05) * diff(range(l2))
 # axis labels
  pct <- round(100* (obj$sv^2) / sum(obj$sv^2), 1)
  pct <- paste0(" (", pct[dim], "%)")
  if (xlab == "_auto_"){
    xlab = paste0("Dimension ", dim[1], pct[1])
    }  
  if (ylab == "_auto_"){
    ylab = paste0("Dimension ", dim[2], pct[2])
    }
  pty.backup <- par()$pty
 # plot:
  plot(c(x[,1],y[,1]), c(x[,2],y[,2]), xlab = xlab, ylab = ylab, type = "n", axes = FALSE, asp = 1, ...)
  box()
  abline(h = 0, v = 0, lty = 3)
  axis(1)
  axis(2)
 # rows
  if (!is.na(x[1]) & labels[1] != 1) {
    if (arrows[1]) {
      .arrows(rep(0, length(x[,1])), rep(0, length(x[,1])), x[,1], x[,2], col = col.x, lwd=lwd, length = 0.1) 
      } else {
      points(x[,1], x[,2], cex = cex.x, col = col.x, pch = x.pch)
      }
    }
  if (labels[1] > 0) {
    xoff1 <- if(labels[1]>1) .5 * strwidth(x.names, cex = .75) + .5 * strwidth("o", cex = .75) else 0
    xoff2 <- if(labels[1]>1) .5 * strheight(x.names, cex = .75) + .5 * strheight("o", cex = .75) else 0
    text(x[,1] + xoff1, x[,2] + xoff2, x.names, cex = 0.75, xpd = TRUE, col = col.lab[1])
    }
 # columns
  if (!is.na(y[1]) & labels[2] != 1 ) {
    if (arrows[2]) {
      .arrows(rep(0, length(y[,1])), rep(0, length(y[,1])), y[,1], y[,2], col = col.y, lwd=lwd, length = 0.1) 
      } else {
      points(y[,1], y[,2], cex = cex.y, col = col.y, pch = y.pch)
      }
    }
  if (labels[2] > 0) {
    yoff1 <- if(labels[2]>1) .5 * strwidth(y.names, cex = 0.75) + .5 * strwidth("o", cex = .75) else 0
    yoff2 <- if(labels[2]>1) .5 * strheight(y.names, cex = 0.75) + .5 * strheight("o", cex = .75) else 0
    text(y[,1] + yoff1, y[,2] + yoff2, y.names, cex = 0.75, xpd = TRUE, col = col.lab[2])
    }
 # plot connecting lines (sorted by X value)
  if (lines[1]) lines(x[order(x[,1]),], col = col.x, lwd=lwd)
  if (lines[2]) lines(y[order(y[,1]),], col = col.y, lwd=lwd)
  par(pty = pty.backup)

 # return a result for further plot annotation
  rownames(x) <- x.names; colnames(x) <- paste0("Dim", dim)
  rownames(y) <- y.names; colnames(y) <- paste0("Dim", dim)
  result <- list(rows = x, cols = y)
  invisible(result)
  }
################################################################################

# the following function isn't exported
# Provides a simple way to make more attractive arrows

.arrows <- function(..., angle=15){
  angles <- seq(1, angle, by=2)
  for (ang in angles) arrows(..., angle=ang)
}

