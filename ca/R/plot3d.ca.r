################################################################################
# plot3d.ca(): 3D-Plotting of ca objects (ca package 0.70)
################################################################################
plot3d.ca <- function(x, 
                      dim       = c(1,2,3), 
                      map       = "symmetric", 
                      what      = c("all","all"), 
                      contrib   = c("none","none"), 
                      col       = c("#6666FF","#FF6666"), 
                      labcol    = c("#0000FF","#FF0000"), 
                      pch       = c(16,1,18,9), 
                      labels    = c(2,2), 
                      sf        = 0.00001,
                      arrows    = c(FALSE,FALSE),
                      axiscol   = "#333333",
                      axislcol  = "#333333",
                      laboffset = list(x = 0, y = 0.075, z = 0.05),
                      ...){
  #require(rgl)
  requireNamespace("rgl")
 ########## RGLPLOT0: temporary solution for using 'pch' within RGL ##########
  rglplot0 <- function(x = 0, y = 0, z = 0, v = 1, pch = 1, segments = 16, size = 1, ...) {
    # coordinates for circles/tetraeders etc
    i  <- seq(0, 360, length = segments)
    xc <- rep(sin(pi*i/180), each = 2)
    yc <- rep(cos(pi*i/180), each = 2)
    xc <- xc[c(-1,-length(xc))]
    yc <- yc[c(-1,-length(yc))]
    i0 <- c(0,109.5,109.5,0,109.5,109.5,0,109.5,109.5,109.5,109.5,109.5)
    j0 <- c(0,0,-120,0,-120,120,0,120,0,-120,0,120)+90
    i1 <- c(0,109.5,0,109.5,0,109.5,109.5,109.5,109.5,109.5,109.5,109.5)
    j1 <- c(0,0,0,-120,0,120,120,0,0,-120,-120,120)+90
    i2 <- c(rep(c(180,81.5), 4), rep(81.5, 8))
    j2 <- c(0,-45,0,45,0,135,0,-135,-135,-45,-45,45,45,135,135,-135)
    i3 <- c(i1, 180-i1)
    j3 <- c(j1, j1-60)
    r  <- v
   # 'base-matrices' for primitives
    c01 <- list(x = c(xc, yc, rep(0, length(xc))), 
                z = c(yc, rep(0, length(xc)), xc),
                y = c(rep(0, length(xc)), xc, yc))
    c02 <- list(x = r*sin(pi*i1/180)*cos(pi*j1/180), 
                z = r*cos(pi*i1/180), 
                y = r*sin(pi*i1/180)*sin(pi*j1/180))
    c03 <- list(x = c(-1,1,0,0,0,0), z = c(0,0,-1,1,0,0), y = c(0,0,0,0,-1,1))
    c04 <- list(x = c(-1,1,1,-1,1,-1,-1,1), 
                z = c(1,-1,-1,1,1,-1,-1,1), 
                y = c(1,-1,1,-1,-1,1,-1,1) )
    c05 <- list(x = c(-1,0,0,0,1,0,0,0,-1,0,-1,0,1,0,1,0,0,0,-1,0,0,0,1,0),
                z = c(0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,-1,0,0,-1,-1,0,0,-1),
                y = c(0,0,0,1,0,0,0,-1,0,-1,0,1,0,1,0,-1,0,-1,0,0,0,1,0,0) )
    c06 <- list(x = r*sin(pi*i2/180)*cos(pi*j2/180), 
                z = r*cos(pi*i2/180), 
                y = r*sin(pi*i2/180)*sin(pi*j2/180))
    c11 <- list(x = r*sin(pi*i3/180)*cos(pi*j3/180), 
                z = r*cos(pi*i3/180), 
                y = r*sin(pi*i3/180)*sin(pi*j3/180))
    c15 <- list(x = c(-1, -1, rep(1, 8), rep(-1,8), 1, 1, -1, 1, 1, -1), 
                z = c(1, 1, 1, 1, rep(c(-1,1,1,-1,1,-1,-1,1), 2), rep(-1, 4)), 
                y = c(-1,1,1,-1,-1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1))
    c16 <- list(x = 0, z = 0, y = 0)
    c17 <- list(x = r*sin(pi*i0/180)*cos(pi*j0/180), 
                z = r*cos(pi*i0/180), 
                y = r*sin(pi*i0/180)*sin(pi*j0/180))
    c18 <- list(x = 0.8*c(-1,0,1,1,0,1,1,0,-1,-1,0,-1,1,0,-1,1,0,1,-1,0,1,-1,0,-1),
                z = c(rep(c(0,1,0), 4), rep(c(0,-1,0), 4)),
                y = 0.8*c(-1,0,-1,-1,0,1,1,0,1,1,0,-1,-1,0,-1,1,0,-1,1,0,1,-1,0,1))
    c19 <- list(x = 0, z = 0, y = 0)
    c20 <- list(x = 0, z = 0, y = 0)
    c22 <- list(x = c(-1,1,1,1,-1,1,rep(c(-1,1), each=4),rep(-1,5),1,1,1,-1,1),
                z = c(rep(1,9),-1,1,-1,1,-1,1,rep(-1,9)),
                y = c(-1,-1,-1,1,1,1,-1,1,rep(c(-1,1),each=4),-1,1,-1,-1,-1,1,1,1))
    c07 <- list(x = c(c22$x, c04$x), z = c(c22$z, c04$z), y = c(c22$y, c04$y))
    c08 <- list(x = c(c03$x, c04$x), z = c(c03$z, c04$z), y = c(c03$y, c04$y))
    c09 <- list(x = c(c05$x, c03$x), z = c(c05$z, c03$z), y = c(c05$y, c03$y))
    c10 <- list(x = c(c01$x, c03$x), z = c(c01$z, c03$z), y = c(c01$y, c03$y))
    c12 <- list(x = c(c22$x, c03$x), z = c(c22$z, c03$z), y = c(c22$y, c03$y))
    c13 <- list(x = c(c01$x, c04$x), z = c(c01$z, c04$z), y = c(c01$y, c04$y))
    c14 <- list(x = c(c22$x,1,0,1,0,-1,0,-1,0),
                z = c(c22$z,1,-1,1,-1,1,-1,1,-1),
                y = c(c22$y,-1,0,1,0,1,0,-1,0))
    c21 <- c01; c23 <- c05; c24 <- c02; c25 <- c06
   # rgl-function/scaling "lookup-table" 
    cc <- list(c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, c11, c12, c13,
               c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25)
    cc.scale <- c( 3*v/(4*pi), 4*v/sqrt(32), 3*v/8, v/8, 3*v/8,
                   3*v/4, v/8, v/8, 3*v/8, 3*v/(4*pi),
                   9*v/(4*sqrt(50)), v/8, 3*v/(4*pi), v/8, v/8,
                   3*v/(4*pi), 4*v/sqrt(32), 3*v/8, 3*v/(4*pi), 3*v/(4*pi),
                   3*v/(4*pi), v/8, v*3/8, 4*v/sqrt(32), 3*v/4 )^(1/3)
    pchlevels <- as.numeric(levels(as.factor(pch)))
    for (k in 1:length(pchlevels)) {
      coord <- cc[[pchlevels[k]]]
      coord.scale <- cc.scale[pchlevels[k]]
      fm <- c(rep("rgl::lines3d", 14), "rgl::quads3d", "rgl::spheres3d", 
              rep("rgl::triangles3d", 2), rep("rgl::spheres3d", 2), 
              rep("rgl::lines3d", 5))
      fm.suffix <- c(rep("",15),paste(coord.scale,",",sep=""),"","",
                     rep(paste(coord.scale,",",sep=""),2),rep("",5))
      x0 <- rep(x[pch==pchlevels[k]], each = length(coord$x)) + coord$x * coord.scale
      z0 <- rep(z[pch==pchlevels[k]], each = length(coord$z)) + coord$z * coord.scale
      y0 <- rep(y[pch==pchlevels[k]], each = length(coord$y)) + coord$y * coord.scale
      eval(parse(text = paste(fm[pchlevels[k]], "(x0, y0, z0,", 
        fm.suffix[pchlevels[k]], "size = size, ...)", sep = "")))
      }
    }
  rglarrows0 <- function(x, y, z, col = "white"){
    x.new <- as.vector(rbind(rep(0,length(x)),x))
    y.new <- as.vector(rbind(rep(0,length(y)),y))
    z.new <- as.vector(rbind(rep(0,length(z)),z))
    rgl::rgl.lines(x.new, y.new, z.new, color = col)
    }
 ########## END OF RGLPLOT0

 # recycling input:
  if (length(what) != 2){
    what <- rep(what, length = 2)
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
  obj <- x
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
  mt    <- c("symmetric", "rowprincipal", "colprincipal", "symbiplot", "rowgab", 
             "colgab", "rowgreen", "colgreen")
  mti   <- 1:length(mt)
  mtlut <- list(symmetric    = list(x = rpc, y = cpc),
                rowprincipal = list(x = rpc, y = obj$colcoord),
                colprincipal = list(x = obj$rowcoord, y = cpc),
                symbiplot    = list(x = symrpc, y = symcpc), 
                rowgab       = list(x = rpc, y = obj$colcoord * obj$colmass),
                colgab       = list(x = obj$rowcoord * obj$rowmass, y = cpc), 
                rowgreen     = list(x = rpc, y = obj$colcoord * sqrt(obj$colmass)), 
                colgreen     = list(x = obj$rowcoord * sqrt(obj$rowmass), y = cpc)
                )
  x       <- mtlut[[mti[mt==map]]][[1]]
  y       <- mtlut[[mti[mt==map]]][[2]]
  x.names <- obj$rownames
  y.names <- obj$colnames
  rm(mt, mti, mtlut)

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
                   supplementary = list(x.names = xn.sup, y.names = yn.sup),
                   all           = list(x.names = x.names, y.names = y.names) )
  pchlut <- list(none          = list(x.pch = NA, y.pch = NA),
                 active        = list(x.pch = rep(pch[1],dim(x)[1]), y.pch = rep(pch[3],dim(y)[1])),
                 supplementary = list (x.pch = rep(pch[2],dim(x)[1]), y.pch = rep(pch[4],dim(y)[1])),
                 all           = list(x.pch = pch.x, y.pch = pch.y) )

  x       <- prlut[[pri[pr == what[1]]]][[1]]
  y       <- prlut[[pri[pr == what[2]]]][[2]]
  x.names <- nameslut[[pri[pr == what[1]]]][[1]]
  y.names <- nameslut[[pri[pr == what[2]]]][[2]]
  x.pch   <- pchlut[[pri[pr == what[1]]]][[1]]
  y.pch   <- pchlut[[pri[pr == what[2]]]][[2]]

 # dimensions to plot
  if(is.matrix(x)){
    x <- x[,dim]
    } else {
    x <- matrix(x[dim], ncol = length(dim), nrow = 1)
    }
  if(is.matrix(y)){
    y <- y[,dim]
    } else {
    y <- matrix(y[dim], ncol = length(dim), nrow = 1)
    }
## plot setup
 # radius/mass
  cex.x <- cex.y <- sf
 # contributions/colour intensities
  calpha.x <- 1
  calpha.y <- 1
  if (contrib[1] == "relative") {
    calpha.x <- obj$rowmass*(rpc[,dim[1]]^2 + rpc[,dim[2]]^2) / obj$rowinertia
    } else { 
    if (contrib[1] == "absolute") {
      calpha.x <- obj$rowmass*(rpc[,dim[1]]^2 + rpc[,dim[2]]^2) / (obj$sv[dim[1]]^2 + obj$sv[dim[2]]^2)
      }
    }
  if (contrib[2] == "relative") {
    calpha.y <- obj$colmass*(cpc[,dim[1]]^2 + cpc[,dim[2]]^2) / obj$colinertia
    } else {
    if (contrib[2] == "absolute") {
      calpha.y <- obj$colmass*(cpc[,dim[1]]^2 + cpc[,dim[2]]^2) / (obj$sv[dim[1]]^2 + obj$sv[dim[2]]^2)
      }
    }
 ## plotting:
  rgl::plot3d(0, 0, 0, xlab = "", ylab = "", zlab = "", type = "n", box = FALSE, axes = FALSE, aspect = TRUE)
  xt.1 <- c(x[,1], y[,1])
  xt.2 <- c(x[,2], y[,2])
  xt.3 <- c(x[,3], y[,3])
  x001 <- range(xt.1[!is.na(xt.1)])
  x002 <- range(xt.2[!is.na(xt.2)])
  x003 <- range(xt.3[!is.na(xt.3)])
  laboffset$x <- laboffset$x * abs(diff(range(x001)))
  laboffset$y <- laboffset$y * abs(diff(range(x002)))
  laboffset$z <- laboffset$z * abs(diff(range(x003)))
  rgl::lines3d(x001, c(0, 0), c(0, 0), col = axiscol)
  rgl::lines3d(c(0, 0), x002, c(0, 0), col = axiscol)
  rgl::lines3d(c(0, 0), c(0, 0), x003, col = axiscol)
 # axis labels
  rgl::texts3d(x001[2]+0.01*abs(diff(range(x001))), 0, 0, dim[1], col = axislcol, cex = 0.75)
  rgl::texts3d(0, x002[2]+0.01*abs(diff(range(x001))), 0, dim[2], col = axislcol, cex = 0.75)
  rgl::texts3d(0, 0, x003[2]+0.01*abs(diff(range(x001))), dim[3], col = axislcol, cex = 0.75)
 # rows
  if (!is.na(x[1])) {
    if (labels[1] != 1) {
      if (arrows[1]) {
        rglarrows0(x[,1], x[,2], x[,3], col = col[1]) 
        } else {
        rglplot0(x[,1], x[,2], x[,3], col = col[1], v = cex.x, alpha = calpha.x, pch = x.pch)
        }
      }
    if (labels[1] > 0) {
      rgl::texts3d(x[,1]+laboffset$x, x[,2]+laboffset$y, x[,3]+laboffset$z, x.names, cex = 0.95, col = labcol[1])
      }
    }
 # columns
  if (!is.na(y[1])) {
    if (labels[2] != 1 ) {
      if (arrows[2]) {
        rglarrows0(y[,1], y[,2], y[,3], col = col[2]) 
        } else {
        rglplot0(y[,1], y[,2], y[,3], col = col[2], v = cex.y, alpha = calpha.y, pch = y.pch)
        }
      }
    if (labels[2] > 0) {
      rgl::texts3d(y[,1]+laboffset$x, y[,2]+laboffset$y, y[,3]+laboffset$z, y.names, cex = 0.95, col = labcol[2])
      }
    }
  }
################################################################################
