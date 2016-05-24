simul.simplefield <-
function(gdi="gqi", b=3000, sigma=NULL, clusterthr=0.6, logplot=TRUE,
 savedir=tempdir(), fmask="m1", ang=NULL, ...)
{
  ## glyph field mask 
  tfmask=c("m1","m2","m3","mx1","mx2","mx3")
  ## default angles in masks
  amask =c(60,60,60,60,90,45)
  mfun <- match(fmask, tfmask)
  if(is.null(ang))
    ang <- amask[mfun]
  switch(mfun,
    a <- m1(ang),
    a <- m2(ang),
    a <- m3(ang),
    a <- mx1(ang),
    a <- mx2(ang),
    a <- mx3(ang))
  ## S2 shell grid
  s2 <- s2tessel.zorder(depth=3, viewgrid=FALSE)
  g0 <- s2$pc
  ## swap mask (from top left in row order to column order)
  if(is.matrix(a)){
    nr <- ncol(a); nc <- nrow(a)
    as <- matrix(-1,nr,nc) 
    as[nr:1,1:nc] <- t(a)
   }
  else 
    if(is.array(a)) {
      dm <- dim(a)
      as <- array(0, dim=dm)
      as[dm[1]:1, 1:dm[2], ] <- abind::abind(t(a[,,1]), t(a[,,2]), along=3)
    } else
      stop("Error in mask specification")
  ## field simulation 
	gc()
  field <- myglyph.synthsimul(as, ang=ang, g0, b=b, sigma=sigma,
                              logplot=logplot)
  ## Estimate ODFs
  odfs <- fieldtestodf.gqi(gdi=gdi, g0, field,  b=b, lambda=NULL,
    savedir=savedir)
  ## Visualize grid of glyphs with color
  plotodfvxgrid(g0, field=field, odfsdata=odfs)
  ## Using movMF for mixture estimation and peak detection)
  estimate.vmf.lines(g0, field=field, odfsdata=odfs,
    showglyph=FALSE, clusterthr=clusterthr, savedir=savedir, ...)
}

myglyph.synthsimul <-
function(as, ang=60, g0, b=3000, sigma=NULL, logplot=TRUE)
{
  if(is.matrix(as)){
    nr <- nrow(as); nc <- ncol(as)
    nn <- length(which(as != -1))
  } else 
    if(is.array(as)) {
      dd <- dim(as)
      nr <- dd[2]; nc <- dd[1]
      nn <- length(which(as[,,1] != -1))
    } else
      stop("Error in field mask specification")
  ## synthesis
  ng0 <- dim(g0)[1]
  S <- matrix(0, nn, ng0)
  k <- 0
  open3d()
  for(j in 1:nc) {
    for(i in 1:nr) {
      # single fiber
      if(is.matrix(as)) {
        ang <- as[i,j]
        if(ang == -1) next
          pos <- 2 * c(i,j,0)
          sv <- synthfiberss2z(g0=g0, angles=ang, b=b, sigma=sigma,
            pos=pos, showglyph=TRUE, new=FALSE, logplot=logplot)
          k <- k+1
          S[k,] <- sv
        mask <- as
      }
      else {
        if(is.array(as)) {
          ang1 <- as[i,j,1]
          if(ang1 == -1) next
          ang2 <- as[i,j,2]
          if(ang2 == -1) {
             pos <- 2 * c(i,j,0)
               sv <- synthfiberss2z(g0=g0, angles=ang1, b=b, sigma=sigma,
                      pos=pos, showglyph=TRUE, new=FALSE, logplot=logplot)
              k <- k+1
              S[k,] <- sv
          } else {
             pos <- 2 * c(i,j,0)
            sv <- synthfiberss2z(g0=g0, angles=c(ang1,ang2), b=b, sigma=sigma,
                     pos=pos, showglyph=TRUE, new=FALSE, logplot=logplot)
            k <- k+1
             S[k,] <- sv
            mask <- as[,,1]
          }
        }
      }
    }
  }
  list(S=S, mask=mask)
}

##--------------------------
fieldtestodf.gqi <-
function(gdi="gqi", grad, field, b=3000, lambda=NULL,
 savedir=tempdir())
{
  cat("estimating field odfs ...\n")
  gdimethods <- c("gqi", "gqi2")
  gdimethod <- match(gdi, gdimethods)
  sfield <- field$S
  dv <- dim(sfield)
  mask <- field$mask
  dm <- dim(mask)
  nr <- dm[1]; nc <- dm[2]
  odfs <- matrix(0, dv[1], dv[2])
  bn <- rep(b, dim(grad)[1])
  btable <- as.matrix(cbind(bn, grad))
  ##-----------------------------
  ## "gdimethod" process
  cat("Estimating slice odfs ...\n")
  switch(gdimethod,
      q2odf <- gqifn(odfvert=grad, btable=btable,
                     lambda=lambda),
      q2odf <- gqifn2(odfvert=grad, btable=btable,
                     lambda=lambda) )
  ##-----------------------------
  k <- 0
  for(j in 1:nc) {
    for(i in 1:nr) {
      if(mask[i,j] == -1) next
      k <- k+1
      odf <- as.vector(q2odf%*%sfield[k,])
      odf <- odf - min(odf)
      odfs[k,] <- odf
    }
  }
  f <- file.path(savedir,"simtable.txt")
  write(t(btable), file=f, ncolumns=4)
  cat("wrote",f,"\n")
  invisible(odfs)
}

#--------------------------

plotodfvxgrid <-
function(pc0, field, odfsdata)
{
  mask <- field$mask
  dm <- dim(mask)
  nr <- dm[1]; nc <- dm[2]
  ## GFA
  odfs.reg <- t(apply(odfsdata, 1 , norm01))
  gfas <- apply(odfs.reg, 1, genfa)
  tc <-  geometry::delaunayn(pc0)
  tc.surf <- t( surf.tri(pc0,tc) )
  dt2 <- dim(pc0[tc.surf,])[1]
  d1 <- dim(odfsdata)[1]
  sgrid <- matrix(0, nrow=d1*dt2, ncol=3)
  vcolors <- matrix(0, nrow=d1, ncol=dt2)
  k <- 0
  cat("Running ...\n")
  for(j in 1:nc) {
    for(i in 1:nr) {
      if(mask[i,j] == -1) next
      k <- k+1
      odf <- odfs.reg[k,]
      gk <- gfas[k]
      ## RGB channels
      zch <- pc0*gk
      zch <- t(apply(abs(zch),1,norm01))
      ck <- rgb(zch)
      pc <- pc0 * as.vector(odf) 
      pc <- pc / (2*max(pc))
      pos <- 2 * c(i,j,0)
      pcsurf <- cbind(
        pc[tc.surf,1]+pos[1], pc[tc.surf,2]+pos[2], pc[tc.surf,3]+pos[3])
      b <- (k-1)*dt2; e <- b+dt2
      sgrid[(b+1):e, ] <- pcsurf
      vcolors[k,] <- ck[tc.surf]
    }
  }
  cat("Plotting ...\n")
  # rgl.open()
  open3d()
  rgl.viewpoint(theta=0, phi=0)
  rgl.triangles(sgrid[,1], sgrid[,2], sgrid[,3], col=t(vcolors))
  rgl.viewpoint(0,0)
}


#--------------------------
#
# Estimate vMF mixture and principal distribution directions (PDDs)
#
estimate.vmf.lines <-
function(pc0, field, odfsdata, showglyph=FALSE, clusterthr=0.6, savedir=tempdir(), ...)
{
  normvf <- function(x) { norm(matrix(x,length(x),1),"f") }
  ## control parameters for movMF
  E <- list(...)[["E"]]
  if (is.null(E)) E <- "softmax"
  kappa <- list(...)[["kappa"]]
  if (is.null(kappa)) kappa <- "Newton_Fourier"
  minalpha <- list(...)[["minalpha"]]
  if (is.null(minalpha)) minalpha <- 8
  start <- list(...)[["start"]]
  if (is.null(start)) start <- "s"
  startctl=list(E=E, kappa=kappa, minalpha=minalpha, start=start) ## movMF inits
  ## 
  mask <- field$mask
  ## GFA
  odfs.reg <- t(apply(odfsdata, 1 , norm01))
  gfas <- apply(odfs.reg, 1, genfa)
  tc <-  geometry::delaunayn(pc0)
  tc.surf <- t( surf.tri(pc0,tc) )
  ## ------------
  d1 <- dim(odfsdata)[1]
  nn <- 8*d1
  v <- matrix(0, nrow=nn, ncol=3)
  ck <- numeric(nn)
  q <- 1
  m <- 0
  cat("running ... \n")
  ## z2d <- which(mask != 0, arr.ind=TRUE)
  z2d <- which(mask != -1, arr.ind=TRUE)
  lix <- dim(z2d)[1]
  d <- dim(mask)
  ##
  v1perslice <- matrix(0, nrow=lix,ncol=3) # v1 directions 
  v2perslice <- matrix(0, nrow=lix,ncol=3) # v2 directions 
   volgfa <- array(0, dim=c(dim(mask),1)) ## gfas map
  V1 <- array(0, dim=c(dim(mask),1, 3)) ## V1 vol
  V2 <- array(0, dim=c(dim(mask),1, 3)) ## V2 vol
  for(m in 1:lix) {
    odf <- odfs.reg[m,]
    gk <- gfas[m]
    ith <- which(odf < clusterthr) 
    vx <- pc0[-ith,]
    n <- dim(vx)[1]
    nc <- dim(vx)[2]
    kc <- 1:8
    npar <- nc*kc+kc-1
    bic <- -1.0e+10; nf <- 0; yy <- NULL
    for(k in seq(2,4,by=2)) {
      y2 <- movMF::movMF(vx, k=k, control=startctl) 
      par <- logb(n)*npar[k]
      bic2 <- 2*logLik(y2) - par
      if(bic2 > bic) {
        bic <- bic2
        nf <- k
        yy <- y2
      }  
    }
    np <- dim(yy$theta)[1]
    ##   pcoords <- yy$theta/max(yy$theta)
    pk <- list(np=np , pcoords=t(yy$theta)) # no scaling 
    v1perslice[m,] <- pk$pcoords[,1]
    if(np == 4) {
      ## !!! alternative option: use identical alpha values for selection
      if(all.equal(abs(pk$pcoords[,1]), abs(pk$pcoords[,2]), toler=0.01)
         == TRUE) 
        v2perslice[m,] <- pk$pcoords[,3]
      else 
        v2perslice[m,] <- pk$pcoords[,2]
    }
    ## reorder 
    vref <- c(1,0,0)
    c1 <- crossprod(v1perslice[m,], vref)
    if(c1 < 0)
      v1perslice[m,] <- -v1perslice[m,]
    if(np == 4) {
      vref <- c(1,0,0)
      c1 <- crossprod(v1perslice[m,], vref)                
      c2 <- crossprod(v2perslice[m,], vref)                
      if(abs(c2) > abs(c1)) {
        tmp <- v1perslice[m,]
        v1perslice[m,] <- v2perslice[m,]
        v2perslice[m,] <- tmp
      }
    }
    ## V volume
    for(k in 1:3) { # axial
      mx <- matrix(0, d[1],d[2])
      mx[z2d] <- v1perslice[,k]
      V1[,,1,k] <- mx
      mx <- matrix(0, d[1],d[2])
      mx[z2d] <- v2perslice[,k]
      V2[,,1,k] <- mx
      mx <- matrix(0, d[1],d[2])
      mx[z2d] <- gfas
      volgfa[,,1] <- mx
     }
    ##----------------------------
    ## 1st direction of max odf values for  odfs
    pc <- pc0 * as.vector(odf) 
    pc <- pc / (2*max(pc))
    # pos <- c(j,i,0)
    pos <- c(z2d[m,], 0)
    ## normalize for visualization
    mx <- apply(yy$theta, 1, normvf)
    coords <- t(yy$theta/mx)
    for(k in 1:min(pk$np, 4)) {
      zch <- coords[,k] * gk
      zch <- t(norm01(abs(zch)))
      ck[q] <- rgb(zch)
      ck[q+1] <- ck[q]
      pp <- coords[,k]/2
      v[q,] <- pos
      v[q+1,]  <-  pp/2 + pos
      q <- q+2
    }
  }
  cat("\n")
  f <- paste(savedir,"/data_gfa",sep="")
  writeNIfTI(volgfa, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  f <- paste(savedir,"/data_V1",sep="")
  writeNIfTI(V1, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  f <- paste(savedir,"/data_V2",sep="")
  writeNIfTI(V2, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  ##--
  open3d()
  cat("plotting ... \n")
  segments3d(v[1:(q-1),], col=ck[1:(q-1)], lwd=2, alpha=1)
  rgl.viewpoint(0,0)
  rgl.bringtotop()

  ## 
}

##--------------------------
## Built-in field masks
m1 <- function(ang=60) {
  gn <- c(3,3); a <- matrix(ang,gn[1],gn[2]); invisible(a) }
m2 <- function(ang=60) {
  gn <- c(3,3); a <- matrix(-1, gn[1], gn[2]);
  a[,1] <- 180-ang; a[,2] <- 90; a[,3] <- ang;
  invisible(a) }
m3 <- function(ang=60) {
  gn <- c(5,5); a <- matrix(-1, gn[1], gn[2])
  a[,1] <- 120; a[,2] <- 120; a[,3] <- 90;
  a[,4] <- ang; a[,5] <- ang
  invisible(a) }
mx1 <- function(ang) {
  gn <- c(4,4); a <- array(-1, dim=c(gn[1],gn[2],2))
  a[1,3:4,1] <- ang; a[2,1:4,1] <- 0; a[3,1:4,1] <- 0
  a[4,1:2,1] <- ang; a[2,2:3,2] <- ang; a[3,2:3,2] <- ang
  invisible(a) }
mx2 <- function(ang=90) {
  gn <- c(13,13)
  a <- array(-1, dim=c(gn[1],gn[2],2))
  a[6:8,,1] <- 0; 
  a[1:5,6:8,1] <- ang; 
  a[1:5,6:8,1] <- ang; 
  a[9:13,6:8,1] <- ang; 
  a[6:8,6:8,2] <- ang; 
  invisible(a) }
mx3 <- function(ang=45) {
  gn=c(10,10)
  a <- array(-1, dim=c(gn[1],gn[2],2))
  a[4:7,1:6,1] <- 0; 
  a[5:6,,1] <- 0; 
  a[1,8:9,1] <- ang; 
  a[2,7:8,1] <- ang; 
  a[3,6:7,1] <- ang; 
  a[8,6:7,1] <- 360-ang
  a[9,7:8,1] <- 360-ang
  a[10,8:9,1] <- 360-ang
  a[4,5:6,2] <- ang
  a[5, 4:5,2] <- ang
  a[6,4:5,2] <- 360-ang
  a[7,5:6,2] <- 360-ang
  invisible(a) }
#--------------------------

