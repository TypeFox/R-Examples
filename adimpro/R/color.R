colorize <- function(obj, compress=TRUE) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type %in% c("greyscale","grayscale")) {
     img <- extract.image(obj)
     img <- array(c(img,img,img),c(obj$dim,3))
     obj <- make.image(img,gammatype=obj$gammatype,compress=compress)
}
obj
}

rgb2grey <- function(obj, compress=TRUE) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (!(obj$type %in% c("rgb","xyz","yuv","yiq"))) {
    warning("Error: image type is not implemented \n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    if(obj$gamma) obj <- invgamma.correction(obj,alg=1)
    dm <- obj$dim
    dim(obj$img) <- c(prod(dm),3)
    if(obj$cspace!="sRGB") {
       obj$img <- obj$img%*%t(xyz2rgbmat("sRGB")%*%rgb2xyzmat(obj$cspace))
       obj$cspace <- "greyscale"
       obj$img <- as.integer(pmax(0,pmin(65535,c(0.299, 0.587, 0.114) %*% t(obj$img))))
    } else {
    obj$img <- as.integer(c(0.299, 0.587, 0.114) %*% t(obj$img))
    }
    dim(obj$img) <- dm
    obj$type <- "greyscale"
  }
  invisible(if(compress) compress.image(obj) else obj)
}

rgb2yiq <- function(obj) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (!(obj$type %in% c("rgb","xyz","yuv","yiq"))) {
    warning("Error: image type is not implemented \n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    if(obj$gamma) obj <- invgamma.correction(obj,alg=1)
    dm <- obj$dim
    obj$img <- obj$img / 65535
    dim(obj$img) <- c(prod(dm),3)
    
    conv <- c(0.299,  0.595716,  0.211456,
              0.587, -0.274453, -0.522591,
              0.114, -0.321263,  0.311135)
    dim(conv) <- c(3,3)
    if(obj$cspace!="sRGB") {
       conv <- conv%*%xyz2rgbmat("sRGB")%*%rgb2xyzmat(obj$cspace)
    }
    
    obj$img <- obj$img %*% t(conv)
    dim(obj$img) <- c(dm,3)
    obj$cspace <- "yiq"
    obj$type <- "yiq"
  }
  invisible(obj)
}

yiq2rgb <- function(obj, cspace="Adobe", compress=TRUE) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="yiq") {
    warning("Error: image type is not yiq\n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    dim(obj$img) <- c(prod(dm),3)
    
    conv <- c(1,          1,          1,
              0.9562957, -0.2721221, -1.1069890,
              0.6210244, -0.6473806,  1.7046150)
    dim(conv) <- c(3,3)
    if(cspace!="sRGB") {
       conv <- xyz2rgbmat(cspace)%*%rgb2xyzmat("sRGB")%*%conv
    }
    
    obj$img <- as.integer(65535 * pmax(0, pmin(1, obj$img %*% t(conv))))
    dim(obj$img) <- c(dm,3)
    obj$cspace <- cspace
    if(cspace %in% c("sRGB","Adobe","wGamut","kodak")) obj$type <- "rgb" else obj$type <- cspace
  }
  invisible(if(compress) compress.image(obj) else obj)
}

rgb2yuv <- function(obj) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (!(obj$type %in% c("rgb","xyz","yuv","yiq"))) {
    warning("Error: image type is not implemented \n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    if(obj$gamma) obj <- invgamma.correction(obj,alg=1)
    dm <- obj$dim
    obj$img <- obj$img / 65535
    dim(obj$img) <- c(prod(dm),3)
    
    conv <- c(0.299, -0.147,  0.615,
              0.587, -0.289, -0.515,
              0.114,  0.436, -0.100)
    dim(conv) <- c(3,3)
    if(obj$cspace!="sRGB") {
       conv <- conv%*%xyz2rgbmat("sRGB")%*%rgb2xyzmat(obj$cspace)
    }
    
    obj$img <- obj$img %*% t(conv)
    dim(obj$img) <- c(dm,3)
    obj$cspace <- "yuv"
    obj$type <- "yuv"
  }
  invisible(obj)
}

rgb2yuvhist <- function(obj) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (!(obj$type %in% c("rgb","xyz","yuv","yiq"))) {
    warning("Error: image type is not implemented \n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    obj$img <- obj$img / 65535
    dim(obj$img) <- c(prod(dm),3)
    
    conv <- c(0.299, -0.147,  0.615,
              0.587, -0.289, -0.515,
              0.114,  0.436, -0.100)
    dim(conv) <- c(3,3)
    if(obj$cspace!="sRGB") {
       conv <- conv%*%xyz2rgbmat("sRGB")%*%rgb2xyzmat(obj$cspace)
    }
    
    obj$img <- obj$img %*% t(conv)
    dim(obj$img) <- c(dm,3)
    obj$cspace <- "yuv"
    obj$type <- "yuv"
  }
  invisible(obj)
}

yuv2rgb <- function(obj, cspace="Adobe", compress=TRUE) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="yuv") {
    warning("Error: image type is not yuv\n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    dim(obj$img) <- c(prod(dm),3)
    
    conv <- c( 1,         1,         1,
              -0.000039, -0.394610,  2.032000,
              1.139828, -0.580500, -0.000481)
    dim(conv) <- c(3,3)
    if(cspace!="sRGB") {
       conv <- xyz2rgbmat(cspace)%*%rgb2xyzmat("sRGB")%*%conv
    }
    
    obj$img <- as.integer(65535 * pmax(0, pmin(1, obj$img %*% t(conv))))
    dim(obj$img) <- c(dm,3)
    obj$cspace <- cspace
    if(cspace %in% c("sRGB","Adobe","wGamut","kodak")) obj$type <- "rgb" else obj$type <- cspace
  }
  invisible(if(compress) compress.image(obj) else obj)
}

hsi2rgb <- function (obj, cspace="Adobe", compress=TRUE) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="hsi") {
    warning("Error: image type is not hsi\n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    proddm <- prod(dm)
    dim(obj$img) <- c(proddm, 3)
    
    obj$img[, 1] <- obj$img[, 1] * 2 * pi
    obj$img[obj$img[,1]<0,1] <- 0
    obj$img[obj$img[,1]>2*pi,1] <- 2*pi
    
    r <- rep(0, proddm)
    g <- rep(0, proddm)
    b <- rep(0, proddm)
    
    ind1<-(1:proddm)[(obj$img[,1] < 2 * pi/3)]
    ind2<-(1:proddm)[(obj$img[,1] >= 2 * pi/3) & (obj$img[,1] < 4 * pi/3)] # use & operator here to compare vectors!! not &&
    ind3<-(1:proddm)[(obj$img[,1] >= 4 * pi/3)]
    
    b[ind1] <- obj$img[ind1,3] * (1 - obj$img[ind1,2])
    r[ind1] <- obj$img[ind1,3] * (1 + obj$img[ind1,2] * cos(obj$img[ind1,1])/cos(pi/3 - obj$img[ind1,1]))
    g[ind1] <- 3 * obj$img[ind1,3] - (b[ind1] + r[ind1])
    
    r[ind2] <- obj$img[ind2,3] * (1 - obj$img[ind2,2])
    g[ind2] <- obj$img[ind2,3] * (1 + obj$img[ind2,2] * cos(obj$img[ind2,1] - 2 * pi/3)/cos(pi - obj$img[ind2,1]))
    b[ind2] <- 3 * obj$img[ind2,3] - (r[ind2] + g[ind2])
    
    g[ind3] <- obj$img[ind3,3] * (1 - obj$img[ind3,2])
    b[ind3] <- obj$img[ind3,3] * (1 + obj$img[ind3,2] * cos(obj$img[ind3,1] - 4 * pi/3)/cos(5 * pi/3 - obj$img[ind3,1]))
    r[ind3] <- 3 * obj$img[ind3,3] - (g[ind3] + b[ind3])
    
    if(cspace!="sRGB") {
       conv <- xyz2rgbmat(cspace)%*%rgb2xyzmat("sRGB")
       obj$img <- cbind(r, g, b) %*% t(conv)
       if(cspace %in% c("sRGB","Adobe","wGamut","kodak")) obj$img <- as.integer(65535 * pmax(0, pmin(1,obj$img )))
    } else {
       obj$img <- as.integer(65535 * pmax(0, pmin(1, cbind(r, g, b))))    
    }
    dim(obj$img) <- c(dm, 3)
    obj$cspace <- cspace
    if(cspace %in% c("sRGB","Adobe","wGamut","kodak")) obj$type <- "rgb" else obj$type <- cspace
  }
  invisible(if(compress) compress.image(obj) else obj)
}

rgb2hsi <- function (obj) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (!(obj$type %in% c("rgb","xyz","yuv","yiq"))) {
    warning("Error: image type is not implemented \n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    if(obj$gamma) obj <- invgamma.correction(obj,alg=1)
    dm <- obj$dim
    if(obj$cspace %in% c("sRGB","Adobe","wGamut","kodak")) obj$img <- obj$img / 65535
    dim(obj$img) <- c(prod(dm), 3)
    if(obj$cspace!="sRGB") {
       conv <- xyz2rgbmat("sRGB")%*%rgb2xyzmat(obj$cspace)
#       if(obj$type != "rgb") conv <- conv*65535
       obj$img <- obj$img%*%t(conv)
       obj$cspace <- "sRGB"
    }
    
    z <- obj$img[,1] - 0.5 * obj$img[,2] - 0.5 * obj$img[,3]
    n <- (obj$img[,1] - obj$img[,2])^2 + (obj$img[,1] - obj$img[,3]) * (obj$img[,2] - obj$img[,3])
    frac <- z[n != 0]/sqrt(n[n != 0])
    frac <- signif(frac, 5)
    
    h <- rep(0, prod(dm))
    h[n != 0] <- acos(frac) # h not defined for n==0
    h[obj$img[, 3] > obj$img[,2]] <- 2 * pi - h[obj$img[,3] > obj$img[,2]]
    h <- h/(2 * pi)
    
    i <- (obj$img[,1] + obj$img[,2] + obj$img[,3])/3
    i <- pmax(0,pmin(1,i))
    s <- rep(0, prod(dm))
    s[i != 0] <- 1 - 1/i[i != 0] * pmin(obj$img[,1], obj$img[,2], obj$img[,3])[i != 0] # s not defined for i==0
    s <- pmax(0,pmin(1,s))
    obj$img <- c(h, s, i)
    dim(obj$img) <- c(dm, 3)
    obj$cspace <- "hsi"
    obj$type <- "hsi"
  }
  invisible(obj)
}



rgb2xyzmat <- function(cspace){
#
#  provides transfer matrix from RGB to XYZ
#
z <- matrix(switch(EXPR=cspace,
              sRGB=c(0.412424,0.212656,0.0193324,
                     0.357579,0.715158,0.119193,
                     0.180464,0.0721856,0.950444),
              Adobe=c(0.576700,0.297361,0.0270328,
                      0.185556,0.627355,0.0706879,
                      0.188212,0.0752847,0.991248),
              wGamut=c(0.716105,0.258187,0.000000,
                       0.100930,0.724938,0.0517813,
                       0.147186,0.0168748,0.773429),
              kodak=c(0.797675,0.288040,0.000000,
                      0.135192,0.711874,0.000000,
                      0.0313534,0.000086,0.825210),
              xyz=diag(c(1,1,1)),
              c(0.412424,0.212656,0.0193324,
                     0.357579,0.715158,0.119193,
                     0.180464,0.0721856,0.950444)),3,3)
if(cspace %in% c("yuv","yiq")) {
z <- z%*%matrix(switch(EXPR=cspace,
        yuv=c( 1,         1,         1,
              -0.000039, -0.394610,  2.032000,
              1.139828, -0.580500, -0.000481),
        yiq=c(1,          1,          1,
              0.9562957, -0.2721221, -1.1069890,
              0.6210244, -0.6473806,  1.7046150)),3,3)
}
z
}

xyz2rgbmat <- function(cspace){
#
#  provides transfer matrix from XYZ to RGB
#
z <- matrix(switch(EXPR=cspace,
              sRGB=c(3.24071,-0.969258,0.0556352,
                     -1.53726,1.87599,-0.203996,
                     -0.498571,0.0415557,1.05707),
              Adobe=c(2.04148,-0.969258,0.0134455,
                      -0.564977,1.87599,-0.118373,   
                      -0.344713,0.0415557,1.01527),
              wGamut=c(1.46281,-0.521793,0.0349342,
                       -0.184062,1.44724,-0.0968931,
                       -0.274361,0.0677228,1.28841),
              kodak=c(1.34594,-0.544599,0.000000,
                      -0.255608,1.50817,0.000000,
                      -0.0511118,0.0205351,1.21181),
              xyz=diag(c(1,1,1)),
              c(3.24071,-0.969258,0.0556352,
                     -1.53726,1.87599,-0.203996,
                     -0.498571,0.0415557,1.05707)),3,3)
if(cspace %in% c("yuv","yiq")) {
z <- matrix(switch(EXPR=cspace,
        yuv=c(0.299, -0.147,  0.615,
              0.587, -0.289, -0.515,
              0.114,  0.436, -0.100),
        yiq=c(0.299,  0.595716,  0.211456,
              0.587, -0.274453, -0.522591,
              0.114, -0.321263,  0.311135)),3,3)%*%z
}
z
}
rgb2xyz <- function(obj) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type %in% c("rgb","xyz","yuv","yiq")) {
     if(obj$compressed) obj <- decompress.image(obj)
     if(obj$gamma) obj <- invgamma.correction(obj,alg=1)
     if(obj$type=="rgb") obj$img <- obj$img / 65535
     if (obj$cspace!="xyz") {
#      if object is already in XYZ not much to do ... 
       dm <- obj$dim
       dim(obj$img) <- c(prod(dm),3)
    
       conv <- rgb2xyzmat(obj$cspace)
    
       obj$img <- obj$img %*% t(conv)
       obj$img[obj$img>1] <- 1
       obj$img[obj$img<0] <- 0
       dim(obj$img) <- c(dm,3)
    } else{
    warning("Error: incorrect image type in rgb2xyz \n")
    } 
    obj$cspace <- "xyz"
    obj$type <- "xyz"
  }
  invisible(obj)
}

xyz2rgb <- function(obj, cspace="Adobe", black=0, exposure=1, compress=TRUE) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (!(cspace%in%c("sRGB","Adobe","wGamut","kodak","xyz","yuv","yiq","hsi"))){
     warning(paste("invalid color space",cspace,"reset to Adobe"))
     cspace <- "Adobe"
     } 
  rgb <- cspace%in%c("sRGB","Adobe","wGamut","kodak")
  if (obj$type!="xyz") {
    warning("Error: image type is not xyz\n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    dim(obj$img) <- c(prod(dm),3)
    if(black!=0) {
       y <- obj$img[,2]
       obj$img[,1] <- obj$img[,1]/y*(y - black)
       obj$img[,2] <- y - black
       obj$img[,3] <- obj$img[,3]/y*(y - black)
       obj$img[is.na(obj$img)|obj$img<0] <- 0
    }
    conv <- xyz2rgbmat(cspace)*exposure
    
    if(rgb) {
       obj$img <- as.integer(65535 * pmax(0, pmin(1, obj$img %*% t(conv))))
       } else {
       obj$img <- obj$img %*% t(conv)
       }
    dim(obj$img) <- c(dm,3)
    obj$cspace <- cspace
    obj$type <- switch(cspace,"sRGB"="rgb","Adobe"="rgb","wGamut"="rgb","kodak"="rgb",cspace)
  }
  invisible(if(compress&&rgb) compress.image(obj) else obj)
}

cam2rgbmat <- function(obj, cspace="sRGB") {
  if(obj$compressed) obj <- decompress.image(obj)
  cam <- extract.info(obj,"Camera")
  cam.xyz <- raw2xyzmat(obj)
  xyz.rgb <- rgb2xyzmat(cspace)
  cam.rgb <- t(cam.xyz)%*%xyz.rgb
  out.cam <- t(solve(cam.rgb/apply(cam.rgb,1,sum)))
  out.cam
}

gamma.correction <- function (img, gammatype="ITU",
                              nbins = 65536, alg = 1, log = FALSE) {
  if(gammatype=="histogram") return(hequalize(img,compress=FALSE))
  ga <- switch(gammatype,None=1,ITU=20/9,sRGB=2.4,CIE=3,1)
  bp <- switch(gammatype,None=0,ITU=0.018,sRGB=0.00304,CIE=0.008856,0)
  sls <- 1/(ga/bp^(1/ga - 1) - ga * bp + bp)
  fs <- ga/bp^(1/ga - 1)    # slope divided by sls to easy computation
  c0 <- fs * bp^(1/ga) - bp # segment offset divided by sls to easy computation
  if (log) {
    cat("Gamma correction - Gamma:",ga,"\n")
    cat("             Break Point:",bp,"\n")
    cat("                   Slope:",sls,"\n")
    cat("   Slope matching factor:",sls * fs,"\n")
    cat("          Segment offset:",sls * c0,"\n")
  }
  
  if(img$type=="rgb"||img$type=="greyscale") img$img <- img$img / 65535 else {
     for(i in 1:dim(img$img)[3])
     img$img[,,i] <- (img$img[,,i]-min(img$img[,,i]))/(max(img$img[,,i])-min(img$img[,,i]))
  }
  di <- dim(img$img)
  
  breaks <- seq(0,1+1/(nbins-1),length=nbins+1)
  nimg <- (nbins-1)*img$img
  iimg <- as.integer(nimg)
  ind <- breaks>bp
  if (alg == 1) {
    midpoints <- (breaks[-1]+breaks[-nbins-1])/2
    gammaofmidpoints <- sls * midpoints
    gammaofmidpoints[ind] <- sls * (fs * midpoints[ind]^(1/ga) - c0)
    img$img <- gammaofmidpoints[iimg+1]
  } else if (alg == 2) {
    aimg <- nimg - iimg
    gammaofbreaks <- sls * breaks
    gammaofbreaks[ind] <- sls * (fs * breaks[ind]^(1/ga) - c0)
    img$img <- (1-aimg)*gammaofbreaks[iimg+1] + aimg*gammaofbreaks[iimg+2]
  } else {
    ind <- (1:length(img$img))[img$img > bp]
    img$img[ind] <- fs * img$img[ind]^(1/ga) - c0
    img$img <- sls * img$img
  }
  
  img$img <- as.integer(65535 * img$img)
if(any(img$img>65535)) cat("large values in gamma\n")
  dim(img$img) <- di
  if(ga!=1) img$gamma <- TRUE
  img$gammatype <- gammatype
  storage.mode(img$img) <- "integer"
  invisible(img)
}

invgamma.correction <- function (img,
                              nbins = 65536, alg = 1, log = FALSE) {
 compress <- img$compressed
 if(img$depth=="8Bit") warning("inverting gamma may lead to image degradation for 8 Bit images")
 if(img$gammatype=="histogram") return(invhequalize(img))
 if(img$compressed) img <- decompress.image(img)
  ga <- switch(img$gammatype,None=1,ITU=20/9,sRGB=2.4,CIE=3,1)
  bp <- switch(img$gammatype,None=0,ITU=0.018,sRGB=0.00304,CIE=0.008856,0)

  sls <- 1/(ga/bp^(1/ga - 1) - ga * bp + bp)
  fs <- ga/bp^(1/ga - 1)    # slope divided by sls to easy computation
  c0 <- fs * bp^(1/ga) - bp # segment offset divided by sls to easy computation
  
  img$img <- img$img / 65535
  di <- dim(img$img)
    ind <- (1:length(img$img))[img$img > sls*bp]
    img$img[ind] <- ((img$img[ind]/sls + c0)/fs)^ga
    img$img[-ind] <- img$img[-ind]/sls
  
  img$img <- as.integer(65535 * img$img)
if(any(img$img>65535)) cat("large values in gamma\n")
  dim(img$img) <- di
  if(ga!=1) {
     img$gamma <- FALSE
     img$gammatype <- "None"
  }
 invisible(if(compress)  compress.image(img) else img)
 }

whitepoint <- function(wp){
#
#  evaluate white point entries to get corresponding x,y in xyY
# Adapted from Wikipedia
# note that www.brucelindbloom.com gives slightly
# different values
if(is.character(wp)){
wp <- switch(EXPR=wp,E=c(1,1)/3,
                D50=c(.34567,.35850),
                D55=c(.33342,.34743),
                D65=c(.31271,.32902),
                D75=c(.29902,.31740),
                A=c(.44757,.40745),
                B=c(.34842,.35161),
                C=c(.31006,.31616),
                F2=c(.37207,.37512),
                F7=c(.31285,.32918),
                F11=c(.38054,.37691),
                NULL)
} else if(!is.numeric(wp)||length(wp)!=2) wp <- NULL
if(is.null(wp)) {
warning("Incorrect whitepoint, set to D65")
wp <- c(.31271,.32902)
}
wp
}

wpofT <- function(temp){
# Approximate whitepoint coordinates for given color temperature 
#    
if(!is.numeric(temp)) {
warning("Non-numeric color temperatur. Using 6500 K")
temp <- 6500
}
temp <- min(10000,max(1700,temp))
cx <- c(8.429971e-01, -2.206539e-01,  3.663518e-02, -2.921696e-03,  9.048406e-05)
cy <- c(0.7164137115, -0.4434534157, 0.0494384240, -0.0031368040, 0.0000825464, 0.5974247702)
tx <- (temp/1000)^(0:4)
ty <- c(tx,log(temp/1000))
c(sum(cx*tx),sum(ty*cy))
}

changewhitepoint <- function(wsource,wdest,kind="Bradford"){
#
#  provides transfer matrix from XYZ(wsource) to XYZ(wdest) 
#
MA <- matrix(switch(EXPR=kind,
             Bradford=c(0.8951,-0.7502,0.0389,
                        0.2664,1.7135,-0.0685,
                        -0.1614,0.0367,1.0296),
             XYZscaling=diag(c(1,1,1)),
             VonKries=c(0.40024,-0.22630,0.00000,
                        0.70760,1.16532,0.00000, 
                       -0.08081,0.04570,0.91822)),3,3)
MAinv <- matrix(switch(EXPR=kind,
             Bradford=c(0.986993,0.432305,-0.008529,
                        -0.147054,0.518360,0.040043,
                        0.159963,0.049291,0.968487),
             XYZscaling=diag(c(1,1,1)),
             VonKries=c(1.859936,0.361191,0.000000,
                        -1.129382,0.638812,0.000000,
                        0.219897,-0.000006,1.089064)),3,3)
ws <- c(wsource,1-sum(wsource))%*%MA
wd <- c(wdest,1-sum(wdest))%*%MA
MA%*%diag(as.vector(wd/ws))%*%MAinv
}

white.balance <- function(img,wdest,kind="Bradford"){
mat <- changewhitepoint(whitepoint(img$whitep),whitepoint(wdest),kind=kind)
dm <- img$dim
dim(img$img) <- c(prod(dm[1:2]),3)
img$img <- img$img%*%mat
dim(img$img) <- c(dm,3)
img$whitep <- wdest
img
}

adjust.image <- function(img, gammatype=NULL, cspace=NULL, whitep=NULL, temp=NULL, black=0, exposure=1, kind="Bradford", alg = 1, compress=TRUE) {
  
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if(abs(black)>1) black <- black/65535
  if(img$compressed) img <- decompress.image(img)
  if(is.null(cspace)) cspace <- img$cspace
  if(cspace %in% c("xyz","yuv","yiq","hsi")) {
     gammatype <- "None"
     rgb <- FALSE
  }
  dogamma <- !is.null(gammatype)&&img$gammatype!=gammatype
  if(img$type=="greyscale"){
     if(exposure<.1) {
        warning("misspecified value of exposure, set value to 1")
        exposure <- 1
     }
     if(dogamma||exposure!=1||black!=0){
        if(is.null(gammatype)) gammatype <- img$gammatype
        gamma <- gammatype!="None"
        if(img$gamma) img <- invgamma.correction(img,alg=alg)
        if(black!=0) img$img <- as.integer(pmax(0,img$img-black))
        if(exposure!=1) img$img <- as.integer(pmin(65535,img$img*exposure)) 
        dim(img$img) <- img$dim
        if(gamma) img <- gamma.correction(img,gammatype=gammatype,alg=alg)
     }
  } else {
  if(is.null(whitep)) if(is.null(temp)) whitep <- img$whitep else whitep <- wpofT(temp)
  dowhite <- !is.null(whitep)&&any(whitepoint(img$whitep)!=whitepoint(whitep))
  docspace <- !is.null(cspace)&&img$cspace!=cspace
  dimg <- dim(img$img)
  if(exposure<.1) {
     warning("misspecified value of exposure, set value to 1")
     exposure <- 1
  }
  if(dogamma||dowhite||docspace||black!=0){
     if(is.null(gammatype)) gammatype <- img$gammatype
     gamma <- gammatype!="None"
     if(img$gamma) img <- invgamma.correction(img,alg=alg)
#  transfer to XYZ
     if(dowhite||docspace||black!=0){
     if(img$type=="hsi"){
        img <- hsi2rgb(img,"xyz",compress=FALSE)
     } else if(img$cspace!="xyz") img <- rgb2xyz(img)
     if(black!=0){
       y <- img$img[,,2]
       img$img[,,1] <- img$img[,,1]/y*(y - black)
       img$img[,,2] <- y - black
       img$img[,,3] <- img$img[,,3]/y*(y - black)
       img$img[is.na(img$img)|img$img<0] <- 0
     }
     if(dowhite) img <- white.balance(img,whitep,kind)
     if(cspace=="hsi"){
        if(exposure!=1) img$img <- img$img*exposure
        img <- rgb2hsi(img)
     } else if(cspace %in% c("grey","gray","grayscale")) {
        if(exposure!=1) img$img <- img$img*exposure
        img <- rgb2grey(img)
     } else img <- xyz2rgb(img,cspace=cspace,exposure=exposure,compress=FALSE)
     }
     if(gamma) img <- gamma.correction(img,gammatype=gammatype,alg=alg)
  } else if(exposure!=1) {
    img$img <- as.integer(img$img*exposure)
    img$img[img$img>65535] <- as.integer(65535)
    dim(img$img) <- dimg
  }
  }
  if(img$type %in% c("greyscale","rgb")) storage.mode(img$img) <- "integer"
  invisible(if(compress) compress.image(img) else img)
}

combine <- function(img1,img2,fun="+",rescale=TRUE,compress=TRUE, gammatype="None", whitep = "D65", cspace="Adobe",xmode="RGB",...){
  ff <- match.fun(fun)
  z <- formals(ff)
  if(!(length(z) == 0||length(z) == 2|| length(z) == 3)) {
    stop(" function specified in argument fun needs two arguments.\n")
  }
  if(!check.adimpro(img1)) {
    stop(" Consistency check for argument img1 failed (see warnings).\n")
  }
  if(!check.adimpro(img2)) {
    stop(" Consistency check for argument img2 failed (see warnings).\n")
  }
  if(any(img1$dim != img2$dim)) {
     stop(" images need to have same dimension .\n")
  }
  img1 <- adjust.image(img1,gammatype="None", whitep = "D65", cspace="Adobe")
  img2 <- adjust.image(img2,gammatype="None", whitep = "D65", cspace="Adobe")
  if(img1$cspace!=img2$cspace) {
        warning("Combining Color and greyvalued image \n")
        if(img1$cspace=="greyscale"){
           img1 <- colorize(img1)
        }
        if(img2$cspace=="greyscale"){
           img2 <- colorize(img2)
        }
  }
  if(length(z) == 0||length(z) == 2){
  img3 <- ff(extract.image(img1),extract.image(img2))
  } else {
  img3 <- ff(extract.image(img1),extract.image(img2),...)
  }
  if(rescale){
  img3 <- 65535*(img3-min(img3))/(max(img3)-min(img3))
  } else {
  img3[img3<0] <- 0
  img3[img3>65535] <- 65535
  }
  if(img1$cspace=="greyscale"&&img2$cspace=="greyscale"){
     dim(img3) <- img1$dim
     cpace <- "greyscale"
  } else {
     dim(img3) <- c(img1$dim,3)
  }
  make.image(img3,compress=compress, gammatype=gammatype, whitep = whitep,
             cspace="Adobe", scale="Combined",xmode=xmode)
}

hequalize <- function(img,compress=TRUE){
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if(img$compress) img <- decompress.image(img)
  gammatype <- img$gammatype
  if(img$gamma) img <- invgamma.correction(img)
  if(img$type=="greyscale"){
     n1 <- img$dim[1]
     n2 <- img$dim[2]
     z <- .Fortran("hequalg",
                    as.integer(img$img),
                    as.integer(n1*n2),
                    img=integer(n1*n2),
                    cumhist=integer(65536),
                    DUPL=FALSE,
                    PACKAGE="adimpro")[c("img","cumhist")]
     img$img <- matrix(z$img,n1,n2)
     img$hequal <- z$cumhist
  } else {
     type <- img$type
     cspace <- img$cspace
     if(type%in%c("rgb","xyz")){
        img1 <- rgb2grey(img,compress=FALSE)
        n1 <- img1$dim[1]
        n2 <- img1$dim[2]
        cumhist <- .Fortran("cumhist",
                    as.integer(img1$img),
                    as.integer(n1*n2),
                    cumhist=integer(65536),
                    DUPL=FALSE,
                    PACKAGE="adimpro")$cumhist
        img$img <- array(.Fortran("hequalc",
                    as.integer(img$img),
                    as.integer(n1*n2),
                    img=integer(n1*n2*3),
                    as.integer(cumhist),
                    DUPL=FALSE,
                    PACKAGE="adimpro")$img,c(n1,n2,3))
        img$hequal <- cumhist
     } else {
      warning(paste("not yet implemented for type",type)) 
     }
  }
  img$gamma <- TRUE
  img$gammatype <- "histogram"
#  if(gammatype!="None") img <- gamma.correction(img,gammatype=gammatype)
  if(compress) compress.image(img) else img
}
hequalize.old <- function(img,compress=TRUE){
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if(img$compress) img <- decompress.image(img)
  gammatype <- img$gammatype
  if(img$gamma) img <- invgamma.correction(img)
  if(img$type=="greyscale"){
     n1 <- img$dim[1]
     n2 <- img$dim[2]
     z <- .Fortran("hequalg",
                    as.integer(img$img),
                    as.integer(n1*n2),
                    img=integer(n1*n2),
                    cumhist=integer(65536),
                    DUPL=FALSE,
                    PACKAGE="adimpro")[c("img","cumhist")]
     img$img <- matrix(z$img,n1,n2)
     img$hequal <- z$cumhist
  } else {
     type <- img$type
     cspace <- img$cspace
     if(type%in%c("rgb","xyz")){
        img <- rgb2yuv(img)
        n1 <- img$dim[1]
        n2 <- img$dim[2]
        z <- .Fortran("hequalg",
                    as.integer(img$img[,,1]*65535),
                    as.integer(n1*n2),
                    img=integer(n1*n2),
                    cumhist=integer(65536),
                    DUPL=FALSE,
                    PACKAGE="adimpro")[c("img","cumhist")]
        img$img[,,1] <- z$img/65535
        img$hequal <- z$cumhist
        img <- yuv2rgb(img,cspace=cspace,compress=FALSE)
     } else {
      warning(paste("not yet implemented for type",type)) 
     }
  }
  img$gamma <- TRUE
  img$gammatype <- "histogram"
#  if(gammatype!="None") img <- gamma.correction(img,gammatype=gammatype)
  if(compress) compress.image(img) else img
}

invhequalize <- function(img){
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if(img$gammatype!="histogram") {
    stop(" gammatype != 'histogram'.\n")
  }
  hist <- img$hequal
  if(img$compress) img <- decompress.image(img)
  if(img$type=="greyscale"){
     n1 <- img$dim[1]
     n2 <- img$dim[2]
     img$img <- matrix(.Fortran("ihequal",
                    as.integer(img$img),
                    as.integer(n1*n2),
                    img=integer(n1*n2),
                    as.integer(hist),
                    DUPL=FALSE,
                    PACKAGE="adimpro")$img,n1,n2)
     img$hequal <- NULL
     img$gamma <- FALSE
     img$gammatype <- "None"
  } else {
     type <- img$type
     cspace <- img$cspace
     if(type%in%c("rgb","xyz")){
        n1 <- img$dim[1]
        n2 <- img$dim[2]
        img$img <- array(.Fortran("ihequalc",
                    as.integer(img$img),
                    as.integer(n1*n2),
                    img=integer(n1*n2*3),
                    as.integer(hist),
                    DUPL=FALSE,
                    PACKAGE="adimpro")$img,c(n1,n2,3))
        img$hequal <- NULL
        img$gamma <- FALSE
        img$gammatype <- "None"
     } else {
      warning(paste("not yet implemented for type",type)) 
     }
  }
  img
}
invhequalize.old <- function(img){
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if(img$gammatype!="histogram") {
    stop(" gammatype != 'histogram'.\n")
  }
  hist <- img$hequal
  if(img$compress) img <- decompress.image(img)
  if(img$type=="greyscale"){
     n1 <- img$dim[1]
     n2 <- img$dim[2]
     img$img <- matrix(.Fortran("ihequal",
                    as.integer(img$img),
                    as.integer(n1*n2),
                    img=integer(n1*n2),
                    as.integer(hist),
                    DUPL=FALSE,
                    PACKAGE="adimpro")$img,n1,n2)
     img$hequal <- NULL
     img$gamma <- FALSE
     img$gammatype <- "None"
  } else {
     type <- img$type
     cspace <- img$cspace
     if(type%in%c("rgb","xyz")){
        img <- rgb2yuvhist(img)
#   no need to revers the gamma correction here,
#   so we can't use rgb2yuv
        n1 <- img$dim[1]
        n2 <- img$dim[2]
        img$img[,,1] <- .Fortran("ihequal",
                    as.integer(img$img[,,1]*65535),
                    as.integer(n1*n2),
                    img=integer(n1*n2),
                    as.integer(hist),
                    DUPL=FALSE,
                    PACKAGE="adimpro")$img/65535
        img$hequal <- NULL
        img$gamma <- FALSE
        img$gammatype <- "None"
        img <- yuv2rgb(img,cspace=cspace,compress=FALSE)
     } else {
      warning(paste("not yet implemented for type",type)) 
     }
  }
  img
}