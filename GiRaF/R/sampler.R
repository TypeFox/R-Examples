exact.mrf <- function(h, w, param, ncolors = 2, nei = 4, pot = NULL, 
                      top = NULL, left = NULL, bottom = NULL, 
                      right = NULL, corner = NULL, view = FALSE){
  
  # checking parameters
  if(missing(h)){
    stop('argument "h" is missing with no default')
  }
  if(missing(w)){
    stop('argument "w" is missing with no default')
  }
  if(w < h){
    stop('incorrect arguments values. "h <= w" is mandatory')
  }
  if(ncolors^h > 2^19){
    stop('incorrect argument value. "ncolors^h < 33 554 432" is mandatory')
  }
  
  if(missing(param)){
    stop('argument "param" is missing with no default')
  } else {
    if (nei == 4){
      if (length(param) != 1 & length(param)!= 2){
        stop('incorrect argument value. "param" has to be scalar (isotropic interaction)
             or of size 2 (anisotropic interaction)')
      }
    } else {
      if (length(param) != 1 & length(param)!= 4){
        stop('incorrect argument value. "param" has to be scalar (isotropic interaction)
               or of size 4 (anisotropic interaction)')
      }
    }
  }
  if (nei != 4 & nei != 8){
    stop('incorrect argument values. Allowed values for "nei" are 4 or 8')
  }
  
  if (nei == 4){
    param <- param*c(1., 1., 0., 0.)
  } else {
    param <- param*c(1., 1., 1., 1.)
  }
  
  # Init. the graph and set the borders
  if (is.null(pot)){
    myBlock <- new(Block, h, w, ncolors, nei, param)
  } else {
    if (length(pot) != ncolors){
      stop('incorrect size for argument "pot". Length of "pot" should be "ncolors" = ', ncolors)
    }
    myBlock <- new(Block, h, w, ncolors, nei, param, pot)
  }
  myBlock$initFactor()
  
  if(is.null(top) & is.null(left) & is.null(bottom) & is.null(right) & is.null(corner)){
    
    Z <- myBlock$recursion.mem()
    myBlock$sample(Z)
    
  } else {
    # setting artificial borders if needed
    if (is.null(top)){
      top <- as.vector(rep(ncolors,w))
    } else {
      if (length(top) != w){
        stop('incorrect size for argument "top". Length of "top" should be "w" = ',w)
      }
    }
    if (is.null(left)){
      left <- as.vector(rep(ncolors,h))
    } else {
      if (length(left) != h){
        stop('incorrect size for argument "left". Length of "left" should be "h" = ', h)
      }
    }
    if (missing(bottom)){
      bottom <- as.vector(rep(ncolors,w))
    } else {
      if (length(bottom) != w){
        stop('incorrect size for argument "bottom". Length of "bottom" should be "w" = ', w)
      }
    }
    if (is.null(right)){
      right <- as.vector(rep(ncolors,h))
    } else {
      if (length(right) != h){
        stop('incorrect size for argument "right". Length of "right" should be "h" = ', h)
      }
    }
    if (is.null(corner)){
      corner <- as.vector(rep(ncolors,4))
    } else {
      if (length(corner) != 4){
        stop('incorrect size for argument "corner". Length of "corner" should be ', 4)
      }
    }
    
    myBorder <- new(Border, h, w, nei, param)
    myBorder$setBorders(top, left, bottom, right, corner)
    myBlock$correctFactor(myBorder)
    
    Z <- myBlock$recursion.cond.mem(myBorder)
    myBlock$sampleCond(Z, myBorder)
    
  }
  
  img <- myBlock$vertices()
  img <- matrix(img, h, w)
  
  if (view){
    par(mar=rep(1,4))
    image(img, col = grey(0:255/255), axes = FALSE)
    box(col = "grey50", lwd = 2)
  }
  
  return(img)
}


sampler.mrf <- function(iter, sampler = "Gibbs", h, w, param, ncolors = 2, nei = 4, 
                        pot = NULL, top = NULL, left = NULL, 
                        bottom = NULL, right = NULL, corner = NULL, 
                        initialise = TRUE, random = TRUE, view = FALSE){
  
  # checking parameters
  if(missing(iter)){
    stop('\n argument "iter" is missing with no default')
  }
  if(!sampler %in% c("Gibbs", "SW")){
    stop('\n argument "sampler" is not valid. The latter must be one of "Gibbs" or "SW"')
  }
  if(missing(h)){
    stop('\n argument "h" is missing with no default')
  }
  if(missing(w)){
    stop('\n argument "w" is missing with no default')
  }
  if(missing(param)){
    stop('argument "param" is missing with no default')
  } else {
    if (length(param) != 1){
      stop('incorrect argument value. "param" has to be scalar')
    }
  }
  if (nei != 4 & nei != 8){
    stop('incorrect argument values. Allowed values for "nei" are 4 or 8')
  }
  
  if (nei == 4){
    param <- param*c(1., 1., 0., 0.)
  } else {
    param <- param*c(1., 1., 1., 1.)
  }
  
  # Init. the graph and set the borders
  if (is.null(pot)){
    myLattice <- new(Lattice, h, w, ncolors, nei, param)
  } else {
    if (length(pot) != ncolors){
      stop('incorrect size for argument "pot". Length of "pot" should be "ncolors" = ', ncolors)
    }
    myLattice <- new(Lattice, h, w, ncolors, nei, param, pot)
  }
  
  
  if(is.null(top) & is.null(left) & is.null(bottom) & is.null(right) & is.null(corner)){
    if(sampler == "Gibbs"){
      myLattice$GibbsSampler(iter, random, initialise)
    } else {
      myLattice$SWSampler(iter, initialise)
    }
  } else {
    # setting artificial borders if needed
    if (is.null(top)){
      top <- as.vector(rep(ncolors,w))
    } else {
      if (length(top) != w){
        stop('incorrect size for argument "top". Length of "top" should be "w" = ',w)
      }
    }
    if (is.null(left)){
      left <- as.vector(rep(ncolors,h))
    } else {
      if (length(left) != h){
        stop('incorrect size for argument "left". Length of "left" should be "h" = ', h)
      }
    }
    if (missing(bottom)){
      bottom <- as.vector(rep(ncolors,w))
    } else {
      if (length(bottom) != w){
        stop('incorrect size for argument "bottom". Length of "bottom" should be "w" = ', w)
      }
    }
    if (is.null(right)){
      right <- as.vector(rep(ncolors,h))
    } else {
      if (length(right) != h){
        stop('incorrect size for argument "right". Length of "right" should be "h" = ', h)
      }
    }
    if (is.null(corner)){
      corner <- as.vector(rep(ncolors,4))
    } else {
      if (length(corner) != 4){
        stop('incorrect size for argument "corner". Length of "corner" should be ', 4)
      }
    }
    
    myBorder <- new(Border, h, w, nei, param)
    myBorder$setBorders(top, left, bottom, right, corner)
    
    if(sampler == "Gibbs"){
      myLattice$GibbsSamplerCond(iter, myBorder, random, initialise)
    } else {
      myLattice$SWSamplerCond(iter, myBorder, initialise)
    }
    
  }
  
  img <- myLattice$vertices()
  img <- matrix(img, h, w)
  
  if (view){
    par(mar=rep(1,4))
    image(img, col = grey(0:255/255), axes = FALSE)
    box(col = "grey50", lwd = 2)
  }
  
  return(img)
}