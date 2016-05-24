NC.mrf <- function(h, w, param, ncolors = 2, nei = 4,  
                   pot = NULL, top = NULL, left = NULL, 
                   bottom = NULL, right = NULL, corner = NULL){
  
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
  if(ncolors^h > 2^25){
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
    return(myBlock$recursion())
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
    
    return(myBlock$recursion.cond(myBorder))
  }
}
