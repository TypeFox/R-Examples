#' Smooth Using Gaussian Window
#' 
#' @description The specific function for smoothing using the gaussian window function
#' @param x numeric vector of values to smooth, error will be thrown if not provided.
#' @param window the length of the smoothing window, if an integer, represents
#' number of items, else, if a value between \code{0} and \code{1}, represents the 
#' proportion of the input vector
#' @param alpha parameter to determine the breadth of the gaussian window, yielding more or less sensitive 
#' smoothing characteristics
#' @param ... not used
#' @param tails Logical value as to whether the tail regions should be included or not.
#' @name   smth.gaussian
#' @rdname smth.gaussian
#' @examples
#'   y  = runif(100)
#'   ys = smth.gaussian(y)
#' @export
smth.gaussian <- function( x       = stop("Numeric Vector 'x' is Required"),
                           window  = getOption('smoother.window'),
                           alpha   = getOption('smoother.gaussianwindow.alpha'),
                           ...,
                           tails   = getOption('smoother.tails')){
  
  #Check Numeric Arguments
  if(!is.numeric(x) | !is.numeric(alpha)){stop("argument 'x' and 'alpha' must be numeric",call.=FALSE)} 
  
  #Determine the Window Length
  windowLength   = .determineWindowLength(x,window)
  
  #Hidden Gaussian Window Function
  makeWindow  = function(w,a){
    hw  = abs(w/2.0) #halfwidth
    e   = exp(1)     #eulers number
    a   = abs(a)     #alpha
    ret = sapply(c(0:(w - 1)),function(x){
            n = x - as.integer(hw)
            k = -0.5*(a*n/hw)^2
            e^k
          })
    ret
  }
  
  #Continue with the Convolustion
  w = makeWindow(windowLength,alpha[1])
  
  #Length of F & G Vectors
  sizeW = length(w)
  sizeD = length(x)
  
  #Normalize the window
  w = .normalize(w)
  
  #Prepare 
  hkwL = as.integer(sizeW/2) #Half Kernel Width, for left
  hkwR = sizeW - hkwL        #Half Kernel Width, rhs balance to total width
  
  #Now Do the Gaussian Smoothing.
  ret  = sapply(c(1:sizeD),function(i){
    ix.d = c((i-hkwL):(i+hkwR-1))   #The Ideal Range
    ix.w = which(ix.d %in% 1:sizeD) #Indexes of Values, Window
    ix.d = ix.d[ix.w]               #Indexes of Values, Data
    
    #Determine Normalized Final Window and Data
    W.nm = .ifthenelse(length(ix.w) != sizeW,.normalize(w[ix.w]),w) #Re-Normalize if Needed
    D.nm = x[ix.d] 
    
    #Dot Product
    as.numeric(D.nm %*% W.nm)
  })
  
  #Remove Tail Regions
  if(!tails){ret[c(1:hkwL,(sizeD - hkwR + 1):sizeD)] = NA}
  
  #Done
  ret
}

