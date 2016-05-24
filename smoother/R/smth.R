#' Smooth Numerical Data
#' 
#' @description Helper function to smooth numerical data using methods specified by the user. 
#' 
#' @details At this moment in time, the only method is the \code{'gaussian'} window function (similar to the Matlab 
#' Gaussian Window Smoothing Function) and a number of moving averages \code{'sma', 'ema', 'dema'} or \code{'wma'}. 
#' These are functions that allows the user to smooth an input vector, returning vector of the same length as the input. 
#' This can also be achieved using the specific \code{\link{smth.gaussian}} function.
#' @param x numeric vector of values to smooth
#' @param method one of \code{'gaussian', 'sma', 'ema', 'dema'} or \code{'wma'}.
#' @param ... any other arguments to be passed to each specific smoothing methodology.
#' @return a numeric vector of same length as input \code{'x'} vector
#' @references If the \code{'method'} argument is equal to \code{'gaussian'}, then this function is a port of the function 
#' described here: \url{http://goo.gl/HGM47U}, very loosly based of code which has also been ported to c++ here:
#' \url{http://goo.gl/NK79bJ}
#' @seealso Refer to specific man files: \code{\link{smth.gaussian}}, \code{\link[TTR]{SMA}}, \code{\link[TTR]{EMA}}, 
#' \code{\link[TTR]{DEMA}} or \code{\link[TTR]{WMA}}
#' @exportMethod
#' @examples
#' #Prepare Data
#' n  = 1000
#' x  = seq(-pi,pi,length.out=n)
#' y  = sin(x) + (runif(length(x))*0.1) #NOISY DATA
#' ys = smth(y,window = 0.1,method = "gaussian") #SMOOTHING
#' plot(x,y,type="l",col="red")
#' lines(x,ys,col="black",lwd=3)
#' title("Example Smoothing of Noisy Data")
#' @name smth
#' @rdname smth
#' @aliases smth sma ema dema wma gaussian
smth <- function(x       = stop("Numeric Vector 'x' is Required"),
                 method  = getOption('smoother.method'),
                 ...){
  
  #Local Variables
  methods.local = c('gaussian')
  methods.ttr   = c('sma','ema','dema','wma')
  args <- c(list(x),list(...))
  
  #Process
  if(method %in% methods.local){
    return(do.call(paste('smth.',method,sep=""),args))
  } else if(method %in% methods.ttr){
    return(do.call(toupper(method),args)) #IMPORTED FUNCTIONS FROM TTR PACKAGE.
  } else {
    stop(paste("The 'method' argument (",as.character(method),") can only be one of ['",
               paste(c(methods.local,methods.ttr),collapse="','"),"'] at this moment in time.",sep=""),call.=FALSE)
  }
}