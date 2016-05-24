#' @export
#' @title Rolling STA/LTA
#' @param x an \R numeric vector
#' @param n_sta integer STA window size
#' @param n_lta integer LTA window size
#' @param increment integer shift to use when sliding the window to the next location
#' @description Fast rolling STA/LTA using C++/Rcpp.  Additional performance gains can be achieved by 
#' skipping \code{increment} values between calculations.
#' 
#' The STA/LTA ratio method is used for automatic detection of seismic signal arrival times.
#' @details 
#' The \code{roll_stalta} function described here does no preprocessing of the incoming 
#' data and merely calculates the ratio of the average value in the STA window to the average value 
#' in the LTA window. Windows are aligned so that the index is at the left edge of the STA window and 
#' at the right edge of the LTA window.
#'
#' \deqn{ STA(x_i) = \frac{1}{ns}\sum_{j=i}^{i+ns}{x_i} }
#' 
#' \deqn{ LTA(x_i) = \frac{1}{nl}\sum_{j=i-nl}^{i}{x_i} }
#' 
#' \deqn{ r_i = \frac{STA_i}{LTA_i} } 
#' 
#' %% TODO:  use \preformatted or equivalent whent it becomes available in roxygen2
#' \code{[---------- LTA --------*]........}
#' 
#' \code{.......................[*- STA --]}
#'
#' For proper use of this algorithm seismic data should be preprocessed as in the example below with:
#'  \itemize{
#'    \item{demean, detrend and taper the raw signal}
#'    \item{square the processed signal to get power}
#'  }
#'
#' With \code{increment=1}, this function is equivalent to, eg:
#'
#' \code{sta <- roll_mean(x,3,1,"left")}
#' 
#' \code{lta <- roll_mean(x,30,1,"right")}
#' 
#' \code{r <- sta/lta}
#'
#' For increments greater than one, the rolling means above will not align properly,
#' hence the need for a dedicated \code{roll_stalta} function.
#' 
#' Values within \code{n_lta-1} of the beginning or
#' \code{n_sta-1} of the end of \code{x} are set to \code{NA}.
#'
#' Setting \code{increment} to a value greater than one will result in \code{NA}s for
#' all skipped-over indices.
#' @return A vector of values of the same length as \code{x} with each point containing the STA/LTA
#' ratio at that point.
#' @examples
#' # Contrived example
#' x <- rep(c(1,5,3,2,1),each=20)
#' p <- roll_stalta(x,3,6)
#' plot(x, pch=17, cex=0.8, ylim=c(0,max(x)),
#'      main="Test of roll_stalta on artificial data")
#' points(p,cex=1.5,col='red',type='b')
#' legend('topleft',
#'        legend=c('data','STA/LTA'),
#'        pch=c(17,1),
#'        col=c('black','red'))
#'
#' # Real example requiring the 'seismic' package
#' \dontrun{
#' require(seismic)
#'  
#' # Create a new IrisClient
#' iris <- new("IrisClient")
#'   
#' # Seismic data with a large quake
#' starttime <- as.POSIXct("2010-02-27 06:30:00", tz="GMT")
#' endtime <- as.POSIXct("2010-02-27 07:00:00", tz="GMT")
#' st <- getDataselect(iris,"IU","ANMO","00","BHZ",starttime,endtime)
#' tr <- st@@traces[[1]]
#'  
#' # Preprocess the data
#' x <- DDT(tr)@@data
#'  
#' # Calculate the first break 'picker'
#' n_sta <- 3 * tr@@stats@@sampling_rate
#' n_lta <- 10 * n_sta
#' p <- roll_stalta(x^2,n_sta,n_lta)
#'  
#' first_break <- which(p == max(p,na.rm=TRUE))
#'
#' plot(x,type='l',
#'      main='Test of STA/LTA first break picker on raw seismic data')
#' abline(v=first_break,col='red')  
#'}
#' @references \href{http://en.wikipedia.org/wiki/First_break_picking}{First Break Picking}


roll_stalta <- function( x, n_sta, n_lta, increment=1 ) {
  
  if ( !is.vector(x) ) {
    
    stop("the x supplied is not a vector")
    
  } else {
    
    if ( n_sta > length(x) ) {
      stop("n_sta cannot be greater than length(x).")
    } else if ( n_lta > length(x) ) {
      stop("n_lta cannot be greater than length(x).")
    }
    
    # Avoid infinite loop
    if ( increment < 1 ) {
      stop("increment must be >= 1.")
    }
    
    result <- .Call( "seismicRoll_roll_stalta_numeric_vector", 
                     x, as.integer(n_sta), as.integer(n_lta), as.integer(increment),
                     PACKAGE="seismicRoll")
    
    return (as.numeric(result))
    
  }
  
}

