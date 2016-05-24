#' Expand a time indexed data set with additional prediction times
#' 

#' 
#' Expands a covariate data frame (or vector) that has a separate time index by
#' inserting prediction times and duplicating the covariate values for all
#' prediction time between subsequent data times.
#' 
#' 
#' @param x Data to be expanded.
#' @param Time Either a character naming the column which contains original
#' time values, or a numeric vector of original times
#' @param predTime prediction times to expand data
#' @param time.col Logical value indicating whether to attach the new times to
#' the expanded data
#' @return data.frame expanded by \code{predTime}
#' @author Devin S. Johnson
#' @examples
#' 
#' #library(crawl)
#' origTime <- c(1:10)
#' x <- cbind(rnorm(10), c(21:30))
#' predTime <- seq(1,10, by=0.25)
#' expandPred(x, Time=origTime, predTime, time.col=TRUE)
#' 
#' @export
"expandPred" <- function(x, Time='Time', predTime, time.col=FALSE)
{
   if(is.character(Time)) {
      Time.name <- Time
      if(!Time.name%in%colnames(x)) stop(paste(Time.name, 'is not in x. PLease specify correct time indicator'))
   }
   else if(is.numeric(Time) & length(Time)==nrow(as.matrix(x))) {
              x <- cbind(Time=Time, x)
              Time.name <- 'Time'
           }
   else stop("Value given for 'Time' is not a recognized format. See crawl documentation")
   predData <- data.frame(predTime)
   colnames(predData) <- Time.name
   newx <- merge(as.data.frame(x), predData,
                 by=c(Time.name), all=TRUE)
   for(i in 1:ncol(newx)) {
      vec <- newx[,i]
      newx[,i] <- vec[!is.na(vec)][cumsum(!is.na(vec))]
   }
   
   if(time.col) return(newx[!duplicated(newx[,Time.name]),])
   else return(newx[!duplicated(newx[,Time.name]),!colnames(newx)%in%Time.name])
}
