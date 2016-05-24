#' Print results of observer detection tables
#'
#' Provides formatted output for detection tables
#'
#' @aliases print.det.tables
#' @export
#' @param x result of call to ddf
#' @param \dots unused unspecified arguments for generic print
#' @return None
#' @author Jeff Laake
#' @seealso \code{\link{plot.det.tables}}
#' @keywords utility
print.det.tables<-function(x, ...){
  if(!is.null(x$Observer1)){
    cat("\nObserver 1 detections\n")
    print(x$Observer1)
  }
  if(!is.null(x$Observer2)){
   cat("\nObserver 2 detections\n")
   print(x$Observer2)
  }
  if(!is.null(x$Duplicates)){
    cat("\nDuplicate detections\n")
    print(x$Duplicates)
  }
  if(!is.null(x$Pooled)){
    cat("\nPooled detections\n")
    print(x$Pooled)
  }
  if(!is.null(x$Obs1_2)){
    cat("\nObserver 1 detections of those seen by Observer 2\n")
    print(x$Obs1_2)
  }
  if(!is.null(x$Obs2_1)){
    cat("\nObserver 2 detections of those seen by Observer 1\n")
    print(x$Obs2_1)
  }

  invisible()
}
