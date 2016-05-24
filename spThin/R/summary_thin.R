#' @export summaryThin
#' @title Summary method for results of thin function
#' 
#' @description
#' Summarize the results of \code{thin} function.
#' 
#' @param thinned A list of data.frames returned by \code{\link{thin}} function.
#' @param show logical; if \code{TRUE},the summary values are printed at the console.
#'  
#' @return Returns a list with the (1) maximun number of records, (2) number of data frames
#' with maximun number of records and (3) a table with the number of data frames per 
#' number of records.
#' 
#' @seealso \code{\link{thin.algorithm}}
#' @seealso \code{\link{thin}}


summaryThin <- function(thinned, show=TRUE){

  ## Repetition number
  reps <- length(thinned)
  
  ## Look at the number of locs kept in each thinned dataset
  ## by determining the number of rows in each returned data.frame
  lat.long.thin.count <- unlist(lapply(thinned, nrow ))
  max.lat.long.thin.count <- max(lat.long.thin.count)  
  
  ## Number of data.frames with max records
  n.max.data.frame <- sum(lat.long.thin.count==max.lat.long.thin.count)
  
  n.data.frame.records <- table(lat.long.thin.count)
  n.records <- as.numeric(names(n.data.frame.records))
  Frequency <- as.numeric(n.data.frame.records)
  table2 <- as.data.frame(rbind(n.records, Frequency))
  colnames(table2) <- rep(" ", ncol(table2))
  rownames(table2)[1] <- "Number of records"
 
  if(show){ 
    cat(paste("Maximum number of records after thinning:", 
              max.lat.long.thin.count))  
    cat(paste("\nNumber of data.frames with max records:", 
              n.max.data.frame))
    cat("\n-----------------------------------------------")
    print(table2)
  }
  
  invisible( list("Maximun_number_of_records" = max.lat.long.thin.count, 
                  "Number_of_data_frames_with_maximun_number_of_records" = n.max.data.frame, 
                  "Number_of_Data_frames_per_number_of_records"=table2))
}
