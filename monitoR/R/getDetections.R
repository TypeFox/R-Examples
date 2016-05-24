# For extracting the detections from a detectionList object
# Modified: 2014 MAR 29

getDetections <-
function(
   detection.obj,            # A 
   which.one=names(detection.obj@detections),     # Name of template(s) that should be returned
   id=NULL,                  # Additional column that should be added to first position in returned data frame
   output='data frame'       # Format of output, 'data frame' or 'list'
) {

   # Pull out the template names for adding to output
   hits <- detection.obj@detections

   # Make sure which.one has template names
   if(is.numeric(which.one)) which.one <- names(detection.obj@detections)[which.one]

   # Empty list for holding hits
   results <- list()
   for(i in which.one) {
      dat <- hits[[i]]
      if(nrow(dat)>0) {
         dat$id <- id
         dat$template <- i
      } else {
         if(!is.null(id)) dat$id <- do.call(class(id), list())
         dat$template <- character()
      }
      # Change order of columns
      if(is.null(id)) 
         dat <- dat[, c(ncol(dat), 1:(ncol(dat)-1))] else 
         dat <- dat[, c(ncol(dat)-1:0, 1:(ncol(dat)-2))] 
         results[[i]] <- dat
   }
   
   # Collapse list
   if(grepl('data', output)) 
      results <- rbindf(results) else if(!grepl('list', output)) 
      stop('Output option not recognized')

   return(results)
}
