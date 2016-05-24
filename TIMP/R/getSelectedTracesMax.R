"getSelectedTracesMax" <- function (multimodel, t, plotoptions) 
{
  sp <-getSpecList(multimodel, t)[[1]]
  ## this assumes the spectra are same for all datasets 
  
  sorted <- apply(sp, 2, sort, index.return = TRUE, decreasing = TRUE)
  ## now each element of sorted is a sorted component
  
  cnt <- plotoptions@nummaxtraces 
  cc <- ncol(sp)
  
  ll <- 0
  selectedtraces <- vector()
  while(length(selectedtraces) < cnt && ll < nrow(sp)) {
    for(i in 1:cc) {
      ss <- sorted[[i]]$ix
      if(ll==0)
        selectedtraces <- ss[1]
      else {
        addto <- which(! ss %in% selectedtraces)
        if(length(addto) > 0) 
          selectedtraces <- append(selectedtraces, ss[addto[1]])
      }
      
      ll <- ll + 1 
    }
  }
  cat("The following traces will be plotted:\n", toString(selectedtraces),"\n")
  selectedtraces
}
