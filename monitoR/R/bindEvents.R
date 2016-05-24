# Modified: 2015 JULY 21

bindEvents <-
function(
   rec, 
   file, 
   by.species=TRUE, 
   parallel=FALSE, 
   return.times=FALSE
) {

   # Read in the event file
   if(class(file) != "data.frame") events <- read.csv(file=file)
   else events <- file
   events['duration'] <- events[, 'end.time']-events[, 'start.time']
   if(by.species) {
      spp <- events$name
      events <- split(events, events$name)
      events <- lapply(events, function(x) data.frame(x, start.time.collapsed=cumsum(x$duration)-x$duration, end.time.collapsed=cumsum(x$duration)))
   } else {
      events$start.time.collapsed=cumsum(events$duration)-events$duration
      events$end.time.collapsed=cumsum(events$duration)
   }

   # Collapse recording
   # CAN SIMPLIFY THIS A BIT BY USING A SINGLE LAPPLYFUN AS IN SCCCOR ETC
   if(parallel && by.species) {
      collapsed <- parallel::mclapply(X=events, FUN=function(x) collapseClips(rec=rec, start.times=x$start.time, end.times=x$end.time), max(1, parallel::detectCores()-1))
      #events <- unsplit(events, spp)
   } else if(by.species) {
      collapsed <- lapply(X=events, FUN=function(x) collapseClips(rec=rec, start.times=x$start.time, end.times=x$end.time))
      #events <- unsplit(events, spp)
   } else {
      collapsed <- list(collapseClips(rec=rec, start.times=events$start.time, end.times=events$end.time))
   }

   if(return.times) return(list(times=events, wave=collapsed)) else return(collapsed)
}
