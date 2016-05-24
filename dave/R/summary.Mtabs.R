summary.Mtabs <-
function(object,...,range=NULL) {
     o.Mt<- object
     test<- is.null(range) 
     if(test == TRUE) range<- c(1,ncol(o.Mt$centroids))
     print(t(o.Mt$centroids[,range[1]:range[2]]),quote=FALSE)
     }
