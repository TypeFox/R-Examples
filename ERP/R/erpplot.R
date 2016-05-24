erpplot <-
function(dta,frames=NULL,ylim=NULL,...) {

   erpdta = as.matrix(dta)
   if (typeof(erpdta)!="double") stop("ERPs should be of type double")
   
   if (!is.null(frames))
      if (length(frames)!=ncol(erpdta))
         stop(paste("frames should be either null or of length",ncol(erpdta)))
   if (!is.null(frames)) {
      if (any(frames!=sort(frames))) stop("frames should be an ascending sequence of integers")
      tframes = frames    
   }
   if (is.null(frames)) tframes = 1:ncol(erpdta)   
   if (is.null(ylim)) ylim = 1.05*diag(apply(apply(erpdta,2,range),1,range))
   matplot(tframes,t(erpdta),type="l",ylim=ylim,...)
}