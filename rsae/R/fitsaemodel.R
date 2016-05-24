fitsaemodel <-
function(method, model, ...){
   thecall <- match.call()
   # check if model is appropriate
   if (inherits(model, "saemodel")){
      # check if method exists
      m <- match(method, c("ml", "huberm", "tukeys"))
      if (is.na(m)) stop(paste("Method: ", method, " is not supported! \n", sep=""))
      # ml and huberm
      if (m == 1 | m == 2){
	 tmp <- .fitsaemodel.huberm(method, model, ...)
      }
      # tukeys 
      if (m == 3) stop("Tukey S-estimator not implemented, yet!\n")
   }else{
      stop(paste("Model", thecall[3] ,"must be an instance of 'saemodel' class.\n"))
   }
   tmp
}

