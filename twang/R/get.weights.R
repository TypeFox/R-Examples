get.weights <- function(ps1, stop.method = NULL, estimand = NULL, withSampW = TRUE)
{
   if(class(ps1)=="ps"){
   if(is.null(estimand)) estimand <- ps1$estimand
   
   if(!(estimand %in% c("ATT","ATE"))) stop("estimand must be either \"ATT\" or \"ATE\".")
   if(estimand != ps1$estimand){
   	warning("Estimand specified for get.weights() differs from the estimand used to fit the ps object.")
   }
   if(length(stop.method)>1) stop("More than one stop.method was selected.")
   if(!is.null(stop.method))
   {
   	stop.method.long <- paste(stop.method, ps1$estimand, sep=".")
      i <- match(stop.method.long, names(ps1$w))
      if(is.na(i)) stop("Weights for stop.method=",stop.method, " and estimand=", estimand, " are not available. Please a stop.method and used when fitting the ps object.") 
#      Available options: ",names(ps1$ps),".")
   } else
   {
   	warning("No stop.method specified.  Using ", names(ps1$ps)[1], "\n")
      i <- 1
   }
  
   if (estimand == "ATT") {
   	w <- with(ps1,  treat + (1- treat) * ps[[i]]/(1-ps[[i]]))
   	if(withSampW) w <- w * ps1$sampw
   	return(w)
   }
   else if (estimand == "ATE")
   { 
      w <- with(ps1, treat/ps[[i]] + (1-treat)/(1-ps[[i]]))
      if(withSampW) w <- w* ps1$sampw
      return(w)
   }
   }
   if(class(ps1) == "mnps"){
   if(is.null(estimand)) estimand <- ps1$estimand
   
   if(!(estimand %in% c("ATT","ATE"))) stop("estimand must be either \"ATT\" or \"ATE\".")
   if(estimand != ps1$estimand){
   	warning("Estimand specified for get.weights() differs from the estimand used to fit the ps object.")
   }
   if(length(stop.method)>1) stop("More than one stop.method was selected.")
   if(!is.null(stop.method))
   {
   	stop.method.long <- paste(stop.method, ps1$estimand, sep=".")
      i <- match(stop.method.long, names(ps1$psList[[1]]$ps))
      if(is.na(i)) stop("Weights for stop.method=",stop.method, " and estimand=", estimand, " are not available. Please a stop.method and used when fitting the mnps object.") 
#      Available options: ",names(ps1$ps),".")
   } else
   {
   	if(length(ps1$stopMethods) > 1) warning("No stop.method specified.  Using ", names(ps1$psList[[1]]$ps)[1], "\n")
      i <- 1
   }
  
   if (estimand == "ATT") {
   	w <- rep(0, nrow(ps1$data))
   	w[ps1$treatVar == ps1$treatATT] <- 1
   	for(j in 1:length(ps1$levExceptTreatATT)){
   		### in mnps, treatATT is the treatment in the individual fits
   		### fill in weights for categories other than treatATT first, using (1-ps)/ps
   		w[ps1$treatVar %in% c(ps1$levExceptTreatATT[j], ps1$treatATT)] <- with(ps1$psList[[j]],  (1-treat) * ps[[i]]/(1-ps[[i]]))
   		}
   		w[ps1$treatVar == ps1$treatATT] <- 1
   		if(withSampW) w <- w * ps1$sampw
   		return(w)
   }
   else if (estimand == "ATE"){
    	w <- with(ps1$psList[[1]],  treat/ps[[i]])
	   	for(j in 2:length(ps1$psList)){
	   		w <- w + with(ps1$psList[[j]],  treat/ps[[i]])
   		}
   		if(withSampW) w <- w * ps1$sampw
   		return(w)
   		}
  	 }	 
#      w <- with(ps1, treat/ps[[i]] + (1-treat)/(1-ps[[i]]))
#      if(withSampW) w <- w* ps1$sampW
#      return(w)
#   }
   	
   
   if(!(class(ps1) %in% c('ps', 'mnps'))) stop("The object 'ps1' must be of class 'ps' or 'mnps'.")
}
