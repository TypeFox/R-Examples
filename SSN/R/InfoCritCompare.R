InfoCritCompare <-
function(model.list)
{
   IC <- NULL

   for(i in 1:length(model.list)) {

	 if(class(model.list[[i]])!= "glmssn") {
		 stop("All models must be of type glmssn")
	 }

	 model.name <- NULL
	 ind <- !duplicated(attributes(model.list[[i]]$estimates$theta)$terms)
	 terms<- attributes(model.list[[i]]$estimates$theta)$terms[ind]

         model.name <- paste(terms,collapse=" + ")

	 if(model.list[[i]]$args$family != "gaussian") {
		 model.AIC <- NA
		 model.neg2LogL <- NA
	 }
	 if(model.list[[i]]$args$family =="gaussian"){
		 model.AIC <- AIC(model.list[[i]])
		 model.neg2LogL <- model.list[[i]]$estimates$m2LL
	 }


     IC <- rbind(IC, data.frame(
	   formula = deparse(model.list[[i]]$args$formula, width.cutoff = 500),
       EstMethod = model.list[[i]]$args$EstMeth,
       Variance_Components = model.name,
       neg2LogL = model.neg2LogL,
	   AIC = model.AIC,
       CrossValidationStatsSSN(model.list[[i]]))
       )
   }
   IC
}

