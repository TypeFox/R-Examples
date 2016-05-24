"diffChange" <-
function (modellist, diffsadd) 
{
	## diffschange has structure 
	## list(list(what, dataset, spec), ...)
	## e.g.,  list(list(what="kinpar", datatset=2, 
	##             spec=c(1,2,3), type="multifree"), ...)
	## dataset can be a vector of indices into dataset list
	## type is optional; if multifree then parameters in 
	## what are free per-dataset 
	
	for(diffs in diffsadd) {
		if(length(diffs$type) == 0)
		     diffs$type <- "perd"		 
		for(i in  1:length(diffs$dataset)) 
		   slot(modellist[[diffs$dataset[i] ]], 
		   diffs$what) <- diffs$spec
	}
	newl <- list()
	for(i in 1:length(diffsadd)) {
	      if(length(diffsadd[[i]]$type) == 0)
		 diffsadd[[i]]$type <- "perdataset"
	      if(diffsadd[[i]]$type == "multifree") {
		 for(j in 1:length(diffsadd[[i]]$dataset)){ 
		       newl[[length(newl)+1]] <- diffsadd[[i]]
		       newl[[length(newl)]]$dataset <- diffsadd[[i]]$dataset[j]
		 }
	      }
	      else	
			newl[[length(newl)+1]] <- diffsadd[[i]]
        }
        
			       
	list(modellist=modellist, diffsadd = newl)

}
