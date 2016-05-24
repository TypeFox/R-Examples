"diffRemove" <-
function (modellist, diffsremove) 
{
	## diffsremove has structure 
	## list(list(what, ind, dataset), ...)
	## e.g.,  list(list(what="kinpar", ind=c(1,2), dataset=2), ...)
	## dataset can be a vector of indices into dataset list
		
	for(diffs in diffsremove){
	     for(i in 1:length(diffs$dataset)) {
	      slot(modellist[[diffs$dataset[i]]], 
	             diffs$what) <- 
	      slot(modellist[[diffs$dataset[i]]], 
	             diffs$what)[[diffs$ind[1]]][- diffs$ind[2]]  
	      
	     }
	}
	modellist
}

