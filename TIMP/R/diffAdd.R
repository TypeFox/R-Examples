"diffAdd" <-
function (modellist, diffsadd) 
{
	## diffsadd has structure 
	## list(list(what, dataset, ind, start), ...)
	## e.g.,  list(list(what="kinpar", datatset=2, 
	##             start=28.4, ind=c(1,2)), ...)
	## dataset can be a vector of indices into dataset list
		
	for(diffs in diffsadd){
	  for(i in 1:length(diffs$dataset)) {
	      if(length(diffs$ind) == 1)		
		 slot(modellist[[diffs$dataset[i]]], 
	            diffs$what) <- 
	            append(slot(modellist[[diffs$dataset[i] ]], 
	            diffs$what), diffs$start, after = diffs$ind-1)  
	      if(length(diffs$ind) == 2)	
		 slot(modellist[[diffs$dataset[i] ]], 
	            diffs$what)[[diffs$ind[1]]] <- 
	            append(slot(modellist[[diffs$dataset[i] ]], 
	            diffs$what)[[diffs$ind[1]]], diffs$start, 
		    after = diffs$ind[2]-1)  
	  }
	
	}
	modellist
}

