"cat.to.list" <-

function(x, a)



{



	a <- as.character(a)



	label <- unique(a)



	out <- as.list(1.:length(label))



	names(out) <- label



	for(k in 1.:length(label)) {



		out[[k]] <- x[label[k] == a]



	}



	out



}

