	#Now construct variance matrices
	getVMat.onePhase =
	function(Z.Phase1, design.df, var.comp = NA){

	#Now construct variance matrices
		if(all(is.na(var.comp))){
		  lapply(Z.Phase1, function(x) x %*% t(x))
		}else{
			lapply(makeBlkDesMatrix(design.df, var.comp), function(x) x %*% t(x))
		}

  }

