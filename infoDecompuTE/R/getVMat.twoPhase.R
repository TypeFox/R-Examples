	#Now construct variance matrices
	getVMat.twoPhase =
	function(Z.Phase1, Z.Phase2, design.df, var.comp = NA){

		if(all(is.na(var.comp))){
			lapply(c(Z.Phase1, Z.Phase2[-1]), function(x) x %*% t(x))
		}else{
			 lapply(makeBlkDesMatrix(design.df, var.comp), function(x) x %*% t(x))
		}

  }

