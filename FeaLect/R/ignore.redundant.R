ignore.redundant <-
function(F, num.of.values=1){    
### rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	## If the value a feature is the same for all points (e.g. =0), we can ignore that feature. 
	### INPUT
	# F:				contains the features matrix.
	# num.of.values:	(=1 by default) 
	#					If the number of distint values that a feature takes is less than this number,
	# 					that features should be ignored. 
	#___________________________________________________________________________________________________
	
	omited <- F
	colnames(omited) <- 1:dim(F)[2]	# F might contain columns with repeating names.
    features.var <- c()
    for (ith.feature in colnames(omited)){ # Examining all features,
   	   features.var[ith.feature] <- var(omited[,ith.feature])	# not used for now.
   	    if (length(unique(omited[,ith.feature]))<=num.of.values){ # Ignore this feature, it has just one value.    
			# To find out the "index" which is a number; ith.feature is not good because it's a name.
			#message(ith.feature)
   	       	feature.index <- which(colnames(omited) == ith.feature)	
   	       omited <- omited[,-feature.index]
   	    }#End if.        
   	}#End for.
    #X11(); plot(features.var, main="features variance")
    remained.columns <- as.numeric(colnames(omited))
    result <- as.matrix(F[ ,remained.columns])
    return(result)
#rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr    
}#End ignore.redundant <- function(F).

