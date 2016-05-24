extractBEASTstats <- function(file) {
  
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  
    # isolate NEWICK string
    # ------------------------------
    X <- X[grep("tree TREE1[[:space:]]+=", X)]
    X <- gsub("tree TREE1[[:space:]]+= \\[&R\\] ", "", X)
    
    # store stats per node in a list
    # ------------------------------
    tab <- unlist(strsplit(X, "\\["))[-1]
    tab <- gsub("&|;|\\]", "", tab)
    tab <- gsub(":.+$", "", tab)
    foo <- function(x){x <- unlist(strsplit(x, ",")); x}
    tab <- lapply(tab, foo)
    
    # tidy up this list
    # -----------------
    for (i in seq(along = tab)){
    	ind <- grep("[{]", tab[[i]])
    	names <- gsub("=.+$", "", tab[[i]][ind])
    	tab[[i]][ind] <- gsub("[{]", "", tab[[i]][ind])
    	tab[[i]][ind] <- gsub("=", "_MIN=", tab[[i]][ind])
    	tab[[i]][ind + 1] <- gsub("[}]", "", tab[[i]][ind + 1])
    	tab[[i]][ind + 1] <- paste(paste(names, "MAX=", 			sep = "_"), tab[[i]][ind + 1])
    }
    
    ttab <- data.frame()
    stats <- unique(gsub("=.+$", "", unlist(tab)))
    for (i in seq(along = tab)){
    	for (j in seq(along = stats)){
    		ind <- grep(paste("^", stats[j], "=", sep = ""), 				tab[[i]])
    		if (length(ind) > 0){
    			v <- as.numeric(gsub(paste(stats[j], "=", 					sep = ""), "", tab[[i]][ind]))
    			ttab[i, j] <- v
    		}
    	}
    }
    colnames(ttab) <- stats 
    tip <- which(is.na(ttab$posterior))    
    ttab   
}