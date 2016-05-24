niche.overlap.boot <-
function(mat, method = c("pianka", "schoener", "petraitis",  "czech", "morisita", "levins"), times = 999, quant = c(0.025, 0.975)){
     method = match.arg(method)
     result <- data.frame(col1 = rep(NA, 6))
     for(i in 1:(ncol(mat)-1)){
        for(j in (i+1):(ncol(mat))){
            booted <- niche.overlap.boot.pair(mat[,i], mat[,j],method = method, times = times, quant = quant)
            booted1 <- data.frame(booted)
            colnames(booted1) <- paste(colnames(mat)[i], "-", colnames(mat)[j], sep = "")
            result <- cbind(result, booted1)	 
     	}
     }
     result <- result[,-1]
     return(t(result))
}

