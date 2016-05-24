niche.overlap.boot.pair <-
function (vectorA, vectorB, method = c("levins","schoener","petraitis","pianka","czech","morisita"), 
          times = 999, quant = c(0.025, 0.975)) 
{
    method <- match.arg(method)
    if(!length(vectorA)==length(vectorB)){
	    stop("Length of vectorA differs from lengths of vectorB")
	}
    booted <- rep(NA, times)
    obs <- niche.overlap.pair(vectorA, vectorB, method = method)
    
    for (i in 1:times){
        ind <- sample(1:length(vectorA), size = length(vectorA), replace = TRUE)
        booted[i] <- niche.overlap.pair(vectorA[ind], vectorB[ind], method = method)
    }
    
    result <- c(obs, mean(booted), sd(booted), quantile(booted, quant, na.rm = TRUE), times)
	names(result) <- c("Observed","Boot mean","Boot std","Boot CI1", "Boot CI2", "times")
	return(round(result,3))
}

