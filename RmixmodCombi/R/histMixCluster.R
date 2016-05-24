histMixCluster <-
function (mixmodCombiOutput, nbCluster = mixmodCombiOutput@mixmodOutput@bestResult@nbCluster, data = mixmodCombiOutput@mixmodOutput@data, variables = colnames(data), permutIndices = 1:mixmodCombiOutput@mixmodOutput@bestResult@nbCluster, combiM = NULL, xlab = rep("", length(variables)), main = paste("Histogram of", variables), 
    ...) 
{
    if (!is(mixmodCombiOutput, "MixmodCombi")) 
        stop("'x' must be a MixmodCombi object!")
    if (!is.matrix(data) & !is.data.frame(data) & !is.vector(data)) 
        stop("'data' must be a vector, a data.frame or a matrix object!")
    if ((length(variables) == 0) & (ncol(data) > 1)) 
        stop("'variables' is empty!")
    if (length(variables) > ncol(data)) 
        stop("List of variables too long!")
    if (sum(!(variables %in% colnames(data)))) 
        stop("At least one variable is unknown!")
	 if (! (is.vector(permutIndices) & length(permutIndices) == mixmodCombiOutput@mixmodOutput@bestResult@nbCluster & all((1:mixmodCombiOutput@mixmodOutput@bestResult@nbCluster) %in% permutIndices)))  
		  stop("permutIndices must be a permutation vector of 1:nbCluster")
	 if (is.null(combiM))
		{
			combiM = diag(nrow = mixmodCombiOutput@mixmodOutput@bestResult@nbCluster, ncol = mixmodCombiOutput@mixmodOutput@bestResult@nbCluster)
			for (K in mixmodCombiOutput@mixmodOutput@bestResult@nbCluster : nbCluster)
            {
            	combiM <- mixmodCombiOutput@hierarchy[[K]]@combiM %*% combiM
            }
        }
    if (! (is.matrix(combiM) & all(dim(combiM) == c(nbCluster, mixmodCombiOutput@mixmodOutput@bestResult@nbCluster)) & all(combiM %in% c(0,1)) & sum(combiM) == mixmodCombiOutput@mixmodOutput@bestResult@nbCluster & all(colSums(combiM) == rep(1, mixmodCombiOutput@mixmodOutput@bestResult@nbCluster))) )
    	  stop("combiM must be a permutation matrix of dimension c(nbCluster, mixmodCombiOutput@mixmodOutput@bestResult@nbCluster)")
    
    op <- par(no.readonly = TRUE)
    if (ncol(data) == 1) {
        indices <- 1
    }
    else {
        indices <- which(colnames(data) %in% variables)
    }
    nvar <- length(indices)
    if (is(mixmodCombiOutput@mixmodOutput@bestResult@parameters, "GaussianParameter")) {
        if (isQualitative(data)) 
            stop("data must contain only quantitative variables!")
        if (nvar < 4 & nvar > 1) 
            par(mfrow = c(1, nvar))
        else if (nvar >= 4) {
            nrow <- round(sqrt(nvar))
            if (is.wholenumber(sqrt(nvar))) 
                ncol <- sqrt(nvar)
            else ncol <- sqrt(nvar) + 1
            par(mfrow = c(nrow, ncol))
        }
        i <- 1
        for (j in indices) {
            xaxis <- seq(min(data[, j]), max(data[, j]), by = 1e-04)
            density <- matrix(nrow = mixmodCombiOutput@mixmodOutput@bestResult@nbCluster, ncol = length(xaxis))
            for (k in 1:mixmodCombiOutput@mixmodOutput@bestResult@nbCluster) 
            {
                density[k, ] <- mixmodCombiOutput@mixmodOutput@bestResult@parameters["proportions", k] * dnorm(xaxis, mixmodCombiOutput@mixmodOutput@bestResult@parameters["mean", k][j], sqrt(mixmodCombiOutput@mixmodOutput@bestResult@parameters["variance", k][j, j]))
            }
            density <- density[permutIndices, ]
            mixDensity <- combiM %*% density
            mixture <- apply(mixDensity, 2, sum)
            h <- hist(data[, j], xlab = xlab[i], main = main[i], ...)
            ratio <- max(h$counts)/max(mixture)
            mixDensity <- mixDensity * ratio
            mixture <- mixture * ratio
            lines(xaxis, mixture, col = "azure4", lty = 1, lwd = 4)
            for (k in 1:nbCluster) 
            {
                lines(xaxis, mixDensity[k,], col = k + 1, lty = 2, lwd = 2)
            }
            i <- i + 1
        }
    }
    else if (is(mixmodCombiOutput@mixmodOutput@bestResult@parameters, "MultinomialParameter")) {
        stop("x must contain Gaussian parameters. See barplot() to plot multinomial parameters.")
    }
    else {
        stop("Uknown type of parameters!")
    }
    par(op)
}
