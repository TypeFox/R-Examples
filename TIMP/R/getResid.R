"getResid" <-
function (data, modspec = list(), datasetind = vector(), modeldiffs = list(), 
opt = opt()) 
{
    optN <- opt
    result <- fitModel(data, modspec, datasetind, modeldiffs, optN)
    resultlist <- result$currModel@fit@resultlist 
    m <- result$currModel@modellist 
    svdresidlist <- list()
    for (i in 1:length(m)) {
            residuals <- matrix(nrow = m[[i]]@nt, ncol = m[[i]]@nl)
            for (j in 1:length(resultlist[[i]]@resid)) {
                residuals[, j] <- resultlist[[i]]@resid[[j]]
            }
            svdresidlist[[length(svdresidlist) + 1]] <- doSVD(residuals, 5, 5)
	    svdresidlist[[length(svdresidlist)]]$weight <- m[[i]]@weight
	    svdresidlist[[length(svdresidlist)]]$weightM <- m[[i]]@weightM
   }
   svdresidlist
}
