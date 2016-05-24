msc.nn.kd <- function (y, x, knn = ncol(x)*3, pLevelP = 0.2, pLevel, nLevels, 
                       bw, type = 1, smooth = FALSE, eps=0.01) 
{
    obj <- msc.nn(y,x,knn, pLevelP, pLevel, nLevels, type = type, smooth =
        smooth, eps = eps)
    obj$bw = bw
    class(obj) <- paste(class(obj), ".kd", sep="")
    obj
}

