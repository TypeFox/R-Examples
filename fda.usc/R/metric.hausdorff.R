metric.hausdorff=function (fdata1, fdata2 = fdata1) 
{
    if (is.fdata(fdata1)) {
        tt <- fdata1[["argvals"]]
        rtt <- fdata1[["rangeval"]]
        nas1 <- apply(fdata1$data, 1, count.na)
        if (any(nas1)) 
            stop("fdata1 contain ", sum(nas1), " curves with some NA value \n")
        else if (!is.fdata(fdata2)) {
            fdata2 <- fdata(fdata2, tt, rtt)
        }
        nas2 <- apply(fdata2$data, 1, count.na)
        if (any(nas2)) 
            stop("fdata2 contain ", sum(nas2), " curves with some NA value \n")
        DATA1 <- fdata1[["data"]]
        DATA2 <- fdata2[["data"]]
        range.t <- rtt
    }
    else {
        if (is.vector(fdata1)) 
            fdata1 <- as.matrix(t(fdata1))
        if (is.vector(fdata2)) 
            fdata2 <- as.matrix(t(fdata2))
        DATA1 <- fdata1
        DATA2 <- fdata2
        range.t <- c(1, ncol(DATA1))
    }

    testfordim <- sum(dim(DATA1) == dim(DATA2)) == 2
    twodatasets <- TRUE
    if (testfordim) 
        twodatasets <- sum(DATA1 == DATA2) != prod(dim(DATA1))
    numgr1 <- nrow(DATA1)
	numgr2 <- nrow(DATA2)
	Mtt=outer(tt,tt,"-")^2
	mdist=array(0,dim=c(numgr1,numgr2))
	etiq1=rownames(DATA1)
	etiq2=rownames(DATA2)
	predi=TRUE
	if (testfordim) {
                if (twodatasets) {
                  predi <- FALSE
                  for (i in 1:numgr1) {
                    ii = i + 1
                    for (ii in i:numgr2) {
					  Mxy=sqrt(outer(DATA1[i,],DATA2[ii,],"-")^2+Mtt)
					  mdist[i, ii] = max(max(apply(Mxy,1,min)),max(apply(Mxy,2,min)))
					  }
                  }
                  mdist = t(mdist) + mdist
                }
            }
            if (predi) {
                for (i in 1:numgr1) {
                  for (ii in 1:numgr2) {
					  Mxy=sqrt(outer(DATA1[i,],DATA2[ii,],"-")^2+Mtt)
					  mdist[i, ii] = max(max(apply(Mxy,1,min)),max(apply(Mxy,2,min)))
                    }
                }
            }
			
	attr(mdist, "call") <- "metric.hausdorff"
	rownames(mdist)<-etiq1
	colnames(mdist)<-etiq2
    return(mdist)
}
