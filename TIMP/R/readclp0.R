"readclp0" <-
function (filenm) 
{
    cmat <- read.table(filenm, skip = 1)
    clist <- list()
    for(i in 1:nrow(cmat)) {
	clist[[i]] <- list(low=cmat[i,2], high = cmat[i,3], 
	comp = cmat[i,1])
    }
    clist
}
