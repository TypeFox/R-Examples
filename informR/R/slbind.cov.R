slbind.cov<-function(covar, statslist, type = 1, ...) 
{
    if (length(statslist) != length(covar)) {
        stop("\n Both objects must be the same length.")
    }
    if(any(is.na(unlist(covar)))){
    stop("\n NAs are not allowed in covariates. Impute!")
    }
    newstatsl <- statslist
    for (i in 1:length(statslist)) {
        sldim<-dim(statslist[[i]][[type]])
        slnames<-dimnames(statslist[[i]][[type]])
        cvdim<-length(covar[[i]])
        cvs<-do.call("cbind",covar[[i]])
        ndimnames<-list(slnames[[1]],slnames[[2]],colnames(cvs))
        newstatsl[[i]][[type]] <- abind(statslist[[i]][[type]], array(0,dim=c(sldim[[1]],sldim[[2]],cvdim),dimnames=ndimnames),make.names=TRUE,...)
        for(j in 1:cvdim){
        newstatsl[[i]][[type]][,,j+sldim[3]]<-cvs[,j]
        }
    }
    newstatsl
}
