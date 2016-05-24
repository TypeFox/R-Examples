bestmatch <-
function(rsize,size){
        z<-cbind(c(rsize,size),c(rep(1,length(rsize)),rep(0,length(size))))
        z<-z[order(z[,1]),]
        mymatch<-rep(NA,length(size))
        osize<-order(size)
        orsize<-order(rsize)
        mymatch[osize]<-orsize[cumsum(z[,2])[z[,2]==0]]
        return(mymatch)
}
