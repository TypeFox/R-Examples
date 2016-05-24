fill.values <-
function(fromto,values,ldest){
        dest<-rep(0,ldest)
        z<-cbind(c(fromto[,1],fromto[,2]+1),c(values,-values))
        z<-z[order(z[,1]),]
        z[,2]<-cumsum(z[,2])
        z<-z[nrow(z)-rev(match(unique(rev(z[,1])),rev(z[,1])))+1,]
        #if(z[nrow(z),1]>ldest)z<-z[-nrow(z),]
        if(z[nrow(z),1]>ldest)z<-matrix(ncol=ncol(z),data=z[-nrow(z),])
        dest[z[,1]]<-z[,2]
        zz<-rep(0,length(dest))
        zz[z[,1]]<-c(0,z[-nrow(z),2])
        return(cumsum(dest-zz))
}
