"classnum"<-function(x, breaks="Sturges"){
    if (is.numeric(breaks)& length(breaks)==1) ncl<-breaks
    
    if (is.character(breaks)) {
            breaks <- match.arg(tolower(breaks), c("sturges", 
                "fd", "freedman-diaconis", "scott"))
                ncl <- switch(breaks, sturges = nclass.Sturges(x), 
                "freedman-diaconis" = , fd = nclass.FD(x), scott = nclass.scott(x), 
                stop("Unknown breaks algorithm"))
                }
    if ((is.numeric(breaks)& length(breaks)==1)|(is.character(breaks))){
        etendue<-range(x)
        breaks1<-seq(etendue[1],etendue[2],l=ncl+1)[c(-1,-(ncl+1))]
    }
    
    if (is.numeric(breaks)& length(breaks)>1) {
        breaks1<-sort(breaks)
        ncl<-length(breaks)+1
     }
    
    index<-rep(NA,length(x))
    for (i in 1:(length(breaks1))){
        index[x<=breaks1[i]& is.na(index)]<-i
    }
    
    labtmp<-round(breaks1,2)
    lab<-paste(labtmp[1:length(labtmp)-1],labtmp[2:length(labtmp)],sep=", ")
    lab<-c(paste("<=",labtmp[1],sep=""),lab,paste(">",labtmp[length(labtmp)]))

    index[is.na(index)]<-ncl
    class(index)<-"clnum"
    attributes(index)$breaks<-breaks1
    attributes(index)$labels<-lab
    index
}

print.clnum<-function(x,...){
attributes(x)<-NULL
print(x)
}
