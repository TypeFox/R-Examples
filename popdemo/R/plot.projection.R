plot.projection<-function(x,labs=TRUE,...){
if(is.list(x)) N<-x$N else(N<-x)
if(length(dim(N))==0){
    len<-length(N)
    Time.intervals<-0:(len-1)
    Population.density<-N
    plot(Time.intervals,Population.density,type="l",...)
}
if(length(dim(N))==2){
    len<-dim(N)[1]
    Time.intervals<-c(0,(len-1))
    Population.density<-c(min(N),max(N))
    plot(Time.intervals,Population.density,type="n",...)
    for(i in 1:dim(N)[2]){
        options(warn=-1)
        lines(0:(len-1),N[,i],...)
        if(labs) text(len-1,N[len,i],paste("Bias",i),adj=c(1,1),...)
        options(warn=0)
    }
}}
