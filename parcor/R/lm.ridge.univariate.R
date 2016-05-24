lm.ridge.univariate<-function(x,y,lambda=0,scale=TRUE){
    r<-length(lambda)
    xs<-scale(x,scale=scale)
    ys<-scale(y,scale=scale)
    sx<-1
    sy<-1
    if (scale==TRUE){
      sx<-sd(x)
      sy<-sd(y)
    }
    b<-sum(xs*ys)/(sum(xs^2) +lambda)
    b<-b*sy/sx
    inter<-mean(y) - b*mean(x)
    coefficients<-cbind(inter,b)
    return(coefficients)
}