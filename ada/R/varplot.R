"varplot" <-
function(x,plot.it=TRUE,type=c("none","scores"),max.var.show=30,...){
  if(class(x)!="ada"){
    stop("Object must be of type 'ada'")
  }
  if(missing(type)){
    type="none"
  }
  iter<-x$iter
  nm<-x$names
  vec<-rep(0,length(nm))
  p=x$dim[2]
  g1<-function(i,obj){
    if(dim(obj[[i]]$frame)[1]<2){
      return(rep(0,p))
    }
    imp<-obj[[i]]$splits[,3]^2
    vals<-match(row.names(obj[[i]]$splits),nm)
    vec=rep(0,p)
    vec[vals]<-imp
    vec
  }
  vec<-1/iter*sqrt(apply(sapply(1:iter,function(i)g1(i,x$model$trees)),1,sum))
 
  vars<- order(vec,decreasing=TRUE)
  n<-length(vec)
  max.v=max.var.show
  if(p<max.v)
    max.v=p
  if(plot.it==TRUE){
    dotchart(vec[vars[max.v:1]],labels=nm[vars[max.v:1]],xlab="Score",main="Variable Importance Plot")
  }
  if(type=="scores"){
    vars=vars[1:max.v]
    t1<-vec[vars]
    attr(t1,"names")<-nm[vars]
    return(t1)
 }
}

