sevennum <-
function(x){
x<-na.omit(x)
result<-matrix(c(rep(0,7)),ncol=7)
result[,1]<-min(x);result[,2]<-quantile(x,0.1);result[,3]<-quantile(x,0.25)
result[,4]<-quantile(x,0.5);result[,5]<-quantile(x,0.75)
result[,6]<-quantile(x,0.90);result[,7]<-max(x)
dimnames(result)<-list(c(" "),c("Min.","10th Quan.","25th Quan.","50th Quan.","75th Quan.",
"90th Quan.","Max."))
result
}

