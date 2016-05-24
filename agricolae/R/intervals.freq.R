`intervals.freq` <-
function(x){
if(class(x)[1]=="graph.freq" | class(x)[1]=="histogram")breaks<-x$breaks
if(class(x)[1]=="numeric" | class(x)[1]=="integer"  )breaks<-x
n<-length(breaks)-1
classes<-rep(0,2*n)
dim(classes)<-c(n,2)
for (i in 1:n) {
classes[i,1]<-breaks[i]
classes[i,2]<-breaks[i+1]  
}
colnames(classes)<-c("lower","upper")
return(classes)
}

