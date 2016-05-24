Mode<-function(x){
x<-na.omit(x)
step1<-as.factor(x)
m<-summary(step1)
m1<-which(m==max(m))
if(all(m==1)) print("No mode")
else
as.numeric(names(m1))}
