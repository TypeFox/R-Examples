Summarynumeric0<-function(x,digits=max(3, getOption("digits")-3)){
T1<-summary(na.omit(x))
SD<-signif(sd(na.omit(x)),digits)
T2<-c(T1,SD)
names(T2)[1]<-"Min."
names(T2)[2]<-"Q1"
names(T2)[3]<-"Q2"
names(T2)[4]<-"M"
names(T2)[5]<-"Q3"
names(T2)[6]<-"Max."
names(T2)[7]<-"S"
T2}
