CompCase=function(mydata){
sum(is.na(apply(mydata,1,mean)))
mydata<-mydata[!is.na(apply(mydata,1,mean)),]
} 


#mydata=CompCase(pbc[,c(2:4,10:14)])
#head(mydata)
