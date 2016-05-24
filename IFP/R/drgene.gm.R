drgene<-function(fdg,frg,eg,e){

temp1<-drgen(fdg,frg,eg)
temp<-temp1[[2]]

for(i in 1:9){
 temp[i,]<-c(temp[i,1]*(1-e)^2,temp[i,2]*((1-e)^2+e*(1-e))+temp[i,1]*2*e*(1-e),temp[i,3]+temp[i,2]*(e^2+e*(1-e))+temp[i,1]*e^2)
}

temp

}
