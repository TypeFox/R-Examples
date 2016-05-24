drgegne<-function(fdg,frg,fdge,frge,eg,e){

#g+ge+e

temp1<-drgene(fdge,frge,eg,e)

temp<-drgn(fdg,frg)
temp2<-temp[[2]]

temp<-temp[[2]]
for(i in 1:9){
 temp[i,]<-c(temp1[i,1]*temp2[i,1],temp1[i,1]*temp2[i,2]+temp1[i,2]*temp2[i,1]+temp1[i,2]*temp2[i,2]/2,(temp1[i,1]+temp1[i,2])*temp2[i,3]+temp1[i,2]*temp2[i,2]/2+temp1[i,3])
}

temp

}

