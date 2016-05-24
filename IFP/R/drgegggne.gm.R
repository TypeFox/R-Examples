drgegggne<-function(fdg,frg,fdgg,frgg,fdge,frge,eg,e){

temp1<-drgegne(fdg,frg,fdge,frge,eg,e)

temp<-drggn(fdgg,frgg)
temp2<-temp[[2]]

temp<-temp[[2]]
for(i in 1:9){
 temp[i,]<-c(temp1[i,1]*temp2[i,1],temp1[i,1]*temp2[i,2]+temp1[i,2]*temp2[i,1]+temp1[i,2]*temp2[i,2]/2,(temp1[i,1]+temp1[i,2])*temp2[i,3]+temp1[i,2]*temp2[i,2]/2+temp1[i,3])
}

temp

}
