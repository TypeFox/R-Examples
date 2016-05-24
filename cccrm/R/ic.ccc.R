ic.ccc <-
function(ccc,dev,S,alpha){
  
se.ccc <- sqrt(dev%*%S%*%t(dev))
z<-0.5*log((1+ccc)/(1-ccc))
se.z<-sqrt(  (se.ccc^2)/(((1+ccc)^2)*((1-ccc)^2))  )
ic.z=z+c(-1,1)*qnorm(1-alpha/2)*se.z
ic.icc=(exp(2*ic.z)-1)/(exp(2*ic.z)+1)
result<-c(ccc,ic.icc,se.ccc,z,se.z)
conf.lab=paste((1-alpha)*100,"%",sep="")

names(result)<-c("CCC",paste("LL CI",conf.lab),paste("UL CI",conf.lab),"SE CCC","Z","SE Z")
return(result)
}

