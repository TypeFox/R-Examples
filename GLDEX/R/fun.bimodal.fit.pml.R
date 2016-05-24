"fun.bimodal.fit.pml" <-
function(data1,data2,first.fit,second.fit,prop,param1,param2,selc1,selc2){

index1<-switch(selc1,"rs"=1,"fmkl"=2,"star"=3)
index2<-switch(selc2,"rs"=1,"fmkl"=2,"star"=3)

result<-optim(c(first.fit[,index1],second.fit[,index2],prop), optim.fun4 ,data1 = data1, data2=data2,param1 = param1, 
param2=param2, control = list(maxit = 2000))


return(result)
}

