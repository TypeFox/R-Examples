
pbc.sample=function(){
Z=list()
D=CompCase(pbc[,c(2:4,10:14)]); 
Z$time=D$time/365.25;
Z$status=as.numeric(D$status==2);
Z$group=as.numeric(D$trt==2);
Z$covariates=as.matrix(D[,4:8]);
Z
}
