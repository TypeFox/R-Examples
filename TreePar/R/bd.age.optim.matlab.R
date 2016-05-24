bd.age.optim.matlab <- function(x, lambdainit=1, kinit=1, thetainit=2, sampling=1, model="G", root=1, inputformat=0,precision=4,numgridpts=500, path,matfilename="setup"){
    if (root == 1){Lcond<-"C"} else {Lcond<-"S"}	
    Param <- paste("'",lambdainit, kinit, thetainit,"'")
    if (model=="E") 	{    Param <- paste("'",lambdainit, thetainit,"'")}
    	if (inputformat==0){
    		create.mat(x,numgridpts, path,matfilename)
    	}
    runCmd <- paste("./run_MaxLFcn.sh ",path, matfilename, "outputML", as.character(precision), Lcond, model, Param, as.character(sampling), sep=" ")	
    system(runCmd)
    out<-read.table("outputML.txt",sep=",")
    res<-c(out[[length(out)]])
    for (i in 1:(length(out)-1)){res<-c(res,out[[i]])}
    res
}