LikAge <- function(x, lambda, k, theta, sampling, root=1, inputformat=0,precision=4,numgridpts=500, path,matfilename="setup"){
    if (root == 1){Lcond<-"C"} else {Lcond<-"S"}	
    Param <- paste("'",lambda, k, theta,"'")
    if (inputformat==0){
    		  create.mat(x,numgridpts, path,matfilename)
    }
    runCmd <- paste("./run_EvalLFcn.sh ",path, matfilename, "outputLik", as.character(precision), Lcond, Param, as.character(sampling), sep=" ")	
    res<-system(runCmd)
	res<-read.table("outputLik.txt")[[1]]
    res
}