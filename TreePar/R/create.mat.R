create.mat <- function(x,numgridpts=500, path,matfilename="setup"){
	write.table(x,file="dataBranching.txt",row.names=FALSE,col.names=FALSE)	
    runCmd <- paste("./run_ReadTreeFcn.sh ",path, "dataBranching", matfilename, as.character(numgridpts), sep=" ")
    res<-system(runCmd)
}    		