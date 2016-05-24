read.var.txt <- function()
{


tmp<-read.table("wave_var_mat_level_1.txt",skip=0,nrows=1)

version<-tmp[[1]]

tmp<-read.table("wave_var_mat_level_1.txt",skip=1,nrows=1)

if(tmp[[2]]!= "Variance") stop("The files are not of class variance matrix")

tmp<-read.table("wave_var_mat_level_1.txt",skip=2,nrows=1)
method<-tmp[[1]]
wf<-tmp[[2]]
boundary<-tmp[[3]]
n.levels<-tmp[[4]]
if(version==2){ 
	proc.length<-tmp[[5]]
	n.regions<- tmp[[6]]
	
}
if(version>2) stop("please use a newer version of the package")


wave.var.list <- vector("list", (3*n.levels))
    names(wave.var.list) <- c(paste("d", 1:n.levels, sep = ""), paste("lowerd", 1:n.levels, sep = ""),paste("upperd", 1:n.levels, sep = ""))

class(wave.var.list)<-"Wave Variance"
attr(wave.var.list, "version") <- version
attr(wave.var.list, "method") <- method
attr(wave.var.list, "wavelet") <- wf
attr(wave.var.list, "boundary") <- boundary
if(version!=1){
attr(wave.var.list, "proc.length") <- proc.length
attr(wave.var.list, "num.time.series") <- n.regions
}


for(i in 1:n.levels){

name.txt<-paste("wave_var_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")


wave.var.list[[i]]<-as.matrix(read.table(name.txt,skip=3))

name.txt<-paste("wave_var_lower_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")

wave.var.list[[i+n.levels]]<-as.matrix(read.table(name.txt,skip=3))

name.txt<-paste("wave_var_upper_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")

wave.var.list[[i+2*n.levels]]<-as.matrix(read.table(name.txt,skip=3))

}

return(wave.var.list)

}

