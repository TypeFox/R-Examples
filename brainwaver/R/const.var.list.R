
const.var.list <-function(data.mat, names.data = 0, method = "modwt" ,wf = "la8", n.levels = 4, boundary = "periodic", save.wave=FALSE, export.data=FALSE)
{

version<-2

if(is.matrix(data.mat)==FALSE) 
		stop("The format of the data is not a matrix")

n.regions<-dim(data.mat)[2]
proc.length<-dim(data.mat)[1]


wave.var.list <- vector("list", (3*n.levels))
   names(wave.var.list) <- c(paste("d", 1:n.levels, sep = ""), paste("lowerd", 1:n.levels,sep = ""),paste("upperd", 1:n.levels, sep = ""))


wave.var.mat<-array(0,c(n.regions,n.levels))
wave.var.lower.mat<-array(0,c(n.regions,n.levels))
wave.var.upper.mat<-array(0,c(n.regions,n.levels))




if(save.wave==TRUE){
	wave.trans.data<-vector("list",n.regions)

	class(wave.trans.data)<-"Wavelet data decomposition"
	attr(wave.trans.data, "version") <- version
	attr(wave.trans.data, "method") <- method
	attr(wave.trans.data, "wavelet") <- wf
	attr(wave.trans.data, "boundary") <- boundary
	attr(wave.trans.data, "proc.length") <- proc.length
	attr(wave.trans.data, "num.time.series") <- n.regions 



 	for(i in 1:n.regions){
		wave.trans.data[[i]]<-wave.trans(data.mat[,i],method,wf,n.levels,boundary)
                wave.var.mat[i,] <- wave.variance(wave.trans.data[[i]])[1:n.levels,1]
		wave.var.lower.mat[i,] <- wave.variance(wave.trans.data[[i]])[1:n.levels,2]
                wave.var.upper.mat[i,] <- wave.variance(wave.trans.data[[i]])[1:n.levels,3]


	}
}


if(save.wave==FALSE){
 	for(i in 1:n.regions){
	x<-data.mat[,i]  
 	  x.bw<-wave.trans(x,method,wf,n.levels,boundary)
          tmp <- wave.variance(x.bw)
          wave.var.mat[i,] <- tmp[1:n.levels,1]
          wave.var.lower.mat[i,] <- tmp[1:n.levels,2]
          wave.var.upper.mat[i,] <- tmp[1:n.levels,3]
	
      }
}

for(n in 1:n.levels){

wave.var.list[[n]]<-wave.var.mat[,n]
wave.var.list[[n+n.levels]]<-wave.var.lower.mat[,n]
wave.var.list[[n+2*n.levels]]<-wave.var.upper.mat[,n]



}


class(wave.var.list)<-"Wave Variance"
attr(wave.var.list, "version") <- version
attr(wave.var.list, "method") <- method
attr(wave.var.list, "wavelet") <- wf
attr(wave.var.list, "boundary") <- boundary
attr(wave.var.list, "proc.length") <- proc.length
attr(wave.var.list, "num.time.series") <- n.regions 


if(export.data==TRUE){
save.var.txt(wave.var.list)
}

if(save.wave==TRUE){
list(wave.var.list = wave.var.list, wave.trans.data = wave.trans.data)
}else{
return(wave.var.list)
}

}