
const.cor.list <-function(data.mat, names.data = 0, method = "modwt" ,wf = "la8", n.levels = 4, boundary = "periodic", p.corr = 0.975, save.wave=FALSE, export.data=FALSE)
{

version<-2

if(is.matrix(data.mat)==FALSE) 
		stop("The format of the data is not a matrix")

n.regions<-dim(data.mat)[2]
proc.length<-dim(data.mat)[1]

wave.cor.list <- vector("list", (3*n.levels))
    names(wave.cor.list) <- c(paste("d", 1:n.levels, sep = ""), paste("lowerd", 1:n.levels, sep = ""),paste("upperd", 1:n.levels, sep = ""))

wave.cor.mat<-array(0,c(n.regions,n.regions,n.levels))
wave.cor.lower.mat<-array(0,c(n.regions,n.regions,n.levels))
wave.cor.upper.mat<-array(0,c(n.regions,n.regions,n.levels))


if(save.wave==TRUE){
	wave.trans.data<-vector("list",n.regions)

	class(wave.trans.data)<-"Wavelet data decomposition"
	attr(wave.trans.data, "version") <- version
	attr(wave.trans.data, "method") <- method
	attr(wave.trans.data, "wavelet") <- wf
	attr(wave.trans.data, "boundary") <- boundary
	attr(wave.trans.data, "proc.length") <- proc.length
	attr(wave.trans.data, "num.time.series") <- n.regions 




	Nx<-dim(data.mat)[1]
 	for(i in 1:n.regions){
		wave.trans.data[[i]]<-wave.trans(data.mat[,i],method,wf,n.levels,boundary)
			     }
	for(i in 1:n.regions){
		for(j in i:n.regions){


		 tmp <- wave.correlation(wave.trans.data[[i]],wave.trans.data[[j]],Nx,p=p.corr)
   		 wave.cor.mat[i,j,] <- tmp[1:n.levels,1]
   		 wave.cor.lower.mat[i,j,] <- tmp[1:n.levels,2]
   		 wave.cor.upper.mat[i,j,] <- tmp[1:n.levels,3]
				
				}}
}

if(save.wave==FALSE){

 for(i in 1:n.regions){
	
	x<-data.mat[,i]  
 	  x.bw<-wave.trans(x,method,wf,n.levels,boundary)
   
		for(j in i:n.regions){
    
  #compute all the coefficients of the wave correlation
  #matrix (only the upper diagonal terms)
  

  y<-data.mat[,j]

Nx<-length(x)
Ny<-length(y)

if(Ny != Nx) stop("The two time series must have the same length")

 
   y.bw<-wave.trans(y,method,wf,n.levels,boundary)
   tmp <- wave.correlation(x.bw,y.bw,Nx,p=p.corr)
   wave.cor.mat[i,j,] <- tmp[1:n.levels,1]
   wave.cor.lower.mat[i,j,] <- tmp[1:n.levels,2]
   wave.cor.upper.mat[i,j,] <- tmp[1:n.levels,3]

       }
    }
}

# To produce the symmetric part (so far only upper diagonal terms written)

for(n in 1:n.levels){
          wave.cor.mat[,,n] <- wave.cor.mat[,,n] + t(wave.cor.mat[,,n])
	  wave.cor.lower.mat[,,n] <- wave.cor.lower.mat[,,n] + t(wave.cor.lower.mat[,,n])
	  wave.cor.upper.mat[,,n] <- wave.cor.upper.mat[,,n] + t(wave.cor.upper.mat[,,n])
	  for(i in 1:n.regions) wave.cor.mat[i,i,n] <- 1 

## putting 1 on leading diagonal
        }

for(n in 1:n.levels){
wave.cor.list[[n]]<-wave.cor.mat[,,n]
wave.cor.list[[n+n.levels]]<-wave.cor.lower.mat[,,n]
wave.cor.list[[n+2*n.levels]]<-wave.cor.upper.mat[,,n]
}

class(wave.cor.list)<-"Wave Correlation"
attr(wave.cor.list, "version") <- version
attr(wave.cor.list, "method") <- method
attr(wave.cor.list, "wavelet") <- wf
attr(wave.cor.list, "boundary") <- boundary
attr(wave.cor.list, "proc.length") <- proc.length
attr(wave.cor.list, "num.time.series") <- n.regions 


if(export.data==TRUE){
save.cor.txt(wave.cor.list)
}




if(save.wave==TRUE){
list(wave.cor.list = wave.cor.list, wave.trans.data = wave.trans.data)
}else{
return(wave.cor.list)
}

}






