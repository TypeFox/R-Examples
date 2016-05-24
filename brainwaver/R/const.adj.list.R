const.adj.list <- function(wave.cor.list, wave.var.ind = 0, n.ind = 0, thresh = 0.05, sup = 0, test.method="gaussian",proc.length=-1,use.tanh=FALSE)
{

  if(class(wave.cor.list)!="Wave Correlation") stop("The type of the data is not correct, please check!")

  n.levels<-length(wave.cor.list)/3
  n.regions<-dim(wave.cor.list[[1]])[1]

  method <- attr(wave.cor.list, "method")
  wf <- attr(wave.cor.list, "wavelet")
  boundary <- attr(wave.cor.list, "boundary")

  if(proc.length==-1) proc.length<-attr(wave.cor.list, "proc.length")
  if(proc.length<0) stop("Error : the length of the time series is negative")

  wave.adj.list <- vector("list", (n.levels))
  names(wave.adj.list) <- paste("d", 1:n.levels, sep = "")
  class(wave.adj.list)<-"Wave Adjacency matrix"
  attr(wave.adj.list, "method") <- method
  attr(wave.adj.list, "wavelet") <- wf
  attr(wave.adj.list, "boundary") <- boundary
  attr(wave.adj.list, "proc.length") <- proc.length

# Loop on the scales

  for(i in 1:n.levels){
    if(wave.var.ind == 0){
      wave.adj.list[[i]]<-const.adj.mat(cor.mat=wave.cor.list[[i]], n.ind = n.ind, thresh = thresh, sup = sup, test.method= test.method, proc.length = proc.length, use.tanh=use.tanh,num.levels=i)
    }else{
      if(class(wave.var.ind != "Wavelet variance")) stop("Incorrect type of data for the variance")
      wave.adj.list[[i]]<-const.adj.mat(cor.mat=wave.cor.list[[i]], var.ind.mat = wave.var.ind[[i]], n.ind = n.ind, thresh = thresh, sup = sup, test.method= test.method, proc.length = proc.length, use.tanh=use.tanh,num.levels=i)
    }
  }

  return(wave.adj.list)
}












