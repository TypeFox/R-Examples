const.adj.mat <- function(cor.mat, var.ind.mat = 0, n.ind = 0, thresh = 0.05, sup = 0, test.method="gaussian",proc.length,num.levels,use.tanh=FALSE)
{

  if(is.matrix(cor.mat)==FALSE) stop("The format of the data is not a matrix")

  n.regions<-dim(cor.mat)[1]

  adj.mat<-array(0,dim=c(n.regions,n.regions))

  pvalue.cor<-p.value.compute(cor.mat,var.ind.mat,test.method,proc.length=proc.length, sup=sup, num.levels=num.levels,use.tanh=use.tanh)
  pvalue.thresh<-compute.FDR(pvalue.cor,thresh)

  test.sign<-(pvalue.cor<=pvalue.thresh)
  l<-1
  for(k in 2:(n.regions)) {
    for(q in 1:(k-1)) {
      if((test.sign[l])==TRUE) {
        adj.mat[k,q]<-1
      }
      l<-l+1
    }
  }

  mat <- adj.mat
  adj.mat <- mat+t(mat)
  return(adj.mat)
}












