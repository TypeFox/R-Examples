choose.thresh.nbedges <- function(cor.mat,  var.ind.mat = 0, n.ind = 0,  thresh = 0.05,  nb.edges = 405, test.method = "gaussian",  proc.length = 518, num.levels, use.tanh = FALSE, max.iter = 10)
{
  n.regions<-dim(cor.mat)[1]

  fin<-2*nb.edges/n.regions

  minR<-0
  maxR<-1
  iter<-0

  n.sup<-(minR+maxR)/2

  adj.mat<-const.adj.mat(cor.mat, thresh = thresh, sup = n.sup, proc.length = proc.length, n.ind = n.ind,use.tanh=use.tanh,var.ind.mat = var.ind.mat,num.levels=num.levels)

  mean.deg.tmp<-sum(adj.mat)/n.regions

  arret<-0

  while((abs(mean.deg.tmp-fin)>(3*10^{-2}))&&(arret==0)){
    iter<-iter+1
    if((mean.deg.tmp-fin)>0){
      minR<-n.sup
      n.sup<-(minR+maxR)/2
      adj.mat<-const.adj.mat(cor.mat, thresh = thresh, sup = n.sup, proc.length = proc.length, n.ind = n.ind,use.tanh=use.tanh,var.ind.mat = var.ind.mat,num.levels=num.levels)
      mean.deg.tmp<-sum(adj.mat)/n.regions
    }
    if((mean.deg.tmp-fin)<0){
      maxR<-n.sup
      n.sup<-(minR+maxR)/2
      adj.mat<-const.adj.mat(cor.mat, thresh = thresh, sup = n.sup, proc.length = proc.length, n.ind = n.ind,use.tanh=use.tanh,var.ind.mat = var.ind.mat,num.levels=num.levels)
      mean.deg.tmp<-sum(adj.mat)/n.regions
    }
    if(iter>=max.iter) arret<-1
  }

  sup.fin<-n.sup
  return(sup.fin)
}
