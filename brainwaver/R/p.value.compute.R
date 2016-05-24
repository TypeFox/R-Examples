p.value.compute <- function(test.mat,var.ind.mat=0,n.ind = 0, test.method="gaussian",proc.length,sup,num.levels,use.tanh=FALSE)
{
  n.regions<-dim(test.mat)[1]

  l<-1
  pvalue.cor<-rep(0,(n.regions-1)*n.regions/2)

  N<-trunc(proc.length/2^(num.levels))

  for(k in 2:(n.regions)){
    for(q in (1):(k-1)){
      if(use.tanh==FALSE){
        f<-atanh(test.mat[k,q])
      }else{
        f<-test.mat[k,q]
      }

      if(test.method=="t.test") var.f<-var.ind.mat[k,q]

      if(test.method=="gaussian"){
        if(sup == 0){
          if(f >= 0){
            pvalue.cor[l]<-2*pnorm(-(f-atanh(sup))*sqrt(N-3))
          }else{
            pvalue.cor[l]<-2*pnorm((f-atanh(-sup))*sqrt(N-3))
          }
        }else{
          if(f>=0){
            pvalue.cor[l]<-pnorm(-(f-atanh(sup))*sqrt(N-3))
          }else{
            pvalue.cor[l]<-pnorm((f-atanh(-sup))*sqrt(N-3))
          }
        }
      }
      if(test.method=="t.test"){
        stat.val<-sqrt(n.ind)*(f-atanh(sup))/sqrt(var.f)
        stat.val.minus<-sqrt(n.ind)*(f-atanh(-sup))/sqrt(var.f)
        if(sup == 0){
          if(f>=0){
            pvalue.cor[l]<-2*pt(-stat.val,df=(n.ind-1))
          }else{
            pvalue.cor[l]<-2*pt(stat.val.minus,df=(n.ind-1))
          }
        }else{
          if(f>=0){
            pvalue.cor[l]<-pt(-stat.val,df=(n.ind-1))
          }else{
            pvalue.cor[l]<-pt(stat.val.minus,df=(n.ind-1))
          }
        }
      }
      l<-l+1
    }
  }
  return(pvalue.cor)
}