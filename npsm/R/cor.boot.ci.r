cor.boot.ci<-function(x,y,method='spearman',conf=0.95,nbs=3000) {
  quantile(
    boot::boot(cbind(x,y),function(d,i) cor(d[i,1],d[i,2],method=method),R=nbs)$t, 
    probs=c((1-conf)/2,1-(1-conf)/2) 
  )
}
