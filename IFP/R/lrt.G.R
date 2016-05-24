
lrtG<-function(n.fp,genoT,genoC){  ## true order

 n<-ncol(genoC)
 ncon<-ncol(genoT)
 n.poly<-nrow(genoT)
 choices<-0
 for(i in 1:n.fp) choices<-choices+choose(n.poly,i)
 af<-allele.freq.G(genoT)
 caf<-allele.freq.G(genoC)
  genoT[which(af>caf),][genoT[which(af>caf),]==1]<-2
  genoT[which(af>caf),][genoT[which(af>caf),]==0]<-1
  genoT[which(af>caf),][genoT[which(af>caf),]==2]<-0
  genoC[which(af>caf),][genoC[which(af>caf),]==1]<-2
  genoC[which(af>caf),][genoC[which(af>caf),]==0]<-1
  genoC[which(af>caf),][genoC[which(af>caf),]==2]<-0

 temp<-.C("LRT",as.integer(n.fp),as.integer(n.poly),as.integer(as.vector( as.vector(t(as.matrix(genoT))))),as.integer( as.vector(t(as.matrix(genoC)))),as.integer(ncon),as.integer(n[1]/2),pMMcV=double(choices),df=double(choices))
 lr2<-temp$pMMcV[1:choices]
 df<-temp$df 
 rm(temp)

info<-array(NA,c(choices,n.fp))
id<-1

for(j in 1:n.fp){

no<-1:j
no.t<-no

for (i in 1:choose(n.poly,j)){

 info[id,1:j]<-no.t
 id<-id+1

 z<-sum(no.t==n.poly-j+no)
 if (z!=j) {
  no.t[j-z]<-no.t[j-z]+1
  if (z>0) for (j in (j-z+1):j) no.t[j]<-no.t[j-1]+1
 }

}
} 

 gc()
 cbind(1:choices,info,pchisq((lr2),df),lr2, df)

}

