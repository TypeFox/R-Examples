luo<-function(tree)
{
S<-tree$S
R<-tree$ssr                #R(i) on noden i ssr eli -log likeli

alknodlkm<-length(tree$left)
#leafloc<-findleafs(tree$left,tree$right)

N<-matrix(0,alknodlkm,1)   #number of leaves in the tree whose root is i
p<-matrix(0,alknodlkm,1)   #parent
g<-matrix(0,alknodlkm,1)   #(R(i)-S(i))/(N(i)-1), R(i) on noden i ssr
G<-matrix(0,alknodlkm,1)   #min{g(t),G(l(t)),G(r(t))}

t<-alknodlkm

while (t>=1){
  if (tree$left[t]==0){  #l(t)=0 eli ollaan lehdessa 
     N[t]<-1
     G[t]<-NA                 #\infty
     g[t]<-NA
  }
  else{     #if (!is.na(leafloc[t])){
     p[tree$left[t]]<-t       #p[t+1]<-t
     p[tree$right[t]]<-t      #p[endpoint(tree,t)+1]<-t
     N[t]<-N[tree$left[t]]+N[tree$right[t]]
     g[t]<-(R[t]-S[t])/(N[t]-1)
     G[t]<-omamindelt(g[t],omamindelt(G[tree$left[t]],G[tree$right[t]]))
  }
  t<-t-1
}

return(p=p,G=G,g=g,N=N)
}



