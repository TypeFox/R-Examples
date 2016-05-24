c3mtc <-
function(dataset,null=NULL,mtc=TRUE,adj="bonferroni",alpha=0.05,nullit=NA, estimator="pearson", disc="none",adjacency=FALSE, igraph=TRUE){

      
      mat.mim=mimwrap(dataset, estimator=estimator, disc=disc)      
      net=c3(mat.mim)
      
      net[upper.tri(net)]=0
      edges<-which(net>0,arr.ind=TRUE) # node indices
      
      pval <- matrix(0,nrow=nrow(edges),ncol=1)

      if(mtc==TRUE){
          if(is.null(null)){
            null=makenull(dataset, nullit=nullit, estimator=estimator, disc=disc)
          }
          maxn=max(null)
          lnull=length(null)
          
          for(e in 1:nrow(edges)){                
            G1 <- edges[e,1]
            G2 <- edges[e,2]          
            if(mat.mim[G1,G2]>maxn){
              pval[e]=1/lnull
            }else{
              # pval[e]=1-sum(mat.mim[G1,G2]>null)/length(null)
              pval[e]=.getpval(mat.mim[G1,G2],null)
            }
          }
          
          padj=p.adjust(pval,method=adj)
    
          for(e in which(padj>alpha)){
             G1 <- edges[e,1]
             G2 <- edges[e,2]        
             net[G1,G2]=0
             net[G2,G1]=0      
          }
          
    }else{
       net=c3(mat.mim)
    }

      
     if(igraph==TRUE){
       if(adjacency==TRUE){
          net=.mat2igraph(net, weighted=FALSE)
        }else{
          net=.mat2igraph(net)
        }        
     }else{
          net=as.matrix(forceSymmetric(net,uplo="L"))
          if(adjacency==TRUE){
             net=(net>0)*1
          }  
     }

     return(net)
}

