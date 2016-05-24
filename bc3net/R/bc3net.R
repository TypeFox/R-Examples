
bc3net <- function(dataset, boot=100, estimator="pearson", disc="equalwidth", mtc1=TRUE, alpha1=0.05, nullit=NA, null=c(),adj1="bonferroni", mtc2=TRUE, alpha2=0.05, adj2="bonferroni", weighted=TRUE, igraph=TRUE, verbose=FALSE){

    bnet = matrix(0,nrow=nrow(dataset),ncol=nrow(dataset))
    colnames(bnet)=rownames(dataset)
    rownames(bnet)=rownames(dataset)

    if(verbose==TRUE){
      cat("Bagging Network Ensemble\n",boot," (Bootstrap size)\n")
    }

    if(mtc1==TRUE){
      if(length(null)<1){
        null=makenull(dataset, nullit=nullit, estimator=estimator, disc=disc)         }
      
      maxn=max(null)
    }
    

    for(b in 1:boot)
    {
        # b=1
        if(verbose==TRUE){
           cat("Processing bootstrap dataset ", b,"of ",boot,"\n")
        }          
 
      expdata=dataset[,sample(1:ncol(dataset),size=ncol(dataset),replace=TRUE)]
      mim=mimwrap(expdata,estimator=estimator,disc=disc) 
            
      net=c3(mim)
      net[upper.tri(net)]=0
      edges<-which(net>0,arr.ind=TRUE)


      pval <- rep(0,nrow(edges))

       if(mtc1==TRUE){        

          for(e in 1:nrow(edges)){                
             G1 <- edges[e,1]
             G2 <- edges[e,2]          
             if(mim[G1,G2]>maxn){
                pval[e]=0
             }else{
                pval[e]=.getpval(mim[G1,G2],null)
             }
          }
      
          pval[pval==0]=1/length(null)
          padj=p.adjust(pval,method=adj1)
          
          
          for(e in which(padj>alpha1)){                      
             G1 <- edges[e,1]
             G2 <- edges[e,2]
             net[G1,G2]=0
             net[G2,G1]=0
          }

      } # mtc

         
      net = (net>0)*1
      bnet = bnet + net
      rm(net)
    } # for loop 1:boot

    
    bnet=bnet/boot
   
    # BINOMIAL TEST
    En=sum(bnet*boot) # number of observed gene pairs 
    p=nrow(bnet)      # number of genes
    E=((p*p)-p)/2     # number of possible gene pairs
    Et=E*boot         # number of possible gene pairs of the ensemble
    prob=En/Et
   
    edges <- which(bnet>0,arr.ind=TRUE)
    
    pval=c()
    padj=rep(0,nrow(edges))
   

    if(mtc2==TRUE & nrow(edges)>0){ # Binomial test for ensemble consensus rate
      
        if(verbose==TRUE){
           cat("Binomial tests for ",nrow(edges)," edges","\n")
        }


        for(e in 1:nrow(edges))
        {      
           G1   <- edges[e,1]
           G2   <- edges[e,2]
           freq <- bnet[G1,G2]*boot
           p    <- 1-pbinom(freq, size=boot, prob=prob)    
           pval <- c(pval,p)
        }

        padj=p.adjust(pval,method=adj2)   
        rject=length(which(padj>alpha2))
        en=length(which(padj<=alpha2))
        
        if(verbose==TRUE){
           cat("Deleting ",rject," edges\n")
           cat(en, " inferred edges\n")
        }
      
        for(e in which(padj>alpha2)){        
           G1 <- edges[e,1]
           G2 <- edges[e,2]
           bnet[G1,G2]=0
           bnet[G2,G1]=0      
        }
      
    }
      

     if(igraph==TRUE){
        if(weighted==TRUE){
          bnet=.mat2igraph(bnet)
        }else{
          bnet=.mat2igraph(bnet, weighted=FALSE)
        }
     }else{
          bnet=as.matrix(forceSymmetric(bnet,uplo="L"))
     }


     return(bnet)
}
