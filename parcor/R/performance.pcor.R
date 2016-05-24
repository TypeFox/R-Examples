performance.pcor <- function(inferred.pcor, true.pcor=NULL, fdr=TRUE, cutoff.ggm=0.8,verbose=FALSE,plot.it=FALSE)
{ 
    
    p <- ncol(inferred.pcor) # number of variables 
    if (fdr==TRUE){
    test.results <-ggm.test.edges(inferred.pcor,verbose=verbose,plot=plot.it)
    # extract network containing edges with prob > cutoff.ggm
    prob.results<-diag(p)
    for (i in 1:length(test.results$prob)){
        prob.results[test.results$node1[i],test.results$node2[i]]<-test.results$prob[i]
        prob.results[test.results$node2[i],test.results$node1[i]]<-test.results$prob[i]
        dim(prob.results)
    }
    net <- extract.network(test.results, cutoff.ggm=cutoff.ggm)
    adj <- diag(p)
    if (nrow(net)!=0) 
    { 
      for (j in 1:nrow(net))
      {
        adj[net[j,2], net[j,3]] <- 1
        adj[net[j,3], net[j,2]] <- 1
      }
    }        
  }
    if (fdr==FALSE)
  {
    adj <- inferred.pcor
    adj[adj!=0] <- 1
  }
     #
    num.selected <- ( sum(abs(adj)>0)-p )/2 # number of selected edges
    connectivity<-apply(adj,2,FUN=function(x) return(sum(x!=0)))-1
    positive.cor<-NULL
    if (num.selected>0){
        if (fdr==TRUE){
            positive.cor<-sum(net[,1]>0)/num.selected
        }
        if (fdr==FALSE){
            positive.cor<-((sum(inferred.pcor>0)-ncol(inferred.pcor))/2)/num.selected
        }
    }
    num.true<-power<-ppv<-NULL
    if (is.null(true.pcor)==FALSE){
        num.true <- ( sum(abs(true.pcor)>0)-p )/2 # number of true edges 
        num.true.positives <-  ( sum(abs(true.pcor*adj)>0)-p )/2
        num.false.positives<-num.selected-num.true.positives
        ## power
        power <- -Inf
        if (num.true>0) {
            power <- num.true.positives/num.true
        }
        ## true discovery rate (positive predictive value)
        ppv <- -Inf
        if (num.selected>0) {
            ppv <- num.true.positives/num.selected
        }
  }
  ## ROC curves and AUC score
  auc<-NA
  tpr<-fpr<-NA
  TPR<-FPR<-NA
  if ((fdr==TRUE) & (is.null(true.pcor)==FALSE)){
    xx<-seq(0,1,length=500) # scale of the x and y axis of the ROC curve
    true.binary=(abs(true.pcor)>0) # logical matrix indicating the true partial correlations
    predicted<-sym2vec(prob.results)
    if (var(predicted)>0){ # there is an error message if all probabilities are zero, weird
    if (plot.it==TRUE){
        plot.roc="ROC"
    }
    if (plot.it==FALSE){
        plot.roc=NA
    }
    myroc<-ROC(predicted,sym2vec(true.binary),plot=plot.roc)
    auc<-myroc$AUC # area under the curve
    TPR<-myroc$res[,1] # sensitivity =TPR
    #cat(paste("length of TPR ", length(TPR),"\n"))
    FPR<-1-myroc$res[,2] # 1-specifity =FPR
    #cat(paste("length of TPR ", length(TPR),"\n"))
    #TPR<-c(1,TPR,0)
    #FPR<-c(1,FPR,0)
    #yy<-approx(TPR,FPR,xout=xx,yleft=0,yright=1)
    #tpr<-xx
    #fpr<-yy$y
    # rest later
    }
  }
  if ((is.null(true.pcor)==FALSE)){
    tpr<-power
    #fpr<--Inf
    if (num.true>0){
    fpr<-num.false.positives/((p^2-p)/2 - num.true) # false positives/true false
   }
    }
    
  # also return the number of positive and negative correlations
  return(list(num.selected=num.selected, power=power,TPR=TPR,FPR=FPR, tpr=tpr,fpr=fpr,ppv=ppv,adj=adj,connectivity=connectivity,positive.cor=positive.cor,auc=auc))
}
 # fpr and tpr are for the optimal parameters
 # FPR and TPR are vectors
