CoreSetterCombined <-
function(genos,size,criterion=c("HET","MTK"),save=NULL,power=10,permutations=100,print=TRUE,mat=NULL){
  genos<-as.matrix(genos)
  n.genos<-ncol(genos)
  mode(genos)<-"numeric"
  result.list<-rep(0,permutations)
  ##Build a list of random subsets to be optimized
  random.lists<-matrix(0,nrow=size,ncol=permutations)
  if(is.null(save)){ 
    for(i in 1:permutations){
      random.lists[,i]<-sample(colnames(genos),size)
    }
  }else{
    random.lists[1:length(save),]<-save
    for(i in 1:permutations){
      random.lists[(length(save)+1):size,i]<-sample(colnames(genos),size-length(save))
    }
  }
  better.randoms<-matrix(nrow=size,ncol=permutations)
  ##Optimize List with HET
  if(criterion=="HET"){
    for(i in 1:permutations){
      better.randoms[,i]<-CoreSetOptimizer(genos,random.lists[,i],criterion="HET",save=save,print=FALSE) #Use CoreSetOptimizer to optimize subsets
      result.list[i]<-HET(genos[,better.randoms[,i]])
    }
    print(c("Best HET:",max(result.list)))
    best<-which.max(result.list)
  }
  ##Optimize list with MTK
  if(criterion=="MTK"){
    if(is.null(mat)){
      mat<-A.mat(t(genos)) #Building Kinship Matrix
    }
    mat.adj<-(2*(mat-min(mat))/(max(mat)-min(mat)))^power #Transforming Kinship Matrix to Emphasize Similar Genotypes
    for(i in 1:n.genos){
      mat.adj[i,i]<-0
    }
    for(i in 1:permutations){
      better.randoms[,i]<-CoreSetOptimizer(genos,random.lists[,i],criterion="MTK",mat=mat,save=save,print=FALSE) #Use CoreSetOptimizer with untransformed matrix to optimize subsets
      result.list[i]<-sum(mat.adj[better.randoms[,i],better.randoms[,i]])/(n.genos^2)
    }
    print(c("Lowest MTK:",min(result.list)))
    best<-which.min(result.list)
  }
  return(better.randoms[,best])
}

