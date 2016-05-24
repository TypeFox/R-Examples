CoreSetOptimizer<-
function(genos,subset,criterion=c("HET","MTK"),mat=NULL,save=NULL,power=10,print=TRUE){
  genos<-as.matrix(genos)
  mode(genos)<-"numeric"
  subset.genos<-length(subset)
  extra<-setdiff(colnames(genos),subset)
  extra.genos<-length(extra)
  n.genos<-ncol(genos)
  m<-nrow(genos)
  i<-0
  j<-1
  if(criterion=="HET"){
    het<-HET(genos[,subset])
    if(print==TRUE) print(c("Starting HET:",het))
    mat.a<-matrix(0,nrow=m,ncol=n.genos)  
    mat.b<-matrix(0,nrow=m,ncol=n.genos)
    mat.a[genos==1]<-1
    mat.a[genos==0]<-0.5
    mat.b[genos==-1]<-1
    mat.b[genos==0]<-0.5
    names<-colnames(genos)
    colnames(mat.a)<-names
    colnames(mat.b)<-names
    mat.a.extras<-mat.a[,extra]
    mat.b.extras<-mat.b[,extra]
    mat.a.short<-mat.a[,subset]
    mat.b.short<-mat.b[,subset]
    while(i<(subset.genos)){
      if(!subset[j]%in%save){
        a<-rowSums(mat.a.short)-mat.a.short[,j]
        b<-rowSums(mat.b.short)-mat.b.short[,j]
        temp.a<-a+mat.a.extras
        temp.b<-b+mat.b.extras
        temp.log<-colMeans(1-(temp.a^2+temp.b^2)/(temp.a+temp.b)^2,na.rm=TRUE)
        if(max(temp.log)>het){
          i<-0
          het<-max(temp.log)
          top<-which.max(temp.log)
          old<-subset[j]
          subset[j]<-extra[top]
          extra[top]<-old
          top.gen.a<-mat.a.extras[,top]
          top.gen.b<-mat.b.extras[,top]
          mat.a.extras[,top]<-mat.a.short[,j]
          mat.b.extras[,top]<-mat.b.short[,j]
          mat.a.short[,j]<-top.gen.a
          mat.b.short[,j]<-top.gen.b
          if(print==TRUE) print(het)
        }
      }
      if(j==subset.genos){
        j<-0
        if(print==TRUE) print("Cycle Completed, Returning to Beginning of Subset")
      }
      j<-j+1
      i<-i+1
    }
  }
  if(criterion=="MTK"){
    if(is.null(mat)){
      mat<-A.mat(t(genos))
    }
    mat.adj<-(2*(mat-min(mat))/(max(mat)-min(mat)))^power
    for(i in 1:n.genos){
      mat.adj[i,i]<-0
    }
    i<-0
    kin<-sum(mat.adj[subset,subset])/(n.genos^2)
    if(print==TRUE) print(c("Starting Value:",kin))
    while(i<(subset.genos)){
      if(!subset[j]%in%save){
        for(k in 1:extra.genos){
          temp.log<-sum(mat.adj[c(subset[-j],extra[k]),c(subset[-j],extra[k])])/(n.genos^2)
          if(temp.log<kin){
            i<-0
            kin<-temp.log
            old<-subset[j]
            subset<-c(subset[-j],extra[k])
            extra[k]<-old
            if(print==TRUE) print(c("Value:",kin))
          }
        }
      }
      if(j==subset.genos){
        j<-0
        if(print==TRUE) print("Cycle Completed, Returning to Beginning of Subset")
      }
      j<-j+1
      i<-i+1
    }
  }
  if(print==TRUE) print("Done")
  return(subset)
}
