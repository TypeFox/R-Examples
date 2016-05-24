sortTree <-
function(tree,data,col.time,col.status,col.id,col.ctg,col.ctg.ordinal){
  leq<-rep(NA,NROW(tree))
  nodesTemp<-tree$node
  #Initialize for root node
  if(any(tree$var[1] %in% col.ctg) & !any(tree$var[1] %in% col.ctg.ordinal)){split<-(data[,tree$var[1]] %in% as.character(tree$cut[1]))
  } else {split<-(data[,tree$var[1]] <= as.numeric(as.character(tree$cut[1])))}  
  leq[1]<-suppressWarnings(coef(coxph(Surv(data[,col.time],data[,col.status])~split+cluster(data[,col.id])))<0)
  for(j in 2:NROW(tree)){
    if(!is.na(tree$score.test[j])){
      if(any(tree$var[j] %in% col.ctg) & !any(tree$var[j] %in% col.ctg.ordinal)){split<-(data[,tree$var[j]] %in% as.character(tree$cut[j]))
      } else {split<-(data[,tree$var[j]] <= as.numeric(as.character(tree$cut[j])))}
      parentNode<-rep(TRUE,NROW(data))
      for(i in 1:(nchar(tree$node[j])-1)){
        #Direction of split depends on next node (1 or 2)
        lastNode<-substr(tree$node[j],nchar(tree$node[j])-i+1,nchar(tree$node[j])-i+1)
        varTemp<-tree$var[substr(tree$node[j],1,nchar(tree$node[j])-i)==tree$node]
        cutTemp<-tree$cut[substr(tree$node[j],1,nchar(tree$node[j])-i)==tree$node]
        if(lastNode==1){
          if(any(varTemp %in% col.ctg) & !any(varTemp %in% col.ctg.ordinal)){
            parentNode[!(data[,varTemp] %in% as.character(cutTemp))]=FALSE
          } else {parentNode[data[,varTemp] > as.numeric(as.character(cutTemp))]=FALSE}
        } else {
          if(any(varTemp %in% col.ctg) & !any(varTemp %in% col.ctg.ordinal)){
            parentNode[data[,varTemp] %in% as.character(cutTemp)]=FALSE
          } else {parentNode[data[,varTemp] <= as.numeric(as.character(cutTemp))]=FALSE}
        }
      }
      #print(c(tree$size[j],sum(parentNode)))
      leq[j]<-suppressWarnings(coef(coxph(Surv(data[,col.time][parentNode],data[,col.status][parentNode])~split[parentNode]+cluster(data[,col.id][parentNode])))<0)
    }
    lastTempNode<-nodesTemp[substr(tree$node[j],1,nchar(tree$node[j])-1)==tree$node]
    endNode<-substr(nodesTemp[j],nchar(tree$node[j]),nchar(tree$node[j]))
    lastLEQ<-ifelse(leq[substr(tree$node[j],1,nchar(tree$node[j])-1)==tree$node],endNode,3-as.numeric(endNode))
    nodesTemp[j]<-paste0(lastTempNode,lastLEQ)
  }
  operator=ifelse(leq,"<=",">")
  operator[tree$var %in% col.ctg]=ifelse(leq[tree$var %in% col.ctg],"in","not in")
  tree$operator<-operator
  tree$node<-nodesTemp
  tree<-tree[order(tree$node),c(1,2,3,4,9,5,6,7,8)]
  return(tree)
}
