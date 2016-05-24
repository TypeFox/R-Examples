autoprune <- function(formula, data, subset=1:length(data[,1]), ...){
  data = data[subset,]
  tree = rpart(formula, data=data, method="class",  
               cp=-1, minsplit=1, maxdepth=30,xval=max(10,length(subset)))
  minerr = which(tree$cptable[,4]==min(tree$cptable[,4]))
  xerr = sum(tree$cptable[minerr[1],4:5])	# min(xerror+1sd)
  cps = which(tree$cptable[,4]<= xerr)
  tree.prune = prune(tree, cp=tree$cptable[cps[1],1])
  class(tree.prune) <- "rpart"
  return(tree.prune)
}