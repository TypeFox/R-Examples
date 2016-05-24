store <-
function(experiment,outputEntireTable=T,path="./"){
  filename <- paste(path,experiment$Name,".txt",sep="")
  if (outputEntireTable==F){write.table(x=experiment$HMMtable,file=filename,quote=F,sep="\t",row.names=F)}
  else if (!is.null(experiment$genes)){write.table(experiment$genes,file=filename,quote=F,sep="\t",row.names=F)}
  cat("Results saved in file",filename,"\n")
  
}
