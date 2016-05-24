### ----------- dataprep1 ---------------------
dataprep1<-function(X){
  X<-as.matrix(X)
  if(length(colnames(X))==0){
    Iname<-nchar(paste(dim(X)[2]))
    colnames(X)<-paste("I",formatC(1:dim(X)[2], width = Iname, format = "d", flag = "0"),sep="") 
    # cat("no item names found in data" ,"\n", "items are named", colnames(X)[1], "(first item) to",  colnames(X)[dim(X)[2]],"(last item)","\n")
  }
  #if(length(rownames(X))==0){  # personen werden immer neu nummeriert
  Pname<-nchar(paste(dim(X)[1]))
  rownames(X)<-paste("P",formatC(1:dim(X)[1], width = Pname, format = "d", flag = "0"),sep="")  
  # cat("no person names (IDs) found in data" ,"\n", "persons are named", rownames(X)[1], "(first row) to",  rownames(X)[dim(X)[1]],"(last row)","\n")
  #}
  return(X)
}