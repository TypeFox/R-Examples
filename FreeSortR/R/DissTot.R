################################################################
#_________________    Function  DissTot    ____   ______________
#           Generates the overall matrix of dissimilarities
################################################################

DissTot<-function(Part){
  
  if (!class(Part)=="SortingPartition"){
    return("The argument is not an object of class SortingPartition")
  }  
  else
  {
    ListDissim<-Dissimil(Part)
    MatDissimTot<-apply(simplify2array(ListDissim),c(1,2),'sum') 
    rownames(MatDissimTot)<-Part@LabStim
    colnames(MatDissimTot)<-Part@LabStim
    
    return(MatDissimTot)
  }
} 
