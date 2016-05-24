
################################################################
#_____________    Function  Cooccurrences     __________________
#         Returns the occurrencesmatrix between stimuli
################################################################

Cooccurrences<-function(Part){
  
  if (!class(Part)=="SortingPartition"){
    return("The argument is not an object of class SortingPartition")
  }  
  else
  {
    nstim<-Part@nstimuli
    nsubjects<-Part@nsubjects
    Parti<-Part@Partition
    
    Cooccur<-matrix(0,nstim,nstim)
    for (subject in 1:nsubjects){
      
      PartSubject<-Parti[[subject]]
      CooccurSubject<-matrix(0,nstim,nstim)
      
      if (!is.list(PartSubject)){
        # unique sorting for the subject
        MatSubject<-PartSubject
        for (i in 1:(nstim-1)){
          for (j in (i+1):nstim){
            if (MatSubject[i]!=0 & MatSubject[j]!=0){
              if (MatSubject[i]==MatSubject[j]){
                CooccurSubject[i,j]<-1
                CooccurSubject[j,i]<-1
              }
            }    
          }
        }
        CooccurSubject<-CooccurSubject+diag(nstim)
      } else {
        # multiple sorting for the subject
        for (rep in (1:length(PartSubject))){
          MatSubject<-PartSubject[[rep]]
          for (i in 1:(nstim-1)){
            for (j in (i+1):nstim){
              if (MatSubject[i]!=0 & MatSubject[j]!=0){
                if (MatSubject[i]==MatSubject[j]){
                  CooccurSubject[i,j]<-CooccurSubject[i,j]+1
                  CooccurSubject[j,i]<-CooccurSubject[j,i]+1
                }
              }    
            }
          }
          CooccurSubject<-CooccurSubject+diag(nstim)
        }  
      }
      Cooccur<-Cooccur+CooccurSubject
    } #end for subject  
    colnames(Cooccur)<-Part@LabStim
    rownames(Cooccur)<-Part@LabStim
    return(Cooccur)
  }
} 
