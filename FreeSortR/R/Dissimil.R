################################################################
#__________________    Function  Dissimil     __________________
#               Generates a list of dissimilarities
################################################################

Dissimil<-function(Part){
  
  if (!class(Part)=="SortingPartition"){
    return("The argument is not an object of class SortingPartition")
  }  
  else
  {
    nstim<-Part@nstimuli
    nsubjects<-Part@nsubjects
    Parti<-Part@Partition
    
    ListDissim<-vector("list")
    for (subject in 1:nsubjects){
      
      PartSubject<-Parti[[subject]]
      DistSubject<-matrix(0,nstim,nstim)
      
      if (!is.list(PartSubject)){
        # unique sorting for the subject
        MatSubject<-PartSubject
        for (i in 1:(nstim-1)){
          for (j in (i+1):nstim){
            if (MatSubject[i]!=0 & MatSubject[j]!=0){
              if (MatSubject[i]!=MatSubject[j]){
                DistSubject[i,j]<-1
                DistSubject[j,i]<-1
              }
            }    
          }
        }
      } else {
        # multiple sorting for the subject
        for (rep in (1:length(PartSubject))){
          MatSubject<-PartSubject[[rep]]
          for (i in 1:(nstim-1)){
            for (j in (i+1):nstim){
              if (MatSubject[i]!=0 & MatSubject[j]!=0){
                if (MatSubject[i]!=MatSubject[j]){
                  DistSubject[i,j]<-DistSubject[i,j]+1
                  DistSubject[j,i]<-DistSubject[j,i]+1
                }
              }    
            }
          }   
        }  
      }
      ListDissim[[subject]]<-DistSubject
    } #end for subject  
    return(ListDissim)
  }
} 
