################################################################
#__________    Function  DescriptionPartition     _____________
#               Description of a partition
################################################################

DescriptionPartition<-function(Part,subject=1,replicate=1,Labels=NULL){
  
  
  
  
  if (!class(Part)=="SortingPartition"){
    if (is.vector(Part)){
      return(Description(Part,Labels))
    } else {
      return("This is not a partition vector or an object of class SortingPartition")
    }    
  } else { 
    nstim<-Part@nstimuli
    nsubjects<-Part@nsubjects
    Labels<-Part@LabStim
    Labsubject<-unique(Part@LabSubj)
    
    # Search for the subject to describe
    if(is.character(subject)){
      num<-which(Labsubject==subject)
      if(sum(num)==0){
        cat("Subject ",subject," does not exist.\n","Description of the subject 1.\n",sep="")
        num<-1
      }
    } else {
      if((subject<=0) | (subject>nsubjects)){
        cat("Subject ",subject," does not exist.\n","Description of the subject 1.\n",sep="")
        num<-1
      } else {
        num<-subject
      }
    }
    
    # Search for the partition to describe
    Parti<-Part@Partition[[num]]
    if (is.list(Parti)){
      # Unique partition for this subject
      if (length(Parti)>=replicate){
        Parti<-Parti[[replicate]]
      } else {
        # Multiple subject
        cat("Replicate ",replicate," doesn't exist for subject ",subject," !\n",sep="")
        cat("Description of the partition 1 of the subject ",subject," : \n",sep="")
        cat("----------------------\n")
        Parti<-Parti[[1]]
      }  
    } else {
      if (replicate>1){
        cat("Replicate ",replicate," doesn't exist for subject ",subject," !\n",sep="")
        cat("Description of the partition 1 of the subject ",subject," : \n",sep="")
        cat("----------------------\n")
      }
    }
    Labels=Labels[Parti!=0]  
    Parti=Parti[Parti!=0]
    res<-Description(Parti,Labels=Labels)
    return(res)
  }
}
