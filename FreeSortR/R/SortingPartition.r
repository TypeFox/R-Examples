#####################################################################
#_________________    class  SortingPartition     ______________________
#               definition of the class SortingPartition
#####################################################################
setClass("SortingPartition",
         representation(type="character",
                        nstimuli="numeric",
                        nsubjects="numeric",
                        LabStim="vector",
                        LabSubj="vector",
                        Partition="list"),
         prototype     (type="Free",
                        nstimuli=1,
                        nsubjects=1,
                        LabStim="S",
                        LabSubj="X",
                        Partition=list(1)),
         validity=function(object){
           if (!is.character(object@LabStim))
             return(FALSE)
           else
             return(TRUE)        
         }
)

#####################################################################
#______________    Function  SortingPartition   _____________________
#           create an object of class SortingPartition
#      DonSort is a data frame with stimuli as rows
#                 and subjects as columns.
#        Rownames of DonSort are labels of stimuli 
#           and colnames are labels of subjects
#####################################################################

SortingPartition <- function(DataSort)
  
{ # test for multiple names of subjects
  # takes the part of the name of subject at the left of the dot
  cutdot<-function(text){unlist(strsplit(text,split="[.]"))[1]}
  cutnames<-unlist(lapply(colnames(DataSort),FUN=cutdot))
  
  if (length(colnames(DataSort))==length(unique(cutnames))){
    type<-"Free"
    # unique names of subjects
    nstimuli<-dim(DataSort)[1]
    nsubjects<-dim(DataSort)[2]
    Labels_stim<-rownames(DataSort)
    Labels_subj<-colnames(DataSort)  
  
    ListPart<-vector("list")
  
    for (subject in 1:nsubjects){
      ListPart[[subject]]<-DataSort[,subject]
    }
  } else {
    # Multiple names of subjects
      type="Multiple"
      cat("Multiple names observed in data\n")
      cat("Data will be considered as multiple sorting data\n")
      
      nstimuli<-dim(DataSort)[1]
      Labels_stim<-rownames(DataSort)
      
      namessubjects<-unique (cutnames)
      nsubjects<-length(namessubjects)
      Labels_subj<-cutnames  
      
      # each subject has a list of partitions 
      ListPart<-vector("list")
      
      for (subject in 1:nsubjects){
        #number of partitions for the subject
        num<-which(cutnames==namessubjects[subject])
        # case of one partition
        if (length(num)==1){
          ListPart[[subject]]<-DataSort[,num]
        } else {
          # case of multiple partitions
          DonSubject<-DataSort[,num]
          ListSubject<-vector("list")
          for (rep in (1:length(num))){
            ListSubject[[rep]]<-DonSubject[,rep]
          }
          ListPart[[subject]]<-ListSubject
        }
      }
  }
  return(new("SortingPartition",type=type,nstimuli=nstimuli,nsubjects=nsubjects,LabStim=Labels_stim,LabSubj=Labels_subj,Partition=ListPart))
  
}



#####################################################################
#_________________     method summary          ______________________
#           method summary for class SortingPartition
#####################################################################
setMethod(
  f ="summary",
  signature ="SortingPartition",
  definition = function(object){
    cat(object@type," sorting of ",object@nstimuli," stimuli by ",object@nsubjects," subjects.",sep="")
  }
)


#####################################################################
#_________________     method show          ______________________
#           method show for class SortingPartition
#####################################################################
setMethod(
  f ="show",
  signature ="SortingPartition",
  definition = function(object){
    TabPart<-NULL
    for (subject in (1:object@nsubjects)){
      # test for multiple partitions
      if (!is.list(object@Partition[[subject]])){
        TabPart<-cbind(TabPart,as.matrix(object@Partition[[subject]]))
      } else {
        for (rep in (1:length(object@Partition[[subject]]))){
          TabPart<-cbind(TabPart,as.matrix(object@Partition[[subject]][[rep]]))
        }
      }
    }
    
    rownames(TabPart)<-object@LabStim
    colnames(TabPart)<-object@LabSubj
    cat(object@type," sorting of ",object@nstimuli," stimuli by ",object@nsubjects," subjects.\n",sep="")
    cat("_________________________________\n\n")
    # printing the partitions or the 10 first partitions (if more of 10 subjects)
    if (ncol(TabPart)<15){
      print(TabPart)
    } else {
      cat("Partitions given by the first subjects : \n")
      print(TabPart[ ,1:15])
      cat("_______________________________\n")
    }
  }
)

#####################################################################
#_________________     method getPartition          ______________________
#           getter for Partition of class SortingPartition
#####################################################################
setGeneric("getPartition",
           function(object){standardGeneric("getPartition")}
)

setMethod("getPartition","SortingPartition",
          function(object){
            TabPart<-NULL
            for (subject in (1:object@nsubjects)){
              # test for multiple partitions
              if (!is.list(object@Partition[[subject]])){
                TabPart<-cbind(TabPart,as.matrix(object@Partition[[subject]]))
              } else {
                for (rep in (1:length(object@Partition[[subject]]))){
                  TabPart<-cbind(TabPart,as.matrix(object@Partition[[subject]][[rep]]))
                }
              }
            }
            rownames(TabPart)<-object@LabStim
            colnames(TabPart)<-object@LabSubj
            
            return(TabPart)
          }
)  



#####################################################################
#_________________     method nGroups          ______________________
#                returns the number of groups
#####################################################################
setGeneric("nGroups",
           function(object){standardGeneric("nGroups")}
)

setMethod("nGroups","SortingPartition",
          function(object){
            if (!class(object)=="SortingPartition"){
              return("This is not an object of class SortingPartition")
            }else{
              # function ng returns the number of groups given by a subject
              ng<-function(vec){return(sum(unique(vec)!=0))}
              
              if (object@type=="Free"){
                res<-simplify2array(lapply(object@Partition,FUN=ng))
              } else {
                # Multiple sorting
                res<-NULL
                for (subject in (1:object@nsubjects)){
                  # test for multiple partitions
                  if (!is.list(object@Partition[[subject]])){
                    res<-c(res,ng(object@Partition[[subject]]))
                  } else {
                    for (rep in (1:length(object@Partition[[subject]]))){
                      res<-c(res,ng(object@Partition[[subject]][[rep]]))
                    }
                  }  
                }
              }
              names(res)<-object@LabSubj 
              res
              return(res)   
            }  
          }
)  
