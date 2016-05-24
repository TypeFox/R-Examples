
#####################################################################
#_________________    class  SortingMds           ______________________
#               definition of the class SortingMds
#####################################################################
setClass("SortingMds",
         representation(nstimuli="numeric",
                        nsubjects="numeric",
                        LabStim="vector",
                        LabSubj="vector",
                        ndim="numeric",
                        Config="array",
                        Percent="vector",
                        Stress="numeric",
                        ResBoot="list"),
         prototype     (nstimuli=1,
                        nsubjects=1,
                        LabStim="S",
                        LabSubj="X",
                        ndim=1,
                        Config=as.array(0),
                        Percent=1,
                        Stress=0,
                        ResBoot=list(0)),
         validity=function(object){
           if (!is.character(object@LabStim))
             return(FALSE)
           else
             return(TRUE)        
         }
)

#####################################################################
#_________________     method summary          ______________________
#           method summary for class SortingMds
#####################################################################
setMethod(
  f ="summary",
  signature ="SortingMds",
  definition = function(object){
    cat("Mds results with ",object@ndim, " dimensions.\n",sep="")
  }
)


#####################################################################
#_________________     method show          ______________________
#           method show for class SortingMds
#####################################################################
setMethod(
  f ="show",
  signature ="SortingMds",
  definition = function(object){
      cat("Configuration of stimuli : \n")
      print(object@Config)
      cat("Stress = ",object@Stress)
      cat("_____________________")
    
  }
)


#####################################################################
#_________________     method getConfig          ______________________
#           getter for Config of class SortingMds
#####################################################################
setGeneric("getConfig",
  function(object){standardGeneric("getConfig")}
)

setMethod("getConfig","SortingMds",
  function(object){
    return(object@Config)
  }
)  
 
#####################################################################
#_________________     method getStress          ______________________
#           getter for Stress of class SortingMds
#####################################################################
setGeneric("getStress",
           function(object){standardGeneric("getStress")}
)

setMethod("getStress","SortingMds",
          function(object){
            return(object@Stress)
          }
)  

#####################################################################
#_________________     method getPercent          ______________________
#           getter for Percent of class SortingMds
#####################################################################
setGeneric("getPercent",
           function(object){standardGeneric("getPercent")}
)

setMethod("getPercent","SortingMds",
          function(object){
            return(object@Percent)
          }
)  