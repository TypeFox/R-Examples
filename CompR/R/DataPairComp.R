setClass(
  Class="DataPairComp",
  representation=representation(
    Cons="character",
    Crit="character",
    Prod="character",
    Paircomp="list"),
  validity=function(object)
  {
    cat("------- DataPairComp  : validity --------\n")
    if (length(object@Crit)!=length(object@Paircomp))
    {
      print(length(object@Crit))
      print(length(object@Paircomp))
      stop("[DataPairComp : validity] Number of criteria is not equal to the length of Paircomp \n")
    }
    else
    {
      for (k in 1:length(object@Crit))
      {
        
        if (length(object@Paircomp[[k]])!=length(object@Cons))
        {
          
          stop("[DataPairComp : validity] Number of pairwise comparison matrices is not equal to the number of individuals\n")
        }
        else
        {
          
          for (h in 1:length(object@Cons))
          {
            
            if ((ncol(object@Paircomp[[k]][[h]])!=length(object@Prod))|(nrow(object@Paircomp[[k]][[h]])!=length(object@Prod)))
            {
              
              stop("[DataPairComp : validity] Dimension of pairwise comparison matrices is not equal to the number of products * number of products \n")
            }
            else
            {
              a<-as.vector(object@Paircomp[[k]][[h]])
              a<-(a<0)
              a<-as.logical(a)
              if (any(a==TRUE))
              {
                stop("[DataPairComp : validity] some values in pairwise comparison matrices are negative \n")
              }
            }
          }
        }
      }
      cat("------- DataPairComp  : validity    OK --------\n")
    }
    return(TRUE)
    
  }
)

setGeneric("getCons",
           function(object)
           {
             standardGeneric("getCons")
           }
           
)

setMethod("getCons","DataPairComp",
          function(object)
          {
            return(object@Cons)
          }
)

setGeneric("getCrit",
           function(object)
           {
             standardGeneric("getCrit")
           }
           
)

setMethod("getCrit","DataPairComp",
          function(object)
          {
            return(object@Crit)
          }
)

setGeneric("getProd",
           function(object)
           {
             standardGeneric("getProd")
           }
           
)

setMethod("getProd","DataPairComp",
          function(object)
          {
            return(object@Prod)
          }
)

setGeneric("getPaircomp",
           function(object)
           {
             standardGeneric("getPaircomp")
           }
           
)

setMethod("getPaircomp","DataPairComp",
          function(object)
          {
            return(object@Paircomp)
          }
)

setMethod("show","DataPairComp",
          function(object){
            cat("\n*** Class DatPairComp, method Show***\n")
            nrowShow<-min(10,nrow(object@Cons))
            ncolShow<-min(10,length(object@Prod))
            cat("\n* Cons (limited to 10 individuals) = ",formatC(object@Cons[1:nrowShow]),"\n")
            cat("\n*Crit =",object@Crit, "\n")
            cat("\n*Prod (limited to 10 products) = ",formatC(object@Prod[1:ncolShow]),"\n")
            cat("\n*PairComp (for all criteria, limited to 4 individuals and 10 products) = \n")
            for (i in 1:length(object@Crit))
            {
              cat("\n**********Crit = ",object@Crit[i],"****************\n")
              cat("\n")
              for (j in 1:4)
              {
                cat("\n",j,"\n")
                print(object@Paircomp[[i]][[j]][1:min(10,length(object@Prod)),1:min(10,length(object@Prod))])
              }
            }
          }
)