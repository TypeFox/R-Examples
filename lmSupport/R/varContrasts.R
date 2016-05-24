#Returns various types of contrast matrix for coding the categorical factor, TheFactor.
#Types inlcude DUMMY, EFFECTS, POC, HELMERT
#DUMMY, EFFECTS, and HELMERT allow you to specify a RefLevel.  
#Ref Leve indicates: For DUMMY- Control group, for EFFECTS- Excluded group, for HELMERT- 1=Reverse Helmert, otherwise, normal Helmert
#POC requires a list of Contrasts in POCList; e.g., list(c(1,0,-1), c(-1,2,-1)).   Best to provide as whole numbers.  Function will re-scale to unit contrast
#POC allows Labels.  If NULL (Default), contrast labels are POC1, POC2, etc.
#Last updated: 2010-03-05, JJC
#2010-10-11: changed parameter name of type to Type, JJC

varContrasts <- function(TheFactor, Type='DUMMY', RefLevel = length(levels(TheFactor)), POCList = NULL, Labels = NULL)

{
switch(toupper(Type),

   DUMMY =
   {
      ContrastMatrix = contr.treatment(levels(TheFactor), base = RefLevel)
      RefLevelName = levels(TheFactor)[RefLevel]
      
      #set contrast names
      for (i in 1:ncol(ContrastMatrix))
      {
         colnames(ContrastMatrix)[i] = paste(levels(TheFactor)[as.logical(ContrastMatrix[,i])],'_v_', RefLevelName, sep='')      
      }     
   },

   EFFECTS =
   {
      ContrastMatrix = contr.sum(levels(TheFactor))
      
      #change ref level
      if (RefLevel < length(levels(TheFactor)))
      {
        TempMatrix = ContrastMatrix
        TempMatrix[RefLevel,] = ContrastMatrix[length(levels(TheFactor)),]
        TempMatrix[length(levels(TheFactor)),] = ContrastMatrix[RefLevel,]
        ContrastMatrix = TempMatrix      
      }
      
      #set contrast names
      colnames(ContrastMatrix) = rep('colname', ncol(ContrastMatrix))  #change from NULL
      for (i in 1:ncol(ContrastMatrix))
      {   
         colnames(ContrastMatrix)[i] = paste(levels(TheFactor)[(ContrastMatrix[,i])>0],'_v_Mean', sep='')      
      }        
   },
   
   POC =
   {
      ContrastMatrix = matrix(data = NA, nrow = length(levels(TheFactor)), ncol = length(POCList))
      rownames(ContrastMatrix) = levels(TheFactor)
      for (i in 1:length(POCList))
      {
        ContrastMatrix[,i] = POCList[[i]]      
      }
      
      #Checks that coefficients are orthogonal
      if(max(abs(colSums(ContrastMatrix))) > 0)
      {
        print('ERROR: ContrastMatrix columns do not sum to ZERO')
      }
      
      for (i in 1:(ncol(ContrastMatrix)-1))
      {
      
        #NEED COLUMN PRODUCT TEST
      }
      
      #rescale contrast to unit difference 
      for (i in 1:ncol(ContrastMatrix))
      {
        ContrastMatrix[,i] = ContrastMatrix[,i] / (max(ContrastMatrix[,i]) - min(ContrastMatrix[,i]))
      }        
      
      #Add contrast labels
      if (is.null(Labels) | (length(Labels) != ncol(ContrastMatrix)))
      {
        colnames(ContrastMatrix) = paste(rep('POC',ncol(ContrastMatrix)),1:ncol(ContrastMatrix),sep='')
      } else
      {
        colnames(ContrastMatrix) = Labels      
      }      
   },
   
   HELMERT =
   {
      if (RefLevel ==1)  #REVERSE HELMERT
      {
        ContrastMatrix = contr.helmert(levels(TheFactor))   #this provides what I call reverse Helmert/SPSS Difference contrast
        colnames(ContrastMatrix) = rep('colname', ncol(ContrastMatrix))  #change from NULL
        for (i in ncol(ContrastMatrix):1)    #set contrast names
        {
          colnames(ContrastMatrix)[i] = paste(levels(TheFactor)[i+1],'_v_Earlier', sep='')          
        }
      } else   #for normal (forward) Helmert (DEFAULT)
      {
        ContrastMatrix = contr.helmert(levels(TheFactor))  #start with reverse Helmert
        TempMatrix = ContrastMatrix[nrow(ContrastMatrix):1,]  #reverse order contrasts
        rownames(TempMatrix) = rownames(ContrastMatrix)  #restore orignal row names
        ContrastMatrix = TempMatrix
        colnames(ContrastMatrix) = rep('colname', ncol(ContrastMatrix))  #change from NULL
        for (i in ncol(ContrastMatrix):1)   #set contrast names
        {
          colnames(ContrastMatrix)[i] = paste(levels(TheFactor)[nrow(ContrastMatrix)-i],'_v_Later', sep='')          
        }            
      }
      
      
      #rescale contrast to unit difference 
      for (i in 1:ncol(ContrastMatrix))
      {
        ContrastMatrix[,i] = ContrastMatrix[,i] / (max(ContrastMatrix[,i]) - min(ContrastMatrix[,i]))
      }  
      ContrastMatrix = ContrastMatrix[,ncol(ContrastMatrix):1]  #reverse order contrast columns
      
   },
   
   #OTHERWISE   
   {print('Valid options for type: DUMMY, EFFECTS, POC, and HELMERT')}
)#end switch

print(ContrastMatrix)
return(ContrastMatrix)

}

