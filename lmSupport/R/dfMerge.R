dfMerge <-
function(DataX, DataY, ByX=0, ByY=0, AllX=TRUE, AllY=TRUE, AddVars = TRUE)
#merges variables from two data frames (dX, dY) by default or merges cases across same dataframes if AddVars = FALSE
#If mergeing variables:
#By default matches on row names other variable name in d1 (by.x) and d2 (by.y)
#By default, includes all cases in dX and dY but can  limit to only matching (all.X=FALSE, all.y=FALSE) or left join (all.y=FALSE) or right join (all.x=FALSE) 
#
#If merging by cases (AddVars=FALSE)
#first adds missing vars to X and Y to have colnames match.  Missing vars set to NA
#Then merges cases
  {
  if (AddVars==TRUE)  #adding variables, matching on cases
  {
   dXY = merge(x=DataX,y=DataY,by.x=ByX, by.y=ByY, all.x=AllX, all.y =AllY, sort=FALSE)
   row.names(dXY) = dXY$Row.names
   dXY$Row.names = NULL
  }
  
  if (AddVars==FALSE)  #adding cases, matching on vars
  {
    
    #add vars from Y to X if needed and set to NA
    YMiss = !(colnames(DataY) %in% colnames(DataX))
    for (i in 1:ncol(DataY)){
      if (YMiss[i]) {
        print(sprintf("Adding %s and setting to NA for all cases in DataX", colnames(DataY)[i]))
        DataX[,colnames(DataY)[i]] = list(NA)
      }
    }
    
    #add vars from X to Y if needed and set to NA
    XMiss = !(colnames(DataX) %in% colnames(DataY))
    for (i in 1:ncol(DataX)){
      if (XMiss[i]) {
        print(sprintf("Adding %s and setting to NA for all cases in DataY ", colnames(DataX)[i]))
        DataY[,colnames(DataX)[i]] = list(NA)
      }
    }
    
    dXY = rbind(DataX,DataY)
  }    
   return(dXY)
}