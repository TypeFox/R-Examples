varRegressors <-
function(Data, VarName, RegressorNames=NULL)
{
#Codes out new variables/columns in dataframe to represent numeric regressors using pre-defined contrasts for factor VarName
#Used for control of a factor 'covariate' when graphing and ignoring this factor and/or other lower-level control
#data is the dataframe
#VarName is string name of factor
#RegressorNames is string vector of names for regressors.  If NULL, use contrast names for VarName
#2010-06-18:  added RegressorNames as argument, JJC
  
  
  
  Regressors = model.matrix( as.formula(paste('~', VarName)), data = Data)  
  
  dRegressors = data.frame(Regressors[,2:ncol(Regressors)])
  nRegressors = ncol(Regressors) - 1
  if ((!is.null(RegressorNames)) & (length(RegressorNames) == nRegressors))
  {
    colnames(dRegressors) = RegressorNames  
  }
  
  
  match = is.element(colnames(Data),colnames(dRegressors))  #find if regressors already in data frame
  Data = subset(Data, select = !match)   #remove them
  
  if (sum(match) > 0)
  {       
    warning(c('The following variables were deleted and then re-added to dataframe: ', paste(colnames(dRegressors), collapse= ', ')))
  } 
  
  dNew = dfMerge(Data,dRegressors)  #merge new versions of regressors with data
  
  return(dNew)
}

