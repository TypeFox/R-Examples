#'Edit and Run an APSIM Simulation
#'
#'This function will edit then run an APSIM simulation.  It can be thought of as a vectorized
#'version of \code{\link[apsimr:apsim]{apsim}}.
#'
#'In order to link between the \code{\link[sensitivity]{sensitivity}} package and 
#'APSIM, there needs to be a function that can edit then run APSIM that produces
#'a univariate output.  The \code{apsim_vector} function satisfies these conditions
#'and can return multivariate output if a multivariate sensitivity analysis is of interest.
#'See the examples section or the package vignette for details.
#'
#' @param X N-by-p matrix of inputs with rows corresponding to runs and columns for variables
#' @param exe where to find the APSIM executable
#' @param wd directory where the .apsim file lives and where the results will be saved
#' @param vars names of the variables, must be of the same length as 'X' has columns
#' @param to.run the .apsim file in \code{wd} to run
#' @param to.edit the .apsim file or .xml file to be edited
#' @param overwrite logical; passed to the \code{edit_apsim} and \code{edit_sim_file} functions
#' @param g a function of the output returned by apsim - must give vector or scalar result
#' @return a vector of length N
#' @export
#' @examples
#' \dontrun{
#' meanCowpea<-function(X){
#'  return(mean(X$lai_cowpea))
#' }
#' 
#' apsimWd <- "~/APSIM"
#' apsimVar <- c(rep("SoilWater/Thickness",11), "SoilOrganicMatter/SoilCN")
#' apsimValue <- matrix(c(rep(200, 2), rep(300, 9), 10,
#'                        rep(350, 2), rep(350, 9), 5),nrow=2,byrow=T)
#' apsimExe <-"C:/Program Files (x86)/Apsim75-r3008/Model/Apsim.exe"
#' apsimFile <- "Canopy.apsim"
#' 
#' #Run APSIM at each new parameter vector specified by apsimVar and only return the mean cowpea value
#' uniRes <- apsim_vector(X = apsimValue, exe = apsimExe, wd = apsimWd, vars = apsimVar, 
#'          to.run = apsimFile, to.edit = apsimFile, g = meanCowpea)
#' uniRes
#' 
#' 
#' allCowpea <- function(x){
#'  return(x$lai_cowpea)
#' }
#' 
#' #Use allCowpea function to return all the cowpea values
#' multiRes <- apsim_vector(X = apsimValue, exe = apsimExe, wd = apsimWd, vars = apsimVar, 
#'              to.run = apsimFile, to.edit = apsimFile, g = allCowpea)
#' multiRes
#' }

apsim_vector<-function(X, exe, wd, vars, to.run, to.edit=to.run, overwrite=FALSE, g){
  #This is a version of the 'apsim' function that can be 
  #used with the sensitivity package
  oldWd<-getwd()
  setwd(wd)
  
  if(ncol(X)!=length(vars)){
    stop("X must have same number of columns as 'var' is long")
  }
  
  N <- nrow(X)  
  
  #Some variables have multiple values that need to be changed, for example
  #soil properties often need to be specified at multiple depths
  unVar<-unique(vars)
  numVar <- length(unVar)
  valueI<-vector("list",numVar)
  
  #When overwrite==FALSE then an new .apsim file is created with the title "___-edited.apsim"
  #Therefore, if to.run==to.edit and overwirte==FALSE then the unchanged version of to.edit is executed
  #This updates the file name so the correct (new) .apsim file is run
  if(!overwrite & to.edit==to.run){
    newName<-paste(gsub(".apsim","",to.run),"-edited",sep="")
    to.run<-paste(newName,".apsim",sep="")
  }
  
  for(i in 1:N){
    
    for(j in 1:numVar){
      valueI[[j]] <- X[i,which(vars==unVar[j])]
    }
    
    if(length(grep(".apsim$",to.edit))>0){
      #edit the .apsim file
      edit_apsim(file=to.edit,wd=wd,var=unVar,value=valueI,overwrite=overwrite)
    }else{
      edit_sim_file(file=to.edit,wd=wd,var=unVar,value=valueI,overwrite=overwrite)
    }

    #Run the edited apsim file and collect the output

    res <- apsim(exe=exe, wd=wd, files = to.run)
    if(i==1){
      g1 <- g(res)
      
      if(any(is.na(g1))){
        setwd(oldWd)
        stop("The provided g function returned NA, check that it is defined correctly.")
      }
      
      if(length(g1)>1){
        multi <- TRUE
        y <- matrix(0,N,length(g1))
        y[1,] <- g1
      }else{
        multi <- FALSE
        y <- rep(0,N)
        y[1] <- g1
      }
      
    }else{ 
      if(multi){
        y[i,] <- g(res)
      }else{
        y[i] <- g(res)
      }
    }
  }
  setwd(oldWd)
  return(y)
}