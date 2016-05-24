#' Reads a data file in a monolix Project format for the simulator  Simulx using mlxtoobox code.
#'  
#' @param dataFile : the name (path) of the data file 
#' @param colTypes : vector of strings, the types of different column present in the data file,
#'                  a column of a type not considered in mlxtoolbox or whose values is not desired can be declared as "IGNORE"
#'   
#' @return  datas: a list of containg:
#' \itemize{
#'   \item \code{sources} :   contains the treatment informations,
#'   \item \code{observation} : contains informations on observations,
#'   \item \code{covariates} : contains the covariates parameters.
#' }
#' 
#' @export       

mlxtoolboxReadData <- function(dataFile,colTypes)
{ 
  myOldENVPATH = Sys.getenv('PATH');
  initMlxLibrary()
  session=Sys.getenv("session.simulx")
  Sys.setenv(LIXOFT_HOME=session)
  argList <- list(TXT_FILE=dataFile, COL_TYPES=colTypes)
  dot_call<- .Call
  datas2 <- dot_call( "mlxDataReaderR", argList, PACKAGE = "mlxDataReaderR");
  
  # set the storage of datas2 into the  format of datas  returned by readdatamlx.R
  obsi = 0
  covi = 0
  idsources = 0
  idcovariate = 0
  idobservation <- c()
  idcovariate <-c()
  for(i in 1:(length(datas2)-1))
  {
    if(datas2[[i]]$label == "dose")
    {
      idsources = i
    }
    if(datas2[[i]]$label == "longitudinal")
    {
      obsi = obsi+1
      idobservation<- c(idobservation,i)      
    }
    if(datas2[[i]]$label == "covariate")
    {
      covi = covi +1
      idcovariate <-c(idcovariate,i)
    }
    
  }
  
  if(idsources)
  {
    sources<-list(label="source",name="doseRegimen", colNames=datas2[[idsources]]$colTypes,
                  value=matrix(unlist(datas2[[idsources]]$values),nrow=length(datas2[[idsources]]$values),byrow = TRUE))  
  } else
  {
    sources =NULL
  }
  observation <-c()
  for(i in  1:obsi)
  {
    obsvalue=matrix(unlist(datas2[[idobservation[i]]]$values),nrow=length(datas2[[idobservation[i]]]$values),byrow = TRUE)
    observation<- c(observation,list(list( label="observation", colNames=c("id", "time"),
                                           value=obsvalue[,1:2])))
  }
  
  datas <-list(sources=sources,observation=observation)
  
  if(covi)
  {
    covariate <-c()
    for(i in  1:covi)
    {
      covariate<-c(covariate,list(list(name=datas2[[idcovariate[i]]]$colNames[2:length(datas2[[idcovariate[i]]]$colNames)],
                                       value=matrix(unlist(datas2[[idcovariate[i]]]$values),nrow=length(datas2[[idcovariate[i]]]$values),byrow = TRUE),
                                       label=datas2[[idcovariate[i]]]$label, colNames=c("id",datas2[[idcovariate[i]]]$colNames[2:length(datas2[[idcovariate[i]]]$colNames)]) )))
    }
    datas <- append(datas,list(covariate=covariate))
  }
  Sys.setenv(LIXOFT_HOME="")
  Sys.setenv('PATH'=myOldENVPATH);
  return(datas)
}
