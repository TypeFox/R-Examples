#' @title make setup_NUCOM
#' @description Code to make the setup structure needed run the model.
#'
#' @author JWM Pullens
#' @source The model can be sent upon request at jeroenpullens[at]gmail[dot]com
#'
#' @param mainDir Working directory
#' @param climate climate input (monthly) format: year, month, air temperature, precipitation, potential evapotranspiration (tab seperated). The potential evapotranspiration needs to be calcluated by using the Penman open water evapotranspiration.
#' @param environment environment input (yearly) format: year, atmospheric co2 values, nitrogen deposition
#' @param inival initial values of biomass
#' @param start year in which the simulation starts
#' @param end year in which the simulation ends
#' @param type Which output is needed? For more information see the help of the getData function.
#' @param numFolders The amount of folders that needs to be created (in case of parallel computing)
#' @param parallel Run the model on parallel cores? TRUE/FALSE, default is FALSE.
#' @param separate Does the model needs to be run for all parameters seperate? Default is FALSE
#' @param startval From which row does the output need to be loaded. Default is 1.
#'
#' @return A list with paths and filenames and parameter values which can be implemented in the runnucom and the runnucomParallel function.
#'
#' @keywords NUCOMBog
#'
#' @examples
#' \dontrun{
#' #Define complete file path in setup
#' for LINUX: ~/home/...../data/ ! pay attention to the last "/"
#' for Windows_ C://..//data// ! pay attention to the last "//"
#'
#' ##Single core setup:
#' test_setup_singlecore <- setupNUCOM(mainDir="/home/jeroen/NUCOMBog_data/",
#'                                      climate="ClimLVMhis.txt",
#'                                      environment="EnvLVMhis.txt",
#'                                      inival="inivalLVMhis.txt",
#'                                      start=1766,
#'                                      end=1999,
#'                                      type=c("NEE","WTD"),
#'                                      parallel=F)
#'
#' ## Multi core setup:
#' names<-c("CO2ref","gram_Beta","eric_MaxGr")
#'
#' nparvector<-50
#' initialParameters <- matrix(runif(n=length(names)*nparvector,
#'                    min=c(300,0.1,40),
#'                    max=c(500,1,80)),
#'                    nrow=length(names))
#' initialParameters<-data.frame(names,initialParameters)
#' names(initialParameters)<-c("names",rep("values",nparvector))
#' initialParameters$names<-as.character(initialParameters$names)
#'
#' test_setup <- setupNUCOM(mainDir="/home/jeroen/NUCOMBog_data/",
#'                           climate="ClimLVMhis.txt",
#'                           environment="EnvLVMhis.txt",
#'                           inival="inivalLVMhis.txt",
#'                           start=1766,
#'                           end=1999,
#'                           type=c("NEE","WTD"),
#'                           parallel=T,
#'                           numFolders=nparvector,
#'                           separate=F,
#'                           startval=1)
#'
#' }
#' @export

setupNUCOM<-function(mainDir,climate,environment,inival,start,end,type,numFolders=1,parallel=F,separate=F,startval=1){
  if(parallel==T){
    setup_parameters<-list()
    for(j in 1:numFolders){
      setup_parameters[[j]]<-list(mainDir=mainDir,runDir=paste(mainDir,"folder",j,"/",sep=""),climate=climate,environment=environment,inival=inival,start=start,end=end,type=type,parameters=NULL,separate=separate,startval=startval)
    }
    return(list(runParameters=setup_parameters))
  }
  if(parallel==F){
    setup_parameters<-list(list(mainDir=mainDir,runDir=mainDir,climate=climate,environment=environment,inival=inival,start=start,end=end,type=type,parameters=NULL,separate=separate,startval=startval))
  }
 }
