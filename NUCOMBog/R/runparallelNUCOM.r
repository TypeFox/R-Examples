#' @title Run parallel NUCOM
#' @description Code to run NUCOMBog parallel on multiple cores.
#'
#' @author JWM Pullens
#' @source The model can be sent upon request at jeroenpullens[at]gmail[dot]com
#'
#' @param setup The setup needs to be made before by running the setup_NUCOM function.
#' @param clustertype Clustertype: The model has only been tested on SOCK cluster, which is the set to default.
#' @param numCores Number of Cores on which are model needs to be run (NOTE: Non-parallel runs can only be run on 1 core). Default is 1.
#' @param parameters The parameters which are used in the model. If no parameter values are given the default values will be used. The parameters have to have the format of a dataframe with colum names: "names" and "values", see example \url{https://github.com/jeroenpullens/NUCOMBog_data}. The default parameters are from Heijmans et al. 2008.
#' @keywords NUCOMBog
#'
#' @references Heijmans, M., Mauquoy, D., van Geel, B., and Berendse, F. (2008). Long-term effects of climate change on vegetation and carbon dynamics in peat bogs. Journal of Vegetation Science, 19(3)
#'
#' @examples
#' \dontrun{
#' !!the variable "test_setup" is from the function setupNUCOM, see the help for more information!!
#'
#' parallel<-runparallelNUCOM(setup = test_setup,
#'                             clustertype = "SOCK",
#'                             numCores = 1,
#'                             parameters=initialParameters)
#' }
#' @export
#' @import snowfall
#' @import utils

runparallelNUCOM<-function(setup,clustertype,numCores=1,parameters){
  if(clustertype=="MPI"){
    stop("MPI cluster type not supported.")
  }
  setwd(setup$runParameters[[1]]$mainDir)
  #we need to make the structure in all the folders
  runParameters<-setup$runParameters


  runParameters<-combine_setup_parameters(runParameters = runParameters,parameters = parameters)
  pb <- txtProgressBar(min = 0, max = length(runParameters), style = 3)
  print("Making Folder Structure")
  for (i in 1:length(runParameters)){
    setTxtProgressBar(pb, i)
    # copy files and folders:
    clim<-readLines(con=paste("input/",runParameters[[i]]$climate,sep=""))
    env<-readLines(con=paste("input/",runParameters[[i]]$environment,sep=""))
    ini<-readLines(con=paste("input/",runParameters[[i]]$inival,sep=""))

    filepath<-paste("folder",i,sep="")
    dir.create(filepath,showWarnings = F)

    # input folder with inival,clim,environm (these files do not change)
    dir.create(paste(filepath,"/input",sep=""),showWarnings = F)
    dir.create(paste(filepath,"/output",sep=""),showWarnings = F)
    writeLines(clim,paste(filepath,"/input/",runParameters[[i]]$climate,sep=""))
    writeLines(env,paste(filepath,"/input/",runParameters[[i]]$environment,sep=""))
    writeLines(ini,paste(filepath,"/input/",runParameters[[i]]$inival,sep=""))



    if(.Platform$OS.type=="unix"){
      file.copy(from = "modelMEE",to = paste(filepath,"/",sep="") )}
      if(.Platform$OS.type=="windows"){
        file.copy(from = "modelMEE.exe",to = paste(filepath,"/",sep="") )}

      }
      close(pb)
      print("Folder Structure Made")
      # then run nucom which creates the "param.txt","Filenames"

      # Create cluster
      print('Create cluster')
      if (clustertype =="SOCK"){
        snowfall::sfInit(parallel=TRUE, cpus=numCores, type="SOCK")
      # }else if (clustertype =="MPI"){
     # snowfall::sfInit(parallel=TRUE, cpus=numCores, type="MPI") # Dont need to specify type
      }

      # Exporting needed data and loading required packages on workers
      snowfall::sfLibrary("NUCOMBog",character.only = TRUE)
      snowfall::sfExport("setup","runParameters")

      # Distribute calculation: will return values as a list object
      cat ("Sending tasks to the cores\n")
      result =  snowfall::sfSapply(runParameters,runnucom_wrapper)

      # Destroy cluster
      snowfall::sfStop()
      return(result)
    }

