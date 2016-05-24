#' @title Run NUCOMBog
#'
#' @description Code to run NUCOMBog on a single core.
#' @author JWM Pullens
#' @source The model can be sent upon request at jeroenpullens[at]gmail[dot]com
#'
#' @param setup The setup structure created by setup_NUCOM function needs to be inserted here, for more information see the setup_NUCOM function help, by typing "?NUCOMBog::setup_NUCOM".
#' @param parameters The parameters which are used in the model. If no parameter values are given the default values will be used. The parameters have to have the format of a dataframe with colum names: "names" and "values".
#' See example data available at \url{https://github.com/jeroenpullens/NUCOMBog_data}. The default parameters are from Heijmans et al. 2008.
#'
#' @references Heijmans, M., Mauquoy, D., van Geel, B., and Berendse, F. (2008). Long-term effects of climate change on vegetation and carbon dynamics in peat bogs. Journal of Vegetation Science, 19(3)
#'
#' @examples
#' \dontrun{
#' names<-c("CO2ref","gram_Beta","eric_MaxGr")
#' initialParameters <- c(380,0.5,65)
#' initialParameters<-data.frame(names,initialParameters)
#' names(initialParameters)<-c("names","values")
#'
#' runNUCOM(setup = test_setup_singlecore,parameters=initialParameters)
#'
#' ## with predefined parameters:
#' runnucom(setup = test_setup_singlecore,parameters=NULL)
#'}
#' @export

runNUCOM<-function(setup,parameters=NULL){

  startval=setup[[1]]$startval
   # print(parameters)
  if(!is.null(parameters)){
    setup<-combine_setup_parameters(runParameters = setup,parameters = parameters,parallel=F)
  }
  setup<-setup[[1]]

  # print(setup)
  setwd(setup$runDir)

  make_filenames(setup$runDir,setup$climate,setup$environment,setup$inival,setup$start,setup$end)

  make_param_file(setup$runDir,setup$parameters)

  if(.Platform$OS.type=="unix"){
  system("./modelMEE")}
  if(.Platform$OS.type=="windows"){
    system("./modelMEE.exe")}

  if(!is.null(setup$type)){
    out<-NUCOMBog::getData(setup,startval)
    return(out)
  }
}

