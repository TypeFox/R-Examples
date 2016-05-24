#'@title Function to retrieve data from the monthly output file created by NUCOMBog
#'
#'@author JWM Pullens
#'@source The model can be sent upon request at jeroenpullens[at]gmail[dot]com
#'
#'@description
#'This function returns the data from the monthly output file created by NUCOMBog.
#'
#'The original model provides net primary production (NPP) as an output, the model has been modified to provide autotrophic respiration aswell. In this way the net ecosystem exchange (NEE) can be calculated, since NEE = NPP - autotrophic respiration. The micrometeorological sign convention is used in this model, e.g. a negative value for NEE means carbon uptake. All fluxes are in gram carbon per square meter per month (gC m-2 month-1).
#'The model gives water table depth (WTD) in meters and positive values mean below ground level.
#'
#' The possible outputs of the model are Net Primary Production (NPP), Net Ecosystem Exchange (NEE), heterotrohpic respiration (hetero_resp) and water table depth (WTD). The desired output needs to be specified in the setup_NUCOM function.
#'
#' The getData function is integrated in all runnucom functions.
#'
#'@param setup setup_structure described in setup_NUCOM
#'@param startval When a spinup is used for the model and not all output is necessary, this "startval" parameter can be used to cut the output off, i.e. the starting row from which the "outmo.txt" file is loaded. Default is 1.
#'
#'@examples
#'\dontrun{
#'getData(setup=test_setup_singlecore,startval=1)
#'}
#' @export

getData<-function(setup,startval){
  out=list()
  NPP=numeric()
  NEE=numeric()
  WTD=numeric()
  hetero_resp=numeric()
  output<-read.csv(paste(setup$runDir,"/output/outmo.txt",sep=""),sep="",header=F,skip = 1,as.is=T)
  output<-output[startval:(nrow(output)-4),]
  outlist=data.frame(as.numeric(output$V1),as.numeric(output$V2))
  names(outlist)=c("year","month")

  if("NEE" %in% setup$type){
    for(i in 1:nrow(output)){
      NPP[i]<-(sum(output[i,4],output[i,8],output[i,12],output[i,16],output[i,20]))
      hetero_resp[i]<-sum(output[i,23:25])
      NEE[i]<- (-1*(NPP[i]-hetero_resp[i]))
    }
    outlist=cbind(outlist,NEE)
  }

  if("WTD" %in% setup$type){
    for(i in 1:nrow(output)){
      WTD[i]<-sum(-1*(output[i,26]/1000))
    }
    outlist=cbind(outlist,WTD)
  }

  if("hetero_resp" %in% setup$type){
    for(i in 1:nrow(output)){
      hetero_resp[i]<-sum(output[i,23:25])
    }
    outlist=cbind(outlist,hetero_resp)
  }

  if("NPP" %in% setup$type){
    for(i in 1:nrow(output)){
      NPP[i]<-(sum(output[i,4],output[i,8],output[i,12],output[i,16],output[i,20]))
    }
    outlist=cbind(outlist,NPP)
  }


  for (l in 1:ncol(outlist)){
  out[l]<-list(outlist[,l])
  }
  names(out) <- colnames(outlist)

  return(out)
}
