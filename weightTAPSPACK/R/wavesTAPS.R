#' Find population statistics for different waves
#'
#' Find the population statistics for each wave in TAPS data to be analyzed. Useful for attrition information.
#'
#' @param outcome A character vector of the names of outcome variables of interest
#'
#' @return A list of population stats
#' @docType methods
#' @author David G. Carlson \email{carlson.david@@wustl.edu}  Michelle Torres: \email{smtorres@@wustl} Taeyong Park \email{t.park@@wustl.edu}
#' @seealso \code{\link{weightTAPSPACK}} \code{\link{weightTAPS}} \code{\link{variablesTAPS}} \code{\link{weightTAPSoutput}} \code{\link{simpleWeight}} \code{\link{attritTAPS}} \code{\link{multipleImp}} \code{\link{hotdeckImp}} \code{\link{subsetTAPS}}
#' @rdname wavesTAPS
#' @aliases wavesTAPS,ANY-method
#' @export
setGeneric(name="wavesTAPS",
           def=function(outcome)
           {standardGeneric("wavesTAPS")}
)

setMethod(f="wavesTAPS",
          definition=function(outcome){
    dataSets<-alply(outcome, 1, subsetTAPS)
    data.subset<-subsetTAPS(outcome=outcome)
 
  if (length(outcome)>1){
    data.waves<-Reduce(function(x, y) merge(x, y, all=TRUE), dataSets)
  }else{
    data.waves<-dataSets
  }
  
  names.waves<-character(0) #Waves names for legend (see loop)
  t.agegend<-t.ethm<-t.educat<-t.regmetro<-t.incomcat<- t.ppnet<-list()
  for (i in 1:length(outcome)){
  names.waves[[i]]<-paste("Outcome",i, sep="")
  t.agegend[[i]]<-round((prop.table(table(dataSets[[i]]$agegend)))*100,2)
  t.ethm[[i]]<-round((prop.table(table(dataSets[[i]]$ethm)))*100,2)
  t.educat[[i]]<-round((prop.table(table(dataSets[[i]]$educat)))*100,2)
  t.regmetro[[i]]<-round((prop.table(table(dataSets[[i]]$regmetro)))*100,2)
  t.incomcat[[i]]<-round((prop.table(table(dataSets[[i]]$incomcat)))*100,2)
  t.ppnet[[i]]<-round((prop.table(table(dataSets[[i]]$ppnet)))*100,2) 
  }
  t.total<-list(t.agegend, t.ethm, t.educat, t.regmetro, t.incomcat, t.ppnet)
  statsWave<-llply(t.total, .fun=function(x) do.call(cbind, x))
  for (i in 1:length(statsWave)){
    colnames(statsWave[[i]])<-names.waves 
  }

  return(statsWave)
  }
)
