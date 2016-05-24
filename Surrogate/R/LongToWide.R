LongToWide <- function(Dataset, OutcomeIndicator, IdIndicator,  TreatIndicator, OutcomeValue){

  OutcomeIndicator <- Dataset[,paste(substitute(OutcomeIndicator))]
  IdIndicator <- Dataset[,paste(substitute(IdIndicator))]
  TreatIndicator <- Dataset[,paste(substitute(TreatIndicator))]
  OutcomeValue <- Dataset[,paste(substitute(OutcomeValue))]
  Dataset <- data.frame(cbind(OutcomeIndicator, IdIndicator, TreatIndicator, OutcomeValue))
  
  reshape(data=Dataset, direction="wide", timevar="OutcomeIndicator", idvar=c("IdIndicator", "TreatIndicator"))

}