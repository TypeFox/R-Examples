#' Filter by RT residuals
#' 
#' Remove those entities with residuals above and below a maximum and minimum specified value.
#' @param MaxResidual Maximum residual allowed for RT+ vs RT- association
#' @param MinResidual Minimum residual allowed for RT+ vs RT- association
#' @param CommonEntitiesImproved Data set resulted from the RemoveMismatch function
#' @examples
#' \dontrun{
#' CommonEntitiesFiltered<-FilterbyRT(CommonEntitiesImproved,MaxResidual=0.5,MinResidual=(-0.5))
#' }
#' @export
#' @return Plot filtered (RT+ vs RT-, regression, "residuals vs predicted", and Q-Q plot)
#' @return New CommonEntities table filtered, obtained after removing entities with very high or low residuals or RT+ vs RT-.
FilterbyRT<-function(CommonEntitiesImproved,MaxResidual,MinResidual) {
  CommonEntitiesPreFiltered<-CommonEntitiesImproved[(CommonEntitiesImproved$Residuals)<MaxResidual,]
  CommonEntitiesFiltered<-CommonEntitiesPreFiltered[(CommonEntitiesPreFiltered$Residuals)>MinResidual,]
  pdf("RTpositive_vs_RTnegative_Filtered.pdf")
  plot(CommonEntitiesFiltered[,'x'],CommonEntitiesFiltered[,'y'],main="RT+ vs RT-")
  dev.off()
  LinearFitFiltered <- lm(y ~ x, data = CommonEntitiesFiltered)
  par(mfrow=c(1,4))
  pdf("LinearRegressionCharacteristics_Filtered.pdf")
  graph3<-plot(LinearFitFiltered)
  dev.off()
  pdf("Histogram_of_residuals_Filtered.pdf")
  hist(residuals(LinearFitFiltered),breaks="FD")
  dev.off()
  colnames(CommonEntitiesFiltered)[colnames(CommonEntitiesFiltered)=="x"] <- "RT+"
  colnames(CommonEntitiesFiltered)[colnames(CommonEntitiesFiltered)=="y"] <- "RT-"
  write.table(CommonEntitiesFiltered,file="CommonEntitiesFiltered.csv",sep=",",row.names=FALSE,na="")
  CommonEntitiesFiltered
}
