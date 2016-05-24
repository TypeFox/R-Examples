#' Study RT differences to lately remove outliers
#' 
#' Study the correlation between RT in positive and negative ionization modes to find those entities that have been associated wrongly.
#' @export
#' @param CommonEntitiesImproved The resultant data set from the function RemoveMismatch
#' @import plyr
#' @importFrom grDevices dev.off
#' @importFrom graphics hist
#' @importFrom stats cor
#' @importFrom utils read.table
#' @importFrom grDevices pdf
#' @importFrom graphics par
#' @importFrom stats lm
#' @importFrom utils write.table
#' @importFrom graphics plot
#' @importFrom stats na.omit
#' @importFrom stats residuals
#' @examples
#' \dontrun{
#' CommonEntitiesImproved<-StudyRTdiff(CommonEntitiesImproved)
#' }
#' @return Plot (RT+ vs RT-, regression, "residuals vs predicted", and Q-Q plot)
#' @return The CommonEntitiesImproved dataset now included a new column with residuals of each entity for the RT+ vs RT- regression.
StudyRTdiff<-function(CommonEntitiesImproved) {
  colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="RT+"] <- "x"
  colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="RT-"] <- "y"
  pdf("RTpositive_vs_RTnegative.pdf")
  plot(CommonEntitiesImproved[,'x'],CommonEntitiesImproved[,'y'],main="RT+ vs RT-")
  dev.off()
  LinearFit <- lm(y ~ x, data = CommonEntitiesImproved)
  par(mfrow=c(1,4))
  pdf("LinearRegressionCharacteristics.pdf")
  graph3<-plot(LinearFit)
  dev.off()
  pdf("Histogram_of_residuals.pdf")
  hist(residuals(LinearFit),breaks="FD")
  dev.off()
  CommonEntitiesImproved$Residuals<-residuals(LinearFit)
  CommonEntitiesImproved
}
