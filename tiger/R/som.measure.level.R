som.measure.level <- function(result, show.measures = 1:num.measures, mfrow=c(4,3)){
#    require(klaR)
    num.measures <- c(1:NCOL(result$som$data))
    newrange<-apply(result$som$code, 2, range)
    newrange
    centers.rer <- result$best.value.location$central.best.reranged
    center.label <- centers.rer
    center.label[is.na(centers.rer)]<-((newrange[1,]+newrange[2,])/2)[is.na(centers.rer)]
    theLabels <- signif((cbind(newrange[1,],center.label,newrange[2,])), digits=3)

    revert <- rep(TRUE, NCOL(result$measures))
    revert[names(result$measures) %in% c("CE", "PI",
    "Rsqr", "IoAd",  "lcs_slope", "GRI")] <- FALSE

level_shardsplot(result$som,
  rows = show.measures,
  mfrow=mfrow, par.names=result$names,
  centers=centers.rer, 
  class.labels=theLabels,
  revert.colors= revert, classcolors="gray")


}
