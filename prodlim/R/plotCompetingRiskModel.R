#' Plotting a competing-risk-model.
#' 
#' Plotting a competing-risk-model.
#' 
#' 
#' @param stateLabels Labels for the boxes.
#' @param horizontal The orientation of the plot.
#' @param \dots Arguments passed to \code{\link{plot.Hist}}.
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{plotIllnessDeathModel}}, \code{\link{plot.Hist}}
#' @keywords survival
#' @examples
#' 
#' plotCompetingRiskModel()
#' plotCompetingRiskModel(labels=c("a","b"))
#' plotCompetingRiskModel(labels=c("a","b","c"))
#'
#' @export
plotCompetingRiskModel <- function(stateLabels,horizontal=TRUE,...){
  if (missing(stateLabels)) stateLabels <- c("Disease\nfree","Cause1","Cause2")
  nTrans <- length(stateLabels)-1
  if (horizontal==TRUE){
    comprisk.model <- data.frame(time=1:3,status=1:3)
    CRHist <- with(comprisk.model,Hist(time,status,cens.code=2))
    plot(CRHist,stateLabels=stateLabels,...)
  }
  else{
    crHist <- Hist(time=1:nTrans,event=list(from=rep("1",nTrans),to=stateLabels[-1]))
    nrow <- 3
    if (nTrans/2==round(nTrans/2)){
      ncol <- nTrans+1
      midCol <- ceiling(ncol/2)
      columns <- c(midCol,(1:ncol)[-midCol])
      names(columns) <- paste("box",1:length(stateLabels),".column",sep="")
      rows <- c(1,rep(3,nTrans))
      names(rows) <- paste("box",1:length(stateLabels),".row",sep="")
    }
    else{
      ncol <- nTrans
      columns <- c(nTrans+1/2,1:nTrans)
      names(columns) <- paste("box",1:length(stateLabels),".column",sep="")
      rows <- c(1,rep(3,nTrans))
      names(rows) <- paste("box",2:length(stateLabels),".row",sep="")
    }
    do.call("plot.Hist",c(list(x=crHist,stateLabels=stateLabels,nrow=nrow,ncol=ncol,...),columns,rows))
  }
}
