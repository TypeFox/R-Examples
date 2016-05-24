#' Plotting an illness-death-model.
#' 
#' Plotting an illness-death-model using \code{plot.Hist}.
#' 
#' 
#' @param stateLabels Labels for the three boxes.
#' @param style Either \code{1} or anything else, switches the orientation of
#' the graph. Hard to explain in words, see examples.
#' @param recovery Logical. If \code{TRUE} there will be an arrow from the
#' illness state to the initial state.
#' @param \dots Arguments passed to plot.Hist.
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{plotCompetingRiskModel}}, \code{\link{plot.Hist}}
#' @keywords survival
##' @examples
##' 
##' plotIllnessDeathModel()
##' plotIllnessDeathModel(style=2)
##' plotIllnessDeathModel(style=2,
##'                       stateLabels=c("a","b\nc","d"),
##'                       box1.col="yellow",
##'                       box2.col="green",
##'                       box3.col="red")
 
#' @export
plotIllnessDeathModel <- function(stateLabels,
                                  style=1,
                                  recovery=FALSE,
                                  ...){
  if (missing(stateLabels)) labels <- c("Disease\nfree","Illness","Death")
  if (recovery==TRUE){
    idHist <- Hist(time=1:4,event=list(from=c(1,1,2,2),to=c(2,3,1,3)))
    if (style==1)
      plot(idHist,
           stateLabels=stateLabels,
           box1.row=2,
           box1.column=1,
           box2.row=1,
           box2.column=3,
           ...)
    else{
      plot(idHist,
           stateLabels=stateLabels,
           ...)
    }
  }
  else{
    idHist <- Hist(time=1:3,event=list(from=c(1,1,2),to=c(2,3,3)))
    if (style==1){
      plot(idHist,
           stateLabels=stateLabels,
           box1.row=2,
           box1.column=1,
           box2.row=1,
           box2.column=3,
           ...)
    }
    else{
      plot(idHist,
           stateLabels=stateLabels,
           ...)
    }
  }
}
