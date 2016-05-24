summary.ellipsefitlist <-
function(object,N=1000,boot=TRUE,seed=NULL,...) {
  g <- object
  summarycall <- match.call()
  if (!is.null(seed)) set.seed(seed)
  if (boot==TRUE) {
res <- mapply(summary.ellipsefit,g$models,MoreArgs=list(N=N,boot=TRUE,...),SIMPLIFY=FALSE)

values <- do.call("rbind",lapply(res,function (x)
  data.frame("Parameter"= rownames(x$values),x$values)))
rownames(values) <- NULL


if (length(dim(g$models))==1) {
thesubjects <- rep(rownames(g$Estimates),each=length(g$models[[1]]$values))
values <- cbind("Subject"=thesubjects,values)
Estimates <- matrix(values[,"Boot.Estimate"],nrow=length(g$Estimates[,1]),byrow=TRUE)
rownames(Estimates) <- rownames(g$Estimates)
colnames(Estimates) <- names(g$models[[1]]$values)
if (g$models[[1]]$method!="harmonic2" & g$models[[1]]$method!="geometric") Estimates <- Estimates[,c("b.x","b.y",
                          "cx","cy","retention","coercion","area",
                          "lag","split.angle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg","rote.rad",
                          "semi.major","semi.minor","focus.x","focus.y","eccentricity")]
else Estimates <- Estimates[,c("b.x","b.y","phase.angle",
                               "cx","cy","retention","coercion","area",
                               "lag","split.angle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg","rote.rad",
                               "semi.major","semi.minor","focus.x","focus.y","eccentricity")]
Std.Error <- matrix(values[,"Std.Error"],nrow=length(g$Estimates[,1]),byrow=TRUE)
rownames(Std.Error) <- rownames(g$Estimates)
colnames(Std.Error) <- names(g$models[[1]]$values)
if (g$models[[1]]$method!="harmonic2" & g$models[[1]]$method!="geometric") Std.Error <- Std.Error[,c("b.x","b.y",
                                    "cx","cy","retention","coercion","area",
                                    "lag","split.angle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg","rote.rad",
                                    "semi.major","semi.minor","focus.x","focus.y","eccentricity")]
else Std.Error <- Std.Error[,c("b.x","b.y","phase.angle",
                               "cx","cy","retention","coercion","area",
                               "lag","split.angle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg","rote.rad",
                               "semi.major","semi.minor","focus.x","focus.y","eccentricity")]}

else {
  length.values <- length(g$models[[1]]$values)
  subjectlen <- length(colnames(g$Estimates))-length.values
  subjectdata <- apply(g$Estimates[,1:subjectlen],2,rep,each=length.values)
  colnames(subjectdata) <- colnames(g$Estimates)[1:subjectlen]
  values <- cbind(subjectdata,values)
  Estimates <- matrix(values[,"Boot.Estimate"],nrow=length(g$Estimates[,1]),byrow=TRUE)
  Estimates <- cbind(g$Estimates[,1:subjectlen],Estimates) 
  colnames(Estimates) <- colnames(g$Estimates)
  Std.Error <- matrix(values[,"Std.Error"],nrow=length(g$Estimates[,1]),byrow=TRUE)
  Std.Error <- cbind(g$Estimates[,1:subjectlen],Std.Error) 
  colnames(Std.Error) <- colnames(g$Estimates)
  dim(res) <- dim(g$models)
  dimnames(res) <- dimnames(g$models)
}
res <- list("models"=res,"values"=values,"Boot.Estimates"=Estimates,"Boot.Std.Errors"=Std.Error,"Delta.Std.Errors"=g$Std.Errors)
res$summarycall <- summarycall  
class(res) <- "ellipsesummarylist"}
  else res <- g
  res
  }
