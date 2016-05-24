##' @export
`PrbMlt.default` <- function(MultiData){
  LME <- try(lme(value~factor(Contrst), random=~1|variable, na.action=na.omit,data=MultiData))
  if(class(LME)!="try-error"){suppressMessages(suppressWarnings(anova(LME)[[4]][2] ))}else{1}  
}