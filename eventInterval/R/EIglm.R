EIglm<-function (event_times,event_vars=list(),formula=NULL) {
 event_times<-sort(event_times)
 event_ints<-get_intervals(event_times)
 event_times<-event_times[-1]
 if(is.null(formula))
 formula<-paste("event_ints",ifelse(length(event_vars), 
  paste(c("event_times",paste("event_vars",names(event_vars), 
   sep="$")),collapse="+"),"event_times"),sep = "~")
 cat("\nformula:", formula)
 formula<-eval(parse(text=formula))
 ei_glm<-glm(formula,family="Gamma")
 print(summary(ei_glm))
 invisible(ei_glm)
}
