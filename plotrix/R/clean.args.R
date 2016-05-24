clean.args<-function(argstr,fn,exclude.repeats=FALSE,exclude.other=NULL,
 dots.ok=TRUE) {

 fnargs<-names(formals(fn))
 if(length(argstr) > 0 && !("..." %in% fnargs && dots.ok)) {
  badargs<-names(argstr)[!sapply(names(argstr),"%in%",c(fnargs,""))]
  for(i in badargs) argstr[[i]]<-NULL
 }
 if(exclude.repeats) {
  ntab<-table(names(argstr))
  badargs<-names(ntab)[ntab > 1 & names(ntab) != ""]
  for (i in badargs) argstr[[i]]<-NULL
 }
 for(i in exclude.other) argstr[[i]]<-NULL
 argstr
}

remove.args<-function(argstr,fn) {
 fnargs <- names(formals(fn))
 argstr[!(names(argstr) %in% fnargs)]
}
