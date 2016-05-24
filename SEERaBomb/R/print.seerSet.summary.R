print.seerSet.summary<-function(x, ...) {
#   cf=function (x) comma_format()(x)
  PY=year=NULL
  cat(x$title)
  print(x$d) #,row.names=FALSE)
  cat("Notes:\n")
  for (i in 1:length(x$notes))
    cat(i,") ",x$notes[i],"\n")
#   N=cf(sum(x$P$PY))
  N=sum(x$P$PY)
  tit=paste(N,"Million",ifelse(x$sex=="male","Male","Female"),"Person-Years")
  p=qplot(year,PY,ylab="Person-Years (Millions)",main=tit,data=x$P)+
    theme(axis.text=element_text(size=rel(2)),axis.title=element_text(size=rel(2)))
  print(p)
}
