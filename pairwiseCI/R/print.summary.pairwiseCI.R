`print.summary.pairwiseCI` <-
function(x,...)
{

pargs<-list(...)

conf.level<-attr(x,"conf.level")

methodI<-attr(x,"methodname")

  cat(" ","\n")
  cat(round(conf.level*100, 4), " %-confidence intervals", "\n")
  cat(" Method: ", methodI, "\n")
  cat(" ","\n")

print.default(x)

}

