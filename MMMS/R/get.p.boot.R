get.p.boot <-
function(chisq,chisq.boot) {
  mean(chisq<=extract.list(chisq.boot,"chisq.best"))
}
