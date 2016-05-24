diffareatriangles <-
function(area.triangle1,area.triangle2) {
  diff.areas=c()
  for (i in 1:length(area.triangle1))
  {
      diff.areas[i]=area.triangle1[i]-area.triangle2[i]
  }
  diff.auc=sum(diff.areas)
  answer=list(diffareas=diff.areas,diffauc=diff.auc)
  return(answer)
}
