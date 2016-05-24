# multcompView : multcompLetters

pairwise.to.groups <-
function(pairwise.test,component="p.value",alpha=0.05) {
  mat <- pairwise.test[[component]]
  comp.mat <- combn(unique(c(colnames(mat),rownames(mat))),2)
  comp.vect <- apply(comp.mat,2,function(x) paste(x,collapse="-"))
  pval <- na.omit(as.vector(mat))
  names(pval) <- comp.vect
  result <- multcompView::multcompLetters(pval,threshold=alpha)
  return(result)
}
