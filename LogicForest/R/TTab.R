TTab <-
function(data, tree, Xs)
{
  if(!is(tree,"logregtree"))
  stop("tree must be an object of class logregtree")
  mod.var<-tree$trees[,3]
  mod.var<-sort(mod.var[mod.var!=0])
  mod.var<-mod.var[!duplicated(mod.var)]
  model.var<-c()
  for (i in 1:length(mod.var))
    {
    modvar<-Xs[mod.var[i]]
    model.var<-append(model.var, modvar)
    }
  mat.perms<-Perms(length(mod.var))
  nms<-colnames(data)
  if (is.null(colnames(data))) {colnames(mat.perms)<-paste("X", model.var, sep="")}
  if (length(nms)>0) {colnames(mat.perms)<-nms[mod.var]}
  mat.bin<-matrix(0, nrow(mat.perms), max(mod.var))
  mat.bin[,mod.var]<-mat.perms
  pred.out<-eval.logreg(tree, mat.bin)
  mat.truth<-cbind(mat.perms, outcome=pred.out)
  truth<-ifelse(tree$coef>0 | is.na(tree$coef),1,0)
  ids.truth<-mat.truth[,"outcome"]==truth
  mat.truth<-mat.truth[ids.truth,-ncol(mat.truth),drop=FALSE]
  mat.truth
}
