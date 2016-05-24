
# take each variable, create the appropriate sized basis
# then nest around matrices with elements 0
# mainly for getting compoenent wise std.Errs

create_component_design<-function(varb, X.list, j, k, sm=T){
  # basis range
  new.varb<-seq(min(varb),  max(varb), length.out = 100)
  # get column dimensions of other components
  l.covs<-lapply(X.list, ncol)
  component.size<-l.covs[[j]]
  if(sm) new.basis<-make_spam(bbase(new.varb, nseg = (component.size-3)))
  if(!sm) new.basis<-as.matrix(new.varb)
  empty.X.list<-lapply(l.covs, make_sparse, nrow = 100)
  empty.X.list[[j]]<-new.basis
  Reduce("cbind", empty.X.list)
}