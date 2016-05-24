matrixplot <-
function(veg,rmember,use,y=1)
  {
  if(use == "columns") vvdat <- t(veg)
  if(use == "rows") vvdat <- veg
#
  rmember<- as.factor(rmember)
  lev<- levels(rmember)
#
  svdat<- cor(t(vvdat^y))
  svdat<- (svdat+1)*0.5
  svdat1<- rowsum(svdat, rmember,reorder=T)
  com<- t(rowsum(t(svdat1), rmember,reorder=T))
  f<- table(rmember)
  nf<- f%*%t(f)
  com<- com/nf
#
  mc<- length(com[,1])
  o.mxplot<- list(order=mc,mmatrix=com,levels=lev)
  }
