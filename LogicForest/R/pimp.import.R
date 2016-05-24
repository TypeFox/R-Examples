pimp.import <-
function(fit, data, testdata, BSpred, pred, Xs)#assumes have test data
{
  n<-nrow(testdata) #number obs in testdata
  tree<-fit$model$trees[[1]] #gives tree in form of class="logreg tree"
  y<-testdata[,(BSpred+1)]
  x.names<-colnames(data[,1:BSpred])
  orig.pred<-predict.logreg(fit, newbin=testdata[,1:BSpred]) #prediction on test data
  orig.miss<-sum(abs(orig.pred-y))/n #misclass rate on test data
  pimpinfo<-prime.imp(tree=tree, data=data, Xs=Xs)#determines prime imps for tree + tmp matrix
  vec.Xvars<-pimpinfo$vec.pimpvars#vector of Xvar ids from tree
  nxvars<-length(vec.Xvars)#number Xvars in tree
  single.vimp<-c()
  single.vimp.nms<-c()
  Xids<-c()
  for (i in 1:nxvars)#checking single var importance for tree
    {
    id<-vec.Xvars[i]
    Xid<-Xs[id]
    permute.ind<-sample(1:n, n, replace=FALSE)#indicators of premuted obs
    permute.col<-testdata[permute.ind, id] #permuted observations for given predictor
    pre.id<-if(id>1) as.matrix(testdata[,1:(id-1)]) 
    post.id<-if(id<pred) as.matrix(testdata[,(id+1):BSpred])
    permute.testdata<-cbind(pre.id, permute.col, post.id)
    perm.pred<-predict.logreg(fit, newbin=as.matrix(permute.testdata[,1:BSpred])) #pred perm
    perm.misclass<-sum(abs(perm.pred-y))/n   #number pred differ bt/ orig & permuted
    vimp<-perm.misclass-orig.miss
    single.vimp<-append(single.vimp, vimp)
    single.vimp.nms<-append(single.vimp.nms, x.names[id])
    Xids<-append(Xids, Xid)
    }
  names(single.vimp)<-single.vimp.nms
  pimpmat<-pimp.mat(pimps.out=pimpinfo, testdata=testdata)[[2]] #transforms data into predictors = pimps
  pimpnames<-pimp.mat(pimps.out=pimpinfo, testdata=testdata)[[1]] #vector of pimp names
  tmp.mat<-pimpinfo$tmp.mat
  zero.ids<-c()
  for(i in 1:ncol(tmp.mat))
    {
    ids<-if(all(tmp.mat[,i]==0)) {ids<-i}
    zero.ids<-append(zero.ids, ids)
    }
  if (length(zero.ids) > 0) {tmp.mat<-tmp.mat[,-zero.ids]}
  if (is.matrix(tmp.mat)) {npimps<-nrow(tmp.mat)}
  if (is.vector(tmp.mat)) {npimps<-1}
  pimp.vimp<-c()
  for (j in 1:npimps) #permuting pimps and checking misclass
    {
    perm.ind<-sample(1:n, n, replace=FALSE)#indicators of premuted obs
    perm.col<-pimpmat[perm.ind, j] #permuted observations for given predictor
    pre.j<-if(j>1) pimpmat[,1:(j-1)] 
    post.j<-if(j<npimps) pimpmat[,(j+1):npimps]
    permute.pimpdata<-cbind(pre.j, perm.col, post.j)
    pimp.pred<-c()
    for (k in 1:n)
      {
      pred<-ifelse(any(permute.pimpdata[k,]==1), 1, 0)
      pimp.pred<-append(pimp.pred, pred)
      }
    permpimp.miss<-sum(abs(pimp.pred-y))/n #permuted pimp misclass for jth pimp
    pvimp<-permpimp.miss-orig.miss #diff bt/ permutation and original
    pimp.vimp<-append(pimp.vimp, pvimp)
    }
  names(pimp.vimp)<-paste(pimpnames)
  out<-list(single.vimp=single.vimp, pimp.vimp=pimp.vimp, vec.Xvars=vec.Xvars, Xids=Xids)
}
