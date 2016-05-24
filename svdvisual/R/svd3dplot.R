svd3dplot <-
function(data, ncomp=3, irow=F, icol=F, isurface=T, iimage=F, xlab='Column', ylab='Row', zlab='', ...){
# this is the svd three dimensional plot


data.rank<-matrixrank(data);
if (ncomp>=data.rank){
	print('The specified rank is larger than the rank of the input matrix');
	print('The parameter will be reset as the rank of the input matrix -1');
	ncomp=data.rank-1;
}
data.dim<-dim(data);
nrow=data.dim[1];
ncol=data.dim[2];
ncell=nrow*ncol;

#make sure which mean will be removed.
meanmat<-NULL
localdata=data
idouble=(irow & icol);
meantitle=''
if (idouble){
  #calculate the double mean
  meanmat<-doublemean(data)
  localdata=data-meanmat;
  meantitle='Mean (Double)'
}

if (irow & (!idouble)){
  #calculate the row mean (mean of the rows)
  meanmat<-rowmean(data)
  localdata=data-meanmat;
  meantitle='Mean (Row)'
}
if (icol & (!idouble)){
  #calculate the column mean (mean of the columns)
  meanmat<-columnmean(data)
  localdata=data-meanmat;
  meantitle='Mean (Column)'
}

#if mean is chosen, reduce the number of component to ncomp-1
if (irow | icol) localncomp=ncomp-1;
if (!irow & !icol) localncomp=ncomp;

data.svd<-svd(localdata, nu=localncomp, nv=localncomp);
umat=data.svd$u;
vmat=data.svd$v;
svec=data.svd$d[1:localncomp];

localmeanmat=meanmat
if (is.null(meanmat)) localmeanmat=rep(0, nrow) %*% t(rep(0, ncol))

app=localmeanmat+umat%*%diag(svec)%*%t(vmat);
res=data-app;
rowmat=seq(1:nrow)%*% t(rep(1, ncol));
colmat=rep(1, nrow) %*% t(seq(1, ncol));

#generate ploting matrix
plotmat<-data.frame(matvalue=as.vector(data), rows=as.vector(rowmat), columns=as.vector(colmat), label=rep("Original data", ncell));
for (i in 1:ncomp){
  if (!irow & !icol){
  tempmat=svec[i]*(umat[, i]%*% t(vmat[, i]));
	plotmat<-rbind(plotmat, data.frame(matvalue=as.vector(tempmat), rows=as.vector(rowmat),  columns=as.vector(colmat), label=rep(paste("SVD", i), ncell)));}
  if (irow | icol){
    if (i==1){
      plotmat<-rbind(plotmat, data.frame(matvalue=as.vector(meanmat), rows=as.vector(rowmat), columns=as.vector(colmat), label=rep(meantitle)))
    }
    if (i>1){
      tempmat=svec[i-1]*(umat[, i-1]%*% t(vmat[, i-1]));
      plotmat<-rbind(plotmat, data.frame(matvalue=as.vector(tempmat), rows=as.vector(rowmat),  columns=as.vector(colmat), label=rep(paste("SVD", i-1), ncell)));      
    }
  }
}
	plotmat<-rbind(plotmat, data.frame(matvalue=as.vector(app), rows=as.vector(rowmat), columns=as.vector(colmat), label=rep("Approximation", ncell)));
	plotmat<-rbind(plotmat, data.frame(matvalue=as.vector(res), rows=as.vector(rowmat), columns=as.vector(colmat), label=rep("Residual", ncell)));

#require('lattice')

if (isurface==iimage){
	print('Can not simultaneously generate image plots and surface plots');
	print('Will only show the image plot');
	isurface=F;
	iimage=T;
}

nplot=ncomp+3;
nplotrow=ceiling(nplot/3);
localcondindex=c(1, nplot-1, nplot, 2:(nplot-2))

if (isurface){
	print(wireframe(matvalue~columns+rows|label, data=plotmat, xlab=xlab, ylab=ylab, zlab=zlab, main='SVD surface plot',  index.cond=list(localcondindex), layout=c(3, nplotrow), ...));
}

if (iimage){
	print(levelplot(matvalue~columns+rows|label, data=plotmat, xlab=xlab, ylab=ylab, zlab=zlab, main='SVD image plot', index.cond=list(localcondindex), layout=c(3, nplotrow), ...));
}
}
