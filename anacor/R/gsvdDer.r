gsvdDer <- function(tab, ndim)
{
# computes first order derivatives

#------------------------- data preparation --------------------------
  data <- tab
  nr<-nrow(data)
  nc<-ncol(data)
  nn<-nr*nc
  mOff<-array(0,c(nr,nc,nn))
  mRow<-array(0,c(nr,nr,nn))
  mCol<-array(0,c(nc,nc,nn))
  freq<-rep(0,nn)
  k<-1
  for (i in 1:nr) for (j in 1:nc) {
  	mOff[i,j,k]<-1
  	mRow[i,i,k]<-1
  	mCol[j,j,k]<-1
  	freq[k]<-data[i,j]                  #frequency vector (row-wise)
  	k<-k+1
  }
  N<-sum(freq)
  par <-freq/N
#----------------------------- end data preparation ---------------------

  p <- matrix(0,nr,nr)
	for (k in 1:nn) p <- p + par[k]*mRow[,,k]

  q <- matrix(0,nc,nc)
	for (k in 1:nn) q <- q + par[k]*mCol[,,k]
  
  r <- matrix(0,nr,nc)
	for (k in 1:nn) r <- r + par[k]*mOff[,,k]

  dr <- mOff
  dp <- mRow
  dq <- mCol

	nrows<-nrow(p)
  ncols<-nrow(q)
  npars<-length(par)

  #--------- call gsvd -----------
  gs <- gsvd(r,p,q)
  gu <- gs$gu
  gv<-gs$gv
  gd<-gs$gd
  ind <- 2:(ndim+1)          #all dimensions except trivial (to be selected in plot)
  neval<-length(ind)

  dl<-matrix(0,neval,npars)
	dx<-array(0,c(nrows,neval,npars))
  dz<-array(0,c(ncols,neval,npars))
	for (i in 1:neval) {
		j<-ind[i]
    x<-gu[,j]
    z<-gv[,j]
    gi<-gd[j]
		xz<-rbind(cbind(x),cbind(z))
		dpx<-apply(dp,3,function(d) x%*%d%*%x)
		dqz<-apply(dq,3,function(d) z%*%d%*%z)
		dl[i,]<-apply(dr,3,function(d) x%*%d%*%z)-gi*(dpx+dqz)/2

    #---------- call ginv --------------
    aux0 <- ginvgsvd(gs,p,j)
		kv<-array(0,c(nrows+ncols,nrows+ncols,npars))
		kv[1:nrows,1:nrows,]<--gi*dp
		kv[1:nrows,nrows+(1:ncols),]<-dr
		kv[nrows+(1:ncols),1:nrows,]<-aperm(dr,c(2,1,3))
		kv[nrows+(1:ncols),nrows+(1:ncols),]<--gi*dq
		aux1<-apply(kv,3,function(d) aux0%*%d%*%xz)
		aux2<-drop(outer(xz,(dpx+dqz)/4))
		dxz<--(aux1+aux2)
		dx[,i,]<-dxz[1:nrows,]
		dz[,i,]<-dxz[nrows+(1:ncols),]
		}
	return(list(gd=gd,gu=gu,gv=gv,dl=dl,dx=dx,dz=dz,ind=ind))
}