rlsqr=function(G=matrix(),u=vector(),wantse=0,damp=0,
               atol=0,btol=0,conlim=0,itnlim=100,nout=0)
{
#m,n: m x n, dimension of G matrix
#damp : damping parameters
#wantse : A logical values to say if want to see standard error
#iw   : an integer array (NA + JA)
#       NA is the number of nonzero coefficients in each row of G
#       JA is the column in which ith nonzero of G lies
#rw   : a real array (non-zero values of sparse matrix)
#u    : input data corresponding d vector, (d = Gx)
#x    : computed solution x
#atol : an estimate of the relative error in the data defining the G matrix
#btol : an estimate of the relative error in the data defining the vector b
#conlim : A upper limit on condition number of A
#itnlim : iteration number
#nout : file number for printed output if positive.

  
m=dim(G)[1]
n=dim(G)[2]

M=Matrix(Matrix::t(G),sparse=TRUE)

rw=attributes(M)$x
iw=c(diff(attributes(M)$p+1),attributes(M)$i+1)

v=rep(0,n);w=rep(0,n);x=rep(0,n)
se=c(1:n)
istop=6;itn=0;anorm=0;acond=0;rrnorm=0;arnorm=0;xnorm=0
lx=.Fortran('LSQR',
      m=as.integer(m),
      n=as.integer(n),
      damp=as.double(damp),
      wantse=as.integer(wantse),
      leniw=as.integer(length(iw)),
      lenrw=as.integer(length(rw)),
      iw=as.integer(iw),
      rw=as.double(rw),
      u=as.double(u),
      v=as.double(v),
      w=as.double(w),
      x=as.double(x),
      se=as.double(se),
      atol=as.double(atol),
      btol=as.double(btol),
      conlim=as.double(conlim),
      itnlim=as.integer(itnlim),
      nout=as.integer(nout),
      istop=as.integer(istop),
      itn=as.integer(itn),
      anorm=as.double(anorm),
      acond=as.double(acond),
      rrnorm=as.double(rrnorm),
      arnorm=as.double(arnorm),
      xnorm=as.double(xnorm), package='PEIP')
return(list(x=lx$x,itn=lx$itn,xnorm=lx$xnorm,rnorm=lx$rrnorm, arnorm=lx$arnorm, acond=lx$acond, anorm=lx$anorm, istop=lx$istop ))
}





