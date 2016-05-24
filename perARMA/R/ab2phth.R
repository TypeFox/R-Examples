#' @export
#' @importFrom matlab rem
#' @importFrom matlab isempty
#' @importFrom stats fft

ab2phth <-
function(a)
{ a=as.matrix(a)
  if (matlab::isempty(a))
    {phi=matrix()
     } else {
    T=nrow(a)
    T=as.numeric(T)
    p=ncol(a)
    p=as.numeric(p)

    A<-matrix(0,T,p)
    phi<-matrix(0,T,p)
    A[1,]=a[1,]
    for ( k in 1:(floor((T-1)/2)) )
        {A[(k+1),]=complex(real=a[2*k,]/2, imaginary= -a[(2*k+1),]/2)
        A[(T-k+1),]=Conj(A[(k+1),])}
    if (!matlab::rem(T,2))
        {A[((T/2)+1),]=a[T,]}

      for ( j in 1:p)
      {phi[,j]=Re(T*fft(A[,j],inverse=TRUE)/T)}

     }
    result = list(phi=phi)
    class(result) = "ab2phth"
    result
}

