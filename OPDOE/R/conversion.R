conversion.low_high.two_levels=function(M)
  { 
    A=M[,-c(1)]
    l=dim(A)[1]
    k=dim(A)[2]
    A=as.matrix(A)
    for (i in 1:l)
            {for (j in 1:k) 
                 {if (A[i,j]==0) A[i,j]=-1}
            }
    A=as.numeric(A)
    A=as.vector(A)
    A=matrix(A,l,k)
    dimnames(A)[[2]]=LETTERS[1:k]
  }

# converts the results of "conf.design" to designs 
# with usual "-1" and "+1" codings

