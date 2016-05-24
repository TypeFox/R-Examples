# Independent rescaling of rows and columns as described in:
#
# Kluger, Y.; Basri, R.; Chang, J.T. & Gerstein, M., 
normalization=function(mat, type="irrc", error=0.000000001, maxit=1000)
  {
  if(type=="irrc")
    {
    irrc(mat)
    }
  else
    {
    if(type=="bistochastization")
      {
      bistochastization(mat, error, maxit)
      }
    else if(type=="logt")
      {
      logt(mat)
      }
    }
  }

# Spectral Biclustering of Microarray Data: Coclustering Genes and Conditions
# Genome Research 2003. 
#
# A - matrix to apply the rescaling
#
# Author: Rodrigo Santamaria (2007)
irrc=function(A)
  {
  R=apply(A,1,sum)^(-.5)
  C=apply(A,2,sum)^(-.5)
  An=t(t(R*A)*C)
  An
  }
  
  
# Bistocastization as described at:
#
# Kluger, Y.; Basri, R.; Chang, J.T. & Gerstein, M., 
# Spectral Biclustering of Microarray Data: Coclustering Genes and Conditions
# Genome Research 2003. 
#
# It is basically a rescaling until convergence (either maximum number of 
# iterations has been reached or change between two bistochastizations is
# lesser than a threshold)
#
# A- matrix to bistochastize
# error - minimum change between two iterations of bistochastization. Default
#          1e-9
# maxit - maximum number of bistochastizations. Default 10^3
#
# Author: Rodrigo Santamaria (2007)

bistochastization=function(A, error=0.000000001, maxit=1000)
  {
  dif=error+1
  numit=0
  An=A
  while(dif>error && numit<maxit)
    {
    Ap=An
    R=apply(A,1,sum)^(-.5)
    C=apply(A,2,sum)^(-.5)
    An=t(t(R*Ap)*C)
    dif=abs(sum(An-Ap))
    numit=numit+1
    }
  An
  }
  
# Logarithmic transform as described at:
#
# Kluger, Y.; Basri, R.; Chang, J.T. & Gerstein, M., 
# Spectral Biclustering of Microarray Data: Coclustering Genes and Conditions
# Genome Research 2003. 
#
# A - matrix to transform. Note that if A values are not greater than 1,
#     log transform will be erroneous.
#
# Author: Rodrigo Santamaria (2007)
logt=function(A)
  {
  n=dim(A)[1]
  m=dim(A)[2]
  L=log(A)
  Li=apply(L,1,sum)/m
  Lj=apply(L,2,sum)/n
  Lij=sum(L)/(n*m)
  K=t(t(L-Li+Lij)-Lj)
  K
  }
  
