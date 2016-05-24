testNumberOfPools <- function(){
  t_start=0
  t_end=2
  tn=100
  tol=.02/tn
  print(tol)
  timestep=(t_end-t_start)/tn
  t=seq(t_start,t_end,timestep)
  nr=3
  nc=3
  # define the transfer functions for the model
  # we could compile them to a matrix valued
  # Function of C and t since they will be 
  # applied in a linear way on the output vector.
  # but we rather store them in an indexed list (
  # as a sparse matrix which also has some 
  # implementational benefits because the single
  # functions are easier to retrieve from the operator
  # if needed.
  alpha=list()
  alpha[["2_to_1"]]=function(C,t){
    1/4*1e-16
  }
  alpha[["3_to_2"]]=function(C,t){
    1/3*1e-16
  }
  alpha[["1_to_2"]]=function(C,t){
    1/2#*1e-16
  }

  k1=-1/2
  k2=-2/3
  k3=-1
  f=function(t,O){
    # in this case the Matrix N is actually completely linear
    N=matrix( 
       nrow=nr,
       ncol=nc,
       c(
          k1,    0,     0,  
          0  ,  k2,     0,  
          0,     0,    k3 
       )
    )
    # so we can write f(O)  as a Matrix product
    # note however that we could anything we like with the components
    # of O here. 
    # The only thing to take care of is that we release a vector of the same
    # size as O
    res=N*O
  }
  c01=3
  c02=2
  c03=2.5
  inputrates=new("TimeMap",t_start,t_end,function(t){return(matrix(
    nrow=3,
    ncol=1,
    c(
       1,  2,  3
    )
  ))})
  A=new("TransportDecompositionOperator",t_start,Inf,nr,alpha,f)
  mod=GeneralNlModel(
   t,
   A,
   c(
      c01,
      c02,
      c03
   ),
  inputrates,
  deSolve.lsoda.wrapper
  )
  # both 
  checkEquals(3,getNumberOfPools(mod))
  checkEquals(3,getNumberOfPools(A))
  

}
