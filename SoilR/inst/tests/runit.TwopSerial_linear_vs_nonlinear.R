#
# vim:set ff=unix expandtab ts=2 sw=2:
# This test usese the non linear model approach for a linear problem to show that the results are consistent
test.TwopSerial_linear_vs_nonlinear=function(){
  require(RUnit)
  t_start=0
  t_end=20
  tn=100
  tol=.02/tn
  print(tol)
  timestep=(t_end-t_start)/tn
  t=seq(t_start,t_end,timestep)
  k1=1/2
  k2=1/3
  a21=1/9
  nr=2
  A=new("ConstLinDecompOp",
    matrix(
      byrow=TRUE,                                                           
      nrow=nr,
      ncol=nr,
      c(
        -k1,   0,  
        a21, -k2  
      )
    )
  )
  
  alpha=list()
  alpha[["1_to_2"]]=function(C,t){
    a21/k1
  }
  pe(quote(names(alpha)),environment())
  N=matrix( 
     nrow=nr,
     ncol=nr,
     c(
        k1,    0,  
        0  ,  k2  
     )
  )
  f=function(C,t){
    # in this case the application of f can be expressed by a matrix multiplication
    # f(C,t)=N C
    # furthermorde the matrix N is actually completely linear and even constant
    # so we can write f(C,t)  as a Matrix product
    # note however that we could anything we like with the components
    # of C here. 
    # The only thing to take care of is that we release a vector of the same
    # size as C
    return(N%*%C)
  }
  Anl=new("TransportDecompositionOperator",t_start,Inf,nr,alpha,f)
 
  
  c01=3
  c02=2
  iv=c(c01,c02)
  inputrates=new("TimeMap",t_start,t_end,function(t){return(matrix(
    nrow=nr,
    ncol=1,
    c( 2,  2)
  ))})
  #################################################################################
  # we check if we can reproduce the linear decomposition operator from the
  # nonlinear one
  Tr=getTransferMatrix(Anl) #this is a function of C and t
  T_00=Tr(matrix(nrow=nr,iv),0)
  af=getFunctionDefinition(A)
  af_0=af(0)
  pp("T_00",environment())
  pp("af_0",environment())
  pe(quote(T_00%*%N),environment())
  checkEquals(af_0,T_00%*%N)
  
          
  #################################################################################
  # build the two models (linear and nonlinear)
  mod=GeneralModel( t, A,iv, inputrates, deSolve.lsoda.wrapper) 
  modnl=GeneralNlModel( t, Anl, iv, inputrates, deSolve.lsoda.wrapper)
  # compare the Cstock
  Y=getC(mod) 
  Ynonlin=getC(modnl) 
  #R=getReleaseFlux(mod) 
  #Rnonlin=getReleaseFlux(modnl) 
#begin plots 
  lt1=2
  lt2=4
  m=matrix(c(1,1),1,1,byrow=TRUE)
  ex=expression(
    layout(m),
    plot(t,Y[,1],type="l",lty=lt1,col=1,ylab="Concentrations",xlab="Time",ylim=c(min(Y),max(Y))),
    lines(t,Ynonlin[,1],type="l",lty=lt2,col=1),
    lines(t,Y[,2],type="l",lty=lt1,col=2),
    lines(t,Ynonlin[,2],type="l",lty=lt2,col=2),
    legend(
    "topright",
      c(
      "linear sol for pool 1",
      "non linear sol for pool 1",
      "linear sol for pool 2",
      "non linear sol for pool 2"
      ),
      lty=c(lt1,lt2),
      col=c(1,1,2,2)
    )
  )  
  plotAndCheck("runit.ThreepSerial_linear_vs_nonlinear.pdf",ex,environment())
  #plot(t,R[,1],type="l",lty=lt1,col=1,ylab="Respirationfluxes",xlab="Time",ylim=c(min(R),max(R)))
  #lines(t,Rnonlin[,1],type="l",lty=lt2,col=1)
  #lines(t,R[,2],type="l",lty=lt1,col=2)
  #lines(t,Rnonlin[,2],type="l",lty=lt2,col=2)
  #lines(t,R[,3],type="l",lty=lt1,col=3)
  #lines(t,Rnonlin[,3],type="l",lty=lt2,col=3)
  #legend(
  #"topright",
  #  c(
  #  "linear sol for pool 1",
  #  "non linear sol for pool 1",
  #  "linear sol for pool 2",
  #  "non linear sol for pool 2",
  #  "linear sol for pool 3",
  #  "non linear sol for pool 3"
  #  ),
  #  lty=c(lt1,lt2),
  #  col=c(1,1,2,2,3,3)
  #)
# end plots 
# begin checks 
  checkEquals(
   Y,
   Ynonlin,
   "test non linear solution for C-Content computed by the ode mehtod against analytical",
   tolerance = tol,
  )
  #checkEquals(
  # R,
  # Rnonlin,
  # "test non linear solution for Respiration computed by the ode mehtod against analytical",
  # tolerance = tol,
  #)

 }
