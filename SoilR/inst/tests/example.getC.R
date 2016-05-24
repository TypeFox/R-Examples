# This test function produces the example for testing the function getC
  t=seq(0,10,by=0.1) 
  k=0.8
  C0=100
  In = 30
  
  #Create two models with same arguments but calling two different model creation functions
  Cmodel=OnepModel(t,k,C0,In)
  C14model=OnepModel14(t,k,C0,In,F0=0,inputFc=IntCal09)
  
  #getC can extract the amount of C from these two type of models
  Ctmodel=getC(Cmodel)
  Ctmodel14=getC(C14model)
  
  #The output is identical because parameter values are the same
  plot(t,Ctmodel,type="l")
  lines(t,Ctmodel14,col=2,lty=2,lwd=2)
  legend("topright",c("OnepModel output","OnepModel14 output"),col=1:2,lty=1:2,bty="n")
