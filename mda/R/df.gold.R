df.gold <-
  function(x,df.target,tol=0.001,trace.it=FALSE){
    dof=function(x,lambda)sum(x^2/(x^2+lambda))
    goldnew =function(x1,x3){
      phi=(1+sqrt(5))/2#golden ratio
      (x3+phi*x1)/(1+phi)
    }
    if(df.target<1)stop("Cannot have less than 1 as L1 df.target")
    lambda1=0
###find a lambda3 that undercuts
    lambda3=0.001
    df3=dof(x,lambda3)
    df1=dof(x,0)
    while(df3>df.target){
      lambda3=10*lambda3
      df3=dof(x,lambda3)
    }
    lambda2=goldnew(0,lambda3)
    df2=dof(x,lambda2)
    crit = abs(df2-df.target)
    while( crit>tol){
      if(trace.it) cat("df2",df2,"lambda",lambda2,"\n")
      lambda0=goldnew(lambda2,lambda3)
      df0=dof(x,lambda0)
      crit0=abs(df0-df.target)
      if(crit0<crit){
        lambda1=lambda2
        df1=df2
        lambda2=lambda0
        df2=df0
        crit=crit0
      }
      else {
        lambda3=lambda1
        df3=df1
        lambda1=lambda0
        df1=df0
      }
    }
    lambda2
  }

