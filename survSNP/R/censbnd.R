censbnd <-
function(lambda,p,crate,rootint=c(0.1,1000))
  {
    eqn<-function(M,lambda,p,crate)
      {
        lhs=p*((1-exp(-M*lambda))/lambda)
        sum(lhs)-M*crate
      }
    uniroot(eqn,rootint,lambda=lambda,p,crate=crate)
  }

