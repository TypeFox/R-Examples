################################################################################
plot.bifd<-function(x,argvals.s,argvals.t,...){
      if (missing(argvals.s)){
         nfine.s = max(c(201,10*x$sbasis$nbasis+1))
         argvals.s = seq(x$sbasis$rangeval[1],x$sbasis$rangeval[2],len=nfine.s)
         }
      if (missing(argvals.t)){
         nfine.t = max(c(201,10*x$tbasis$nbasis+1))
         argvals.t = seq(x$tbasis$rangeval[1],x$tbasis$rangeval[2],len=nfine.t)
         }
      tt<-list(argvals.s,argvals.t)
      rtt<-list(x$sbasis$rangeval,x$tbasis$rangeval)
      plot(fdata(eval.bifd(argvals.s,argvals.t,x),tt,rtt,fdata2d=TRUE),... )
      }
################################################################################