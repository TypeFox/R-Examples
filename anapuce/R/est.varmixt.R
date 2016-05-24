est.varmixt <-
function(VAR,Kmax,dfreedom)
{ 
 nmixt<-1
 bic.old<-10000000
 arret=TRUE
 SSs <- sort(VAR)
  while (arret)
    {
      var.init <- vector(mode = "numeric", length = nmixt)
      var.init <- tapply(SSs, sort(rep(1:nmixt, length.out = length(VAR))),
                         mean)
      pi.init<-rep(1/nmixt,nmixt)
      res.mg<-mixgamma(VAR,dfreedom=dfreedom,var.init,pi.init,nmixt=nmixt,stop.crit=1e-4,
                        display = TRUE, niter.max = 50000,criterion="parameter")
      bic.new <- res.mg$BIC.crit
      cat("\n")
      if (bic.new-bic.old>0)
        {
          nmixt<-nmixt-1
          arret<-FALSE
        } else
        {
          bic.old<-bic.new
          var.old <- res.mg$vars               
          pi.old <- res.mg$p.i
          nmixt<-nmixt+1
          if (nmixt==Kmax)
            arret<-FALSE
        }
    }
  if (nmixt!=Kmax){
       res<-mixgamma(VAR,dfreedom=dfreedom,var.old,pi.old,nmixt=nmixt,stop.crit=1e-6,display=TRUE,
                  niter.max=50000,criterion="parameter")
                  }
  # (c) 2007 Institut National de la Recherche Agronomique
                
}

