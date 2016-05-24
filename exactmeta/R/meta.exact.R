meta.exact <-
function(data, type="risk difference", BB.grdnum=2000, B.sim=20000, cov.prob=0.95, midp=T, print=T, studyCI=T, ratio.upper=1000)
          {if(type=="risk difference")
              fit=meta.RiskD(data, BB.grdnum=BB.grdnum, B.sim=B.sim, cov.prob=cov.prob, print=print, studyCI=studyCI, midp=midp)
           if(type=="rate difference")
              fit=meta.RateD(data, BB.grdnum=BB.grdnum, B.sim=B.sim, cov.prob=cov.prob, print=print, studyCI=studyCI, midp=midp)
           if(type=="risk ratio")
              fit=meta.RiskR(data, BB.grdnum=BB.grdnum, B.sim=B.sim, cov.prob=cov.prob, print=print, studyCI=studyCI, ratio.upper=ratio.upper, 
midp=midp)
           if(type=="rate ratio")
              fit=meta.RateR(data, BB.grdnum=BB.grdnum, B.sim=B.sim, cov.prob=cov.prob, print=print, studyCI=studyCI, ratio.upper=ratio.upper, 
midp=midp)
            return(fit)
            }
