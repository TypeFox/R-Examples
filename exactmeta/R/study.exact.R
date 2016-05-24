study.exact <-
function(x1, x2, n1, n2, type="risk difference", BB.grdnum=2000, cov.prob=0.95, midp=T)
          {if(type=="risk difference")
              fit=ci.RiskD(x1, x2, n1, n2, BB.grdnum=BB.grdnum, cov.prob=cov.prob, midp=midp)
           if(type=="rate difference")
              fit=ci.RateD(x1, x2, n1, n2, BB.grdnum=BB.grdnum, cov.prob=cov.prob, midp=midp)
           if(type=="risk ratio")
              fit=ci.RiskR(x1, x2, n1, n2, BB.grdnum=BB.grdnum, cov.prob=cov.prob, midp=midp)
           if(type=="rate ratio")
              fit=ci.RateR(x1, x2, n1, n2, cov.prob=cov.prob, midp=midp)
            return(fit)
            }
