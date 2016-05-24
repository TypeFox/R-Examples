`vrais.LD.add` <-
function(mu,alpha.Q,s,CD,perf,DL.d)
{
 sd=sqrt(s)
 stde    = sd/CD
 dnorm.Q = dnorm(perf,mu+alpha.Q,stde)
 dnorm.q = dnorm(perf,mu-alpha.Q,stde)
 dnorm.0 = dnorm(perf,mu,stde)
 total  = dnorm.Q*DL.d[,1]+ dnorm.q*DL.d[,2]+dnorm.0*DL.d[,3]+dnorm.0*DL.d[,4]
 total[total <=0] = -750
 total[total>0]   = log(total[total>0])
 logvrais         = sum(total)
 logvrais
}

