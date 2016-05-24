risk.attrib.Copula<-function(tickers, data, start, end, sim.size=1000, df=10){
  
  rt=access(tickers = tickers,
            start = start,
            end = end,
            data = data)
  
  rt<-na.omit(rt)
  weight<-portfolio.optim(rt)$pw
  
  ###################
  ## VaR with t-Copula with empirical Marginal
  ###################
  p<-ncol(rt) ## no of stocks
  p.r<-p*(p-1)/2
  
  tCop <- tCopula(c(runif(p.r)), dim=p
                  , dispstr="un", df=df)
  
  ut<-pobs(rt)
  
  ut<-as.matrix(ut)
  theta.init<-rep(0,p.r)
  theta.init[p.r+1]<-10
  fit<-fitCopula(tCop, ut, method="ml"
                 , start=theta.init)
  
  param_est<-coef(fit)
  tCop_est<-tCopula(param_est[1:p.r],dim=p
                    ,dispstr="un"
                    , df=param_est[p.r+1])
  
  urt_sim<-rCopula(sim.size,tCop_est)
  
  rt_sim<-urt_sim
  
  for(i in 1:p){
    
    d<-sort(as.vector(na.omit(rt[,i])))
    k<-round(urt_sim[,i]*length(d)+1)
    rt_sim[,i]<-d[k]
    
  }
  port_rt_sim<-na.omit(rt_sim%*%weight)
  
  mctr <- mctr(tickers = tickers
               , weights = weight
               , start = start
               , end = end
               , data = data)
  
  cctr <- cctr(tickers = tickers
               , weights = weight
               , start = start
               , end = end
               , data = data)
  vol<-apply(rt,2,sd)
  
  vol_risk<-cbind(weight,vol,mctr,cctr)
  
  port.vol <- portvol(tickers = tickers
                      , weights = weight
                      , start = start
                      , end = end
                      , data = data)
  
  port.VaR<-quantile(port_rt_sim,prob=c(0.01,0.05))*100
  result<-list(vol_risk,port.vol,port.VaR)
  names(result)<-c("Volatility"
                   ,"Portfolio Volatility"
                   ,"Portfilio VaR")
  return(result)
}