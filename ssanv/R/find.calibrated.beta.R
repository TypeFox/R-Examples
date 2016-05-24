"find.calibrated.beta" <-
function(beta,df,alpha=.05){
    if (df==Inf){ calibrated.beta<-beta }
    else{
    root.func<-function(ROOT.BETA,NOMINAL.BETA=beta,ROOT.DF=df,ROOT.ALPHA=alpha){
        Power.given.x<-function(x,POWER.DF=ROOT.DF,POWER.ALPHA=ROOT.ALPHA,POWER.BETA=ROOT.BETA){
          pnorm(   - qnorm( 1- POWER.ALPHA) + (qnorm(1-POWER.ALPHA) + qnorm(1-POWER.BETA) )*( sqrt(x)/sqrt(POWER.DF) ) ) * dchisq(x,POWER.DF)
        }
        REAL.POWER<-integrate(Power.given.x,0,ROOT.DF)$value + integrate(Power.given.x,ROOT.DF,Inf)$value 
        (1-NOMINAL.BETA) - REAL.POWER
    }
    calibrated.beta<-uniroot(root.func,c(10^-10,.5))$root
    }
    calibrated.beta
}

