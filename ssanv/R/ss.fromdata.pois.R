"ss.fromdata.pois" <-
function(xbar0,xbar1,m0,m1, ss.ratio=1, sig.level = 0.05, real.power =.80, nominal.power=NULL, 
    alternative = c("two.sided", "one.sided"),MINN0=1,MAXN0=10^5) 
{
    if (is.null(nominal.power)){
        if (real.power==.8 || real.power==.9){
            if (real.power==.8) nominal.power<-.77
            else nominal.power<-.89 
        }
        else stop("Must supply nominal.power unless real.power=.8 or .9")
    }
    else if (!is.null(real.power) && !is.null(nominal.power)){
            if  (!((real.power==.8 & nominal.power==.77) || (real.power==.9 & nominal.power==.89 ))){
                 real.power<-NULL
                 warning("Real.power ignored, only nominal power used")
            }
          }
    alternative <- match.arg(alternative)
    tside <- switch(alternative, one.sided = 1, two.sided = 2)

    root.func<-function(N0,Xbar0=xbar0,Xbar1=xbar1,M0=m0,M1=m1,R=ss.ratio,ALPHA=sig.level/tside,Power=nominal.power,epsilon=10^-6){
        Xbar0<-Xbar0+1/(2*M0)
        Xbar1<-Xbar1+1/(2*M1)
        maxt<- max(qnbinom(1-epsilon, max(1,Xbar0*M0),M0/(M0+N0)), qnbinom(1-epsilon, max(1,Xbar1*M1),M1/(M1+N0*R)) )
        maxt<-2*maxt
        pout<-0
        Qb<- qbinom(1-ALPHA,1:maxt,1/(1+R)) 
        tindex<- (1:maxt)[ Qb+1 <= (1:maxt) ] 
        dnb0<-dnbinom(0:maxt,M0*Xbar0,M0/(N0+M0))
        dnb1<-dnbinom(0:maxt,M1*Xbar1,M1/(N0*R+M1))
        for (tt in tindex){
            w<- (Qb[tt]+1):tt
	      pout<-pout+ sum( dnb0[w+1]*dnb1[tt-w+1] )
        }

        Power - pout
    }
    #if (root.func(MAXN0)>0){ uout<-MAXN0 ; warning("n0 set to MAXN0")}
    if (root.func(MINN0)<0){ uout<-MINN0; warning("n0 set to MINN0") }
    else { uout<-uniroot.integer(root.func,c(MINN0,MAXN0))$root }

    deltahat<-xbar0+1/(2*m0) - (xbar1+1/(2*m1))
    if (deltahat>0){ 
         if (tside==1) alternative<-paste(alternative,", mu0 > mu1",sep="")   
             } 
    else if (deltahat< 0) {
         if (tside==1) alternative<-paste(alternative,", mu0 < mu1",sep="")   
             }
    else{ stop("mu0hat cannot equal mu1hat") }

    NOTE <- "n0 is number in *control* group\n\tn1 is number in *treatment* group
     \n\tMethod is conservative:\n\tnominal power of .77 gives real power of .80\n\tand nominal power of .89 gives real power of .90" 
    METHOD <- "Sample size for two-sample difference in Poisson means\n\t Means estimated from existing data"

    output<-list(n0 = ceiling(uout), n1=ceiling(ss.ratio*uout), 
        xbar0=xbar0, xbar1=xbar1,
        m0=m0,m1=m1, 
        sig.level = sig.level, 
        nominal.power = nominal.power, 
        real.power=real.power,
        alternative = alternative, note = NOTE, 
        method = METHOD)

     structure(output,class = "power.htest")
}

