power2x2 <-
function(p0,p1,n0,n1=NULL,sig.level=0.05,alternative=c("two.sided","one.sided"),paired=FALSE,strict=FALSE,tsmethod=NULL,nullOddsRatio=1,errbound=10^-6,approx=FALSE){
    prob.reject<-0

    if (p0>1 | p1>1 | p0<0 | p1<0) stop("p0 and p1 should be between 0 and 1")
    if (is.null(n1)) n1<-n0
    if (paired & any(n0 != n1)) stop("n0 must equal n1 when paired=TRUE")

    if (paired) warning("Assumes independent binomials for power calculations, may not be a good assumption for paired=TRUE (McNemar's) tests")
    # use same argument names as in exact2x2, but do not want double naming when calling
    PAIRED<-paired
    TSMETHOD<-tsmethod
    eps<-errbound
    alt<-match.arg(alternative)
    ## SIG.LEVEL and ALT are values that will be used in calculation
    ## not output
    SIG.LEVEL<-sig.level
    ALT<-alt
    if (!strict & ALT=="two.sided"){
        # for two.sided (strict==FALSE) you want to calculate power on one-sided p-values at alpha/2
        # so you do not count rejections in the wrong direction
        SIG.LEVEL<- SIG.LEVEL/2
        ALT<-"one.sided"
    }
    if (ALT=="one.sided"){
        if (p1<p0){ ALT<-"less"
        } else { ALT<-"greater" } 
    }
    if (approx){
        if (paired==TRUE) stop("no approximate method written for McNemar's test yet")
        Crit<- qnorm(1-SIG.LEVEL)
        delta<-p0-p1
        ## see e.g., Chow, Shao and Wang, 2008,p.89. SS Calcs in CLin Research, 2nd edition
        #V0<-p0*(1-p0)/n0
        #V1<-p1*(1-p1)/n1
        #prob.reject<- pnorm( abs(delta)/sqrt(V0+V1) - Crit )


        ## see Fleiss, 1981, p. 44

        r<- n1/n0
        m<-n0
        pbar<- (p0 + r*p1)/(r+1)
        se<- sqrt( pbar*(1-pbar)*(r+1)/(m*r) )
        z<- (abs(delta) - (1/(2*m))*((r+1)/r)  )/se
        prob.reject<- pnorm( z - Crit )

    } else {
        ilow<-max(0,qbinom(eps/2,n0,p0)-1)
        ihigh<-min(n0,qbinom(1-eps/2,n0,p0)+1)
        jlow<-max(0,qbinom(eps/2,n1,p1)-1)
        jhigh<-min(n1,qbinom(1-eps/2,n1,p1)+1)
        for (i in ilow:ihigh){
            for (j in jlow:jhigh){
                # it would be slightly faster not put data in matrix form, but to keep from making programming errors,
                # just call exact2x2
                x<-matrix(c(n0-i,i,n1-j,j),2,2)
                pval<-exact2x2(x,or=nullOddsRatio, alternative=ALT, tsmethod=TSMETHOD, conf.int=FALSE, paired=PAIRED, plot=FALSE)$p.value
                if (pval<=SIG.LEVEL){ 
                    prob.reject<-prob.reject+ dbinom(i,n0,p0)*dbinom(j,n1,p1) 
                }
            }
        }
    }
    ## create pretty output using power.htest class
    if (is.null(TSMETHOD)) TSMETHOD<-"NULL (see exact2x2 help)"
    if (paired==FALSE & strict==FALSE){
        METHOD<-"Power for Fisher's Exact Test"
    } else if (paired==TRUE & strict==FALSE){
        METHOD<-"Power for McNemar's Exact Test"
    } else if (paired==FALSE & strict){
        METHOD<-paste("Power for exact conditional test including power to reject in the wrong direction, tsmethod=",TSMETHOD)
    } else if (paired==TRUE & strict){
        METHOD<-"Power for McNemar's Exact Test, including power to reject in the wrong direction"
    }

    if (approx) METHOD<- paste("Approximate",METHOD)

    ## NOTE: output original sig.level NOT SIG.LEVEL
    ##    e.g., two.sided sig.level=0.05 should be calculated (when strict=FALSE)
    ##    at level SIG.LEVEL=0.025 
    output<-list(power=prob.reject, n0 =n0, n1=n1, 
        p0=p0,p1=p1, sig.level = sig.level, 
        alternative =alt, nullOddsRatio=nullOddsRatio, note =paste("errbound=",errbound), 
        method = METHOD)

    structure(output,class = "power.htest")
}

#pow<-power2x2(.3,.8,12)
#str(pow)
