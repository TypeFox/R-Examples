"ss.fromdata.neff" <-
function(thetahat,m0,m1, ss.ratio=1, thetaB=0, sig.level = 0.05, real.power =.80, nominal.power=NULL, 
    alternative = c("two.sided", "one.sided"),MINN0=2,MAXN0=Inf,subdivisions=1000) 
{

    if (is.null(nominal.power)){
        if (real.power==.8 || real.power==.9){
            if (real.power==.8) nominal.power<-.76
            else nominal.power<-.88 
        }
        else stop("Must supply nominal.power unless real.power=.8 or .9")
    }
    else if (!is.null(real.power) && !is.null(nominal.power)){
            if  (!((real.power==.8 & nominal.power==.76) || (real.power==.9 & nominal.power==.88 ))){
                 real.power<-NULL
                 warning("Real.power ignored, only nominal power used")
            }
          }
    alternative <- match.arg(alternative)
    tside <- switch(alternative, one.sided = 1, two.sided = 2)
    if (tside==2 && thetaB !=0){stop("two-sided test should have thetaB=0") }

    pow<-function(n0,Theta,ThetaB,r,Alpha){
        ncp<- (Theta-ThetaB)*sqrt( (n0*r)/(r+1) )
        nu<- n0+n0*r -2
        1 - pt( qt(1-Alpha,nu), nu, ncp )
    }


    root.func<-function(N0,R=ss.ratio,Thetahat=thetahat,M0=m0,M1=m1,ALPHA=sig.level/tside,Power=nominal.power,THETAB=thetaB,J=subdivisions){
        cc<- 1/sqrt(1/M0+1/M1)
        if (N0==Inf){ powstar.value<- pt(Thetahat*cc,M0+M1-2, THETAB*cc) }
        else {
            thetai<-c(-Inf,qt((1:(J-1))/J,M0+M1-2)/cc + Thetahat,Inf)
            Tvalues<- pt(Thetahat*cc, M0+M1-2, thetai*cc )
            Pows<- pow(N0,thetai,THETAB,R,ALPHA)
            powstar.value<- sum(   (.5*Pows[-1] + .5*Pows[-(J+1)])* (Tvalues[-(J+1)] - Tvalues[-1])   )
        }
        Power - powstar.value
    }

    if (root.func(MAXN0)>0){ uout<-MAXN0 ; warning("n0 set to MAXN0")}
    else if (root.func(MINN0)<0){ uout<-MINN0; warning("n0 set to MINN0") }
    else { uout<-uniroot.integer(root.func,c(MINN0,MAXN0))$root }

    if (thetahat> thetaB){ 
         if (tside==1) alternative<-paste(alternative,", theta > ",thetaB,sep="")   
             } 
    else if (thetahat< thetaB) {
         if (tside==1) alternative<-paste(alternative,", theta < ",thetaB,sep="")   
             }
    else{ stop(paste("thetahat cannot equal", thetaB)) }

    NOTE <- "n0 is number in *control* group\n\tn1 is number in *treatment* group
        \n\tMethod is conservative:\n\tnominal power of .76 gives real power of .80\n\tand nominal power of .88 gives real power of .90" 
    METHOD <- "Sample size for two-sample difference in normal means\n\tStandardized difference in means estimated from existing data"

    output<-list(n0 = ceiling(uout), n1=ceiling(ss.ratio*uout), thetahat = thetahat,
        m0=m0,m1=m1, 
        sig.level = sig.level, 
        nominal.power = nominal.power, 
        real.power=real.power,
        alternative = alternative, note = NOTE, 
        method = METHOD)

     structure(output,class = "power.htest")
}

