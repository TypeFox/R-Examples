"ss.fromdata.nvar" <-
function(delta,  sdhat=NULL, vhat=NULL, df=Inf, ss.ratio=1, var.ratio=1, deltaB=0, sig.level = 0.05, power = .80, 
    alternative = c("two.sided", "one.sided")) 
{

    if (is.null(vhat)) vhat<- sdhat^2
    if (is.null(sdhat)) sdhat<-sqrt(vhat)
    alternative <- match.arg(alternative)
    tside <- switch(alternative, one.sided = 1, two.sided = 2)
    if (tside==2 && deltaB !=0){stop("two-sided test should have deltaB=0") }

    calibrated.power<- 1 - find.calibrated.beta(1-power,df,sig.level/tside)

    Za<-qnorm(1-sig.level/tside)
    Zb<-qnorm(calibrated.power)

       if (delta> deltaB){ 
            if (tside==1) alternative<-paste(alternative,", delta > ",deltaB,sep="")   
             } 
        else if (delta< deltaB) {
            if (tside==1) alternative<-paste(alternative,", delta < ",deltaB,sep="")   
             }
        else{ stop(paste("delta cannot equal", deltaB)) }


    N0<- ( vhat*(ss.ratio+var.ratio)*(Za+Zb)^2 )/(ss.ratio*(delta-deltaB)^2)

    NOTE <- "n0 is number in *control* group\n\tn1 is number in *treatment* group"
    if (df==Inf) NOTE<-paste(NOTE,"\n df=infinity assumes variance is known")
    METHOD <- "Sample size for two-sample difference in normal means\n\tStandard deviation estimated from existing data"


   output<-list(n0 = ceiling(N0), n1=ceiling(ss.ratio*N0), delta = delta, sdhat=sdhat, df=df, var.ratio=var.ratio, 
        #ss.ratio=ss.ratio, 
        sig.level = sig.level, 
        power = power, 
        #calibrated.power=calibrated.power,
        alternative = alternative, note = NOTE, 
        method = METHOD)

    structure(output,class = "power.htest")
}

