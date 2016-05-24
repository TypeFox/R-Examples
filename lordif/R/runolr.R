runolr <-
function(rv,ev,gr,wt) {
    ng<-length(table(gr))
    nobs<-length(rv)
    md1<-lrm(rv ~ ev, weights=wt)
    md2<-lrm(rv ~ ev + gr, weights=wt)
    md3<-lrm(rv ~ ev*gr, weights=wt)
    beta12<-round(abs((md2$coefficients[["ev"]]-md1$coefficients[["ev"]])/md1$coefficients[["ev"]]),4)
    deviance0<-md1$deviance[1]
    deviance1<-md1$deviance[2]
    deviance2<-md2$deviance[2]
    deviance3<-md3$deviance[2]
    pseudo1.CoxSnell<-1-exp(diff(md1$deviance)/nobs)
    pseudo2.CoxSnell<-1-exp(diff(md2$deviance)/nobs)
    pseudo3.CoxSnell<-1-exp(diff(md3$deviance)/nobs)
    pseudo1.Nagelkerke<-pseudo1.CoxSnell/(1-exp(-deviance0/nobs))
    pseudo2.Nagelkerke<-pseudo2.CoxSnell/(1-exp(-deviance0/nobs))
    pseudo3.Nagelkerke<-pseudo3.CoxSnell/(1-exp(-deviance0/nobs))
    pseudo1.McFadden<-1-deviance1/deviance0
    pseudo2.McFadden<-1-deviance2/deviance0
    pseudo3.McFadden<-1-deviance3/deviance0
    df12<-ng-1; df13<-2*(ng-1); df23<-ng-1;
    chi12<-round(1-pchisq(deviance1-deviance2,df12),4)
    chi13<-round(1-pchisq(deviance1-deviance3,df13),4)
    chi23<-round(1-pchisq(deviance2-deviance3,df23),4)
    pseudo12.CoxSnell<-round(pseudo2.CoxSnell-pseudo1.CoxSnell,4)
    pseudo13.CoxSnell<-round(pseudo3.CoxSnell-pseudo1.CoxSnell,4)
    pseudo23.CoxSnell<-round(pseudo3.CoxSnell-pseudo2.CoxSnell,4)
    pseudo12.Nagelkerke<-round(pseudo2.Nagelkerke-pseudo1.Nagelkerke,4)
    pseudo13.Nagelkerke<-round(pseudo3.Nagelkerke-pseudo1.Nagelkerke,4)
    pseudo23.Nagelkerke<-round(pseudo3.Nagelkerke-pseudo2.Nagelkerke,4)
    pseudo12.McFadden<-round(pseudo2.McFadden-pseudo1.McFadden,4)
    pseudo13.McFadden<-round(pseudo3.McFadden-pseudo1.McFadden,4)
    pseudo23.McFadden<-round(pseudo3.McFadden-pseudo2.McFadden,4)
    return(list(chi12=chi12,chi13=chi13,chi23=chi23,beta12=beta12,
                pseudo1.CoxSnell=pseudo1.CoxSnell,pseudo2.CoxSnell=pseudo2.CoxSnell,pseudo3.CoxSnell=pseudo3.CoxSnell,
                pseudo1.Nagelkerke=pseudo1.Nagelkerke,pseudo2.Nagelkerke=pseudo2.Nagelkerke,pseudo3.Nagelkerke=pseudo3.Nagelkerke,
                pseudo1.McFadden=pseudo1.McFadden,pseudo2.McFadden=pseudo2.McFadden,pseudo3.McFadden=pseudo3.McFadden,
                pseudo12.CoxSnell=pseudo12.CoxSnell,pseudo13.CoxSnell=pseudo13.CoxSnell,pseudo23.CoxSnell=pseudo23.CoxSnell,
                pseudo12.Nagelkerke=pseudo12.Nagelkerke,pseudo13.Nagelkerke=pseudo13.Nagelkerke,pseudo23.Nagelkerke=pseudo23.Nagelkerke,
                pseudo12.McFadden=pseudo12.McFadden,pseudo13.McFadden=pseudo13.McFadden,pseudo23.McFadden=pseudo23.McFadden,
                df12=df12,df13=df13,df23=df23))
  }
