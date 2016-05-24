beta.indices <- function(
             ipd.f,
             ad.f)
{

    ipd.terms <- terms(ipd.f)
    ad.terms <- terms(ad.f)

    ipd.labels = attr(ipd.terms,"term.labels")
    ad.labels = attr(ad.terms,"term.labels")

    if(attr(ad.terms,"intercept")) ad.labels = c("intercept",ad.labels)

    p.ipd <- length(ipd.labels) 
    p.ad <- length(ad.labels)

    ipd.index <- 1:p.ipd
    ad.index <- (p.ipd+1):(p.ipd+p.ad)

#SHARED PARAMETERS

    ipd.position <- which.formula(ipd.f,ad.f)
    ad.position <- which.formula(ad.f,ipd.f)+attr(ad.terms,"intercept")

    ad.index[ad.position] <- ipd.position
    ad.index[-ad.position] <- (p.ipd+1):(p.ipd+length(ad.index[-ad.position]))

    shared.index = ad.index<=max(ipd.index)

    labels = beta.construct(ipd.labels,ad.labels,shared.index)

return(
       list(
        ipd = ipd.index,
        ad = ad.index,
        shared.index = shared.index,
        labels = labels)
      )
}

beta.construct <- function(
               beta.ipd,
               beta.ad,
               shared.index)
{
  c(beta.ipd,beta.ad[!shared.index])
}






