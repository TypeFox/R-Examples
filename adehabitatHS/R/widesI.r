"widesI" <- function(u, a, avknown=TRUE, alpha=0.05)
{
    ## Verifications
    if (length(u)!=length(a))
        stop("available and used vector should have the same length")
    if (is.null(names(u)))
        names(u)<-paste("Habitat", 1:length(u), sep="")


    ## Bases
    sorties<-list()
    ui<-u
    ai<-a

    oi<-ui/sum(ui) # proportion used
    pi<-ai/sum(ai) # proportion available
    I<-length(u)   # number of habitats
    wi<-oi/pi      # selection ratios
    bonferroni<-alpha/(I*(I-1)/2) # level of Bonferroni confidence intervals


    ## Preparation of output
    sorties$used.prop<-oi
    sorties$se.used<-sqrt((oi*(1-oi)/sum(ui)))
    sorties$avail.prop<-pi
    sorties$se.avail<-sqrt(pi*(1-pi)/sum(a))
    sorties$wi<-wi
    if (avknown) {
      sorties$se.wi<-sqrt(oi*(1-oi)/(sum(u)*(pi^2)))
    } else {
      sorties$se.wi<-wi*sqrt(1/ui-1/sum(ui)+1/ai-1/sum(ai))
    }

    ## Chi-square of the selection ratios = 1
    testwi<-((sorties$wi-1)/sorties$se.wi)^2
    sorties$chisquwi<-data.frame(testwi=testwi,
                                 p=1-pchisq(testwi, 1))

    ## Standardised selection ratios
    sorties$Bi<-sorties$wi/sum(sorties$wi)

    if (avknown) {

        ## Classical Chi-square
        sorties$Khi2P<-c(tmp<-sum(((ui-sum(u)*pi)^2)/(sum(u)*pi)),
                         length(u)-1,
                         1-pchisq(tmp, I-1))

        ## Chi-square based on likelihood
        tmp<-u*log(u/(sum(u)*pi))
        tmp[is.na(tmp)]<-0
        sorties$Khi2L<-c(tmp<-2*sum(tmp),
                         length(u)-1,
                         1-pchisq(tmp, I-1))
    } else {

        Eui<-(ai+ui)*sum(ui)/(sum(u)+sum(ai))
        Eai<-(ai+ui)*sum(ai)/(sum(u)+sum(ai))

        ## Classical Chi-square
        sorties$Khi2P<-c(tmp<-sum(((ui-Eui)^2)/Eui + ((ai-Eai)^2)/Eai),
                         length(u)-1,
                         1-pchisq(tmp, I-1))

        ## Chi-square based on likelihood
        sorties$Khi2L<-c(tmp<-2*sum(u*log(u/Eui)+ai*log(ai/Eai)),
                         length(u)-1,
                         1-pchisq(tmp, I-1))
    }

    ## output
    names(sorties$Khi2P)<-c("Khi2P", "df", "pvalue")
    names(sorties$Khi2L)<-c("Khi2L", "df", "pvalue")

    ## Bases for the computation of the
    ## ranking matrix of the selection ratios
    diffwi<-matrix(0, nrow=I, ncol=I)
    vardif<-matrix(0, nrow=I, ncol=I)
    sig<-matrix(0, nrow=I, ncol=I)
    ICdiffupper<-matrix(0, nrow=I, ncol=I)
    ICdifflower<-matrix(0, nrow=I, ncol=I)
    sig<-matrix(0, nrow=I, ncol=I)

    ## filling the ranking matrix
    for (i in 1:I) {
        for (j in 1:I) {
            if (i!=j) {

                ## variance of the difference of selection ratios
                vardif[i,j]<-ifelse(avknown,
                                    (oi[i]*(1-oi[i])/(sum(ui)*(pi[i]^2))+
                                     oi[j]*(1-oi[j])/(sum(ui)*(pi[j]^2))-
                                     2*oi[i]*oi[j]/(sum(ui)*pi[i]*pi[j])),
                                    (wi[i]/pi[i]+wi[j]/
                                     pi[j]-((wi[i]-wi[j])^2))/sum(ui)+
                                    ((wi[i]^2)/pi[i]+(wi[j]^2)/
                                     pi[j]-((wi[i]-wi[j])^2))/sum(ai))

                ## difference of selection ratios
                diffwi[i,j]<-(wi[i]-wi[j])

                ## Coonfidence intervals
                ICdiffupper[i,j]<-round(diffwi[i,j]+
                                        sqrt(vardif[i,j])*
                                        qnorm(1-bonferroni/2),4)
                ICdifflower[i,j]<-round(diffwi[i,j]-
                                        sqrt(vardif[i,j])*
                                        qnorm(1-bonferroni/2),4)

                ## The ranking matrix
                if (diffwi[i,j]<0) sig[i,j] <- "-"
                if (diffwi[i,j]>0) sig[i,j] <- "+"
                if (diffwi[i,j]==0) sig[i,j] <- "0"
                if (ICdiffupper[i,j]<0)
                    sig[i,j]<-"---"
                if (ICdifflower[i,j]>0)
                    sig[i,j]<-"+++"
            } else {
                sig[i,j]<-"0"
            }
        }
    }

    ## row and column names
    rownames(diffwi)<-names(u)
    colnames(diffwi)<-names(u)
    rownames(ICdiffupper)<-names(u)
    colnames(ICdiffupper)<-names(u)
    rownames(ICdifflower)<-names(u)
    colnames(ICdifflower)<-names(u)
    rownames(sig)<-names(u)
    colnames(sig)<-names(u)

    ## Output
    sorties$avknown<-avknown
    sorties$comparisons$diffwi<-diffwi
    sorties$comparisons$ICdiffupper<-ICdiffupper
    sorties$comparisons$ICdifflower<-ICdifflower
    sorties$comparisons$signif<-sig
    ## the profile
    sorties$profile<-profilehab(sig, wi)
    sorties$alpha<-alpha
    class(sorties)<-c("wiI", "wi")
    return(sorties)
}

