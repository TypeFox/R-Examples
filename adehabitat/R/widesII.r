"widesII" <- function(u, a, avknown=TRUE, alpha=0.05)
{
    ## Verifications
    u<-as.matrix(u)
    if (ncol(u)!=length(a))
        stop("used and available matrices should have the same number of habitats")

    ## Bases
    uij<-as.matrix(u)
    ai<-a
    pi<-ai/sum(ai)
    n<-nrow(uij)
    I<-ncol(uij)
    uip<-apply(uij,2,sum)
    upj<-apply(uij,1,sum)
    upp<-sum(as.vector(uij))
    oi<-uip/upp
    sorties<-list()
    sorties$used.prop <- oi
    sorties$avail.prop <- pi

    ## Computation of Khi2L1
    Euij<-outer(upj,uip)/upp
    tmp<-log(uij/Euij)
    tmp[abs(tmp)==Inf]<-0
    tmp[is.na(tmp)]<-0
    Khi2L1<-c(tmp<-2*sum(as.vector(uij*tmp)), df<-(I-1)*(n-1),
               1 - pchisq(tmp, df))
    names(Khi2L1)<-c("Khi2L1", "df", "pvalue")
    sorties$Khi2L1<-Khi2L1

    ## Computation of Khi2L2
    if (avknown) {
        Euij<-outer(upj,pi)
        tmp<-log(uij/Euij)
        tmp[abs(tmp)==Inf]<-0
        Khi2L2<-c(tmp<-2*sum(as.vector(uij*tmp)), df<-(I-1)*n,
                  1 - pchisq(tmp, df))
        names(Khi2L2)<-c("Khi2L2", "df", "pvalue")
        sorties$Khi2L2<-Khi2L2


        ## The difference between the two Khi squares
        Khi2L2MinusL1<-c(tmp<-Khi2L2[1]-Khi2L1[1], I-1,
                         1 - pchisq(tmp, I-1))
        names(Khi2L2MinusL1)<-c("Khi2L2MinusL1", "df", "pvalue")
        sorties$Khi2L2MinusL1<-Khi2L2MinusL1

    } else {

        ## Khi2L2 when availability is estimated
        uija<-rbind(uij,ai)
        upja<-apply(uija, 1, sum)
        uipa<-apply(uija, 2, sum)
        Euija<-outer(upja,uipa)/sum(upja)
        tmp<-log(uija/Euija)
        tmp[abs(tmp)==Inf]<-0
        Khi2L2<-c(tmp<-2*sum(as.vector(uija*tmp)), df<-(I-1)*n,
                  1 - pchisq(tmp, df))
        names(Khi2L2)<-c("TrickyKhi2", "df", "pvalue")
        sorties$Khi2L2<-Khi2L2

        ## Computation of the difference between the two Khi square
        Khi2L2MinusL1<-c(tmp<-Khi2L2[1]-Khi2L1[1], I-1,
                         1 - pchisq(tmp, I-1))
        names(Khi2L2MinusL1)<-c("Khi2L2MinusL1", "df", "pvalue")
        sorties$Khi2L2MinusL1<-Khi2L2MinusL1
    }

    ## Matrix of selection ratios
    wij<-t(t(uij/upj)/pi)
    wi<-(uip/upp)/pi
    sorties$wij<-wij
    sorties$wi<-wi

    ## Variance of selection ratios
    if (avknown) {
        varwi<-apply((((t(t(uij)/pi) - outer(upj,wi) )^2)/(n-1)),
                     2, sum)*(n/(upp^2))
        sewi<-sqrt(varwi)
    }  else {
        Vi<-uip/upp
        varVi<-Vi
        for (i in 1:length(Vi)) {
            varVi[i]<-(sum((u[,i]-Vi[i]*upj)^2 )/(n-1))/(n*(mean(upj)^2))
        }
        varpi<-pi*(1-pi)/sum(ai)
        sewi<-sqrt(((Vi/pi)^2)*(varVi/(Vi^2)+varpi/(pi^2)))
    }

    ## Output
    sorties$se.wi<-sewi
    sorties$ICwiupper<-round(wi+sewi*qnorm(1 - alpha/(2*I)), 4)
    sorties$ICwilower<-round(wi-sewi*qnorm(1 - alpha/(2*I)), 4)


    ## Matrix of the standard errors of
    ## the differences of the selection ratios
    diffwi<-outer(wi,wi,"-")
    sediffwi<-diffwi

    if (avknown) {
        for (i in 1:I) {
            for (j in 1:I) {
                tmp<-uij[,i]/pi[i] - uij[,j]/pi[j] - wi[i]*upj + wi[j]*upj
                sediffwi[i,j]<-sqrt(((n/(n-1))/(upp^2))*sum(tmp^2))
            }
        }
    } else {
        for (i in 1:I) {
            for (j in 1:I) {
                tmp<-(sum((uij[,i]/pi[i]-uij[,j]/
                           pi[j]-diffwi[i,j]*upj)^2)/(n-1))*(n/(upp^2))
                tmp<-tmp+((wi[i]^2)/pi[i]+(wi[j]^2)/
                          pi[j]-(diffwi[i,j]^2))/sum(ai)
                sediffwi[i,j]<-sqrt(tmp)
            }
        }
    }

    ## Confidence intervals on these differences
    bonferroni <- alpha/(I * (I - 1)/2)
    ICdiffupper<-round(diffwi+sediffwi*qnorm(1 - bonferroni/2), 4)
    ICdifflower<-round(diffwi-sediffwi*qnorm(1 - bonferroni/2), 4)

    ## The ranking matrix
    sig<-diffwi
    for (i in 1:I) {
        for (j in 1:I) {
            if (i!=j) {
                sig[i, j] <- ifelse(diffwi[i, j] < 0, "-", "+")
                if (ICdiffupper[i, j] < 0)
                    sig[i, j] <- "---"
                if (ICdifflower[i, j] > 0)
                    sig[i, j] <- "+++"
            }
            else {
                sig[i,j]<-"0"
            }
        }
    }

    ## Rownames and colnames
    rownames(diffwi) <- colnames(u)
    colnames(diffwi) <- colnames(u)
    rownames(ICdiffupper) <- colnames(u)
    colnames(ICdiffupper) <- colnames(u)
    rownames(ICdifflower) <- colnames(u)
    colnames(ICdifflower) <- colnames(u)
    rownames(sig) <- colnames(u)
    colnames(sig) <- colnames(u)

    ## Output
    sorties$avknown <- avknown
    sorties$comparisons$diffwi <- diffwi
    sorties$comparisons$ICdiffupper <- ICdiffupper
    sorties$comparisons$ICdifflower <- ICdifflower
    sorties$comparisons$signif <- sig
    sorties$profile <- profilehab(sig, wi)
    sorties$alpha <- alpha
    class(sorties) <- c("wiII", "wi")
    return(sorties)
}

