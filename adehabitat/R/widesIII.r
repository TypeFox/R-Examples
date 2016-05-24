"widesIII" <- function(u, a, avknown = TRUE, alpha = 0.05)
{
    ## Verifications
    u<-as.matrix(u)
    a<-as.matrix(a)
    if (nrow(u) != nrow(a))
      stop("available and used matrix should have the same number of animals")
    if (ncol(u) != ncol(a))
      stop("available and used matrix should have the same number of habitats")

    ## Bases
    sorties<-list()
    pij<-as.matrix(a)

    ## Computation of the availability if not given in percentage
    aip<-apply(a,2,sum)
    apj<-apply(a,1,sum)

    ## Computation of the use if not given in percentage
    uij<-as.matrix(u)
    if (is.null(colnames(u)))
        colnames(uij) <- paste("Habitat", 1:ncol(u), sep = "")
    if (is.null(colnames(a)))
        colnames(pij) <- paste("Habitat", 1:ncol(a), sep = "")

    ## The two matrices
    pij<-as.matrix(a/apj)
    uij<-as.matrix(u)
    I<-ncol(uij)
    J<-nrow(uij)

    ## Calcul de l'IC de Bonferroni
    bonferroni <- alpha/(I * (I - 1)/2)
    upj<-apply(uij,1,sum)
    uip<-apply(uij,2,sum)
    wij<-uij/(upj*pij)
    wi<-uip/apply(pij*upj,2,sum)

    ## Output
    sorties$used.prop <- t(t(uij)/uip)
    sorties$avail.prop <- pij
    sorties$wij<-wij
    sorties$wi<-wi

    ## Computation of the Khi2
    Khi2Lj<-matrix(0, nrow=J, ncol=3)
    colnames(Khi2Lj)<-c("Khi2Lj", "df", "pvalue")
    for (j in 1:J) {
        euij<-uij[j,]*log(uij[j,]/(upj[j]*pij[j,]))
        ddl<-length(euij[!is.na(euij)])-1
        euij<-euij[!is.na(euij)]
        Khi2Lj[j,1]<-sum(euij)
        Khi2Lj[j,2]<-ddl
        Khi2Lj[j,3]<-1 - pchisq(Khi2Lj[j,1], ddl)
    }
    rownames(Khi2Lj)<-rownames(u)

    ## Output
    sorties$Khi2Lj<-Khi2Lj
    Khi2L<-apply(Khi2Lj,2,sum)
    Khi2L[3]<-1 - pchisq(Khi2L[1],Khi2L[2])
    names(Khi2L)<-c("Khi2L", "df", "pvalue")
    sorties$Khi2L<-Khi2L

    ## Variance of the selection ratios
    vwi<-rep(0,I)
    for (i in 1:I) {
        yj<-uij[,i]
        xj<-pij[,i]*upj
        vwi[i]<-(sum((yj-wi[i]*xj)**2)/(J-1))*(1/(J*(mean(xj)**2)))
    }
    sewi<-sqrt(vwi)

    ## Output
    sorties$se.wi<-sewi
    sorties$ICwiupper<-round(wi+sewi*qnorm(1 - alpha/(2*I)), 4)
    sorties$ICwilower<-round(wi-sewi*qnorm(1 - alpha/(2*I)), 4)

    ## Matrix of the differences of the selection ratios, the ranking matrix
    diffwi<-outer(wi,wi,"-")
    vardif<-matrix(0, I, I)
    ICdiffupper<-matrix(0, I, I)
    ICdifflower<-matrix(0, I, I)
    sig<-matrix("0", I, I)

    for (i in 1:I) {
        for (j in 1:I) {
            if (avknown) {
                ## availability known
                spi<-sum(pij[,i]*upj)
                spj<-sum(pij[,j]*upj)

                ## matrix of the variances
                vardif[i,j]<-sum(((uij[,i]-wi[i]*upj)/spi +
                                  (uij[,j]-wi[j]*upj)/spj )**2 )*(J/(J-1))
            } else {
                ## unknown availability
                dftmp<-data.frame(y1=uij[,i],y2=uij[,j],
                                  x1=pij[,i]*upj, x2=pij[,j]*upj)
                vc<-var(dftmp)
                y1<-uip[i]
                y2<-uip[j]
                x1<-sum(pij[,i]*upj)
                x2<-sum(pij[,j]*upj)
                vardif[i,j]<-(1/(y1**2))*vc["x1","x1"] +
                    ((x1**2)/(y1**4))*vc["y1","y1"] +
                        (1/(y2**2))*vc["x2","x2"] +
                            ((x2**2)/(y2**4))*vc["y2","y2"] -
                                2*(x1/(y1**3))*vc["x1","y1"] -
                                    2*(1/(y1*y2))*vc["x1","x2"] +
                                        2*(x2/(y1*(y2**2)))*vc["x1","y2"] +
                                            2*(x1/(y2*(y1**2)))*vc["y1","x2"] -
                                                2*((x1*x2)/
                                                   ((y1**2)*(y2**2))) *
                                                       vc["y1","y2"] -
                                                           2*(x2/(y2**3))*
                                                               vc["x2","y2"]
            }
            vardif[row(vardif)==col(vardif)]<-0

            ## Confidence intervals on the differences of selection ratios
            ICdiffupper[i, j] <- round(diffwi[i, j] +
                                       sqrt(vardif[i,j]) *
                                       qnorm(1 - bonferroni/2), 4)
            ICdifflower[i, j] <- round(diffwi[i, j] -
                                       sqrt(vardif[i,j]) *
                                       qnorm(1 - bonferroni/2), 4)

            ## Ranking matrix
            sig[i, j] <- ifelse(diffwi[i, j] < 0, "-", "+")
            if (ICdiffupper[i, j] < 0)
                sig[i, j] <- "---"
            if (ICdifflower[i, j] > 0)
                sig[i, j] <- "+++"
        }
    }

    ## Row and column names
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
    class(sorties) <- c("wiIII", "wi")
    return(sorties)
  }

