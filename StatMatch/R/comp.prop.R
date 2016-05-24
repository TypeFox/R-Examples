comp.prop <- function(p1, p2, n1, n2=NULL, ref=FALSE){
    if(!ref && is.null(n2)) stop("If p2 is not the reference distribution (ref=FALSE) please \n provide the n2 argument (sample size)")
    if(sum(p1)>1) p1 <- prop.table(p1)
    if(sum(p2)>1) p2 <- prop.table(p2)

    tvd <- 0.5*sum(abs(p1-p2))
    ov <- 1 - tvd
    bhatt <- sum(sqrt(p1*p2))
    hell <-  sqrt(1-bhatt)
    dd <- c("tvd"=tvd,"overlap"=ov, "Bhatt"=bhatt, "Hell"=hell)
#    
    if(ref){
        J <- length(p2)
        not0 <- p2>0
        pp1 <- p1[not0]
        pp2 <- p2[not0]
        J <- length(pp2)
        chi.p <- n1*sum((pp1-pp2)^2/pp2)
        pe <- p2
    }
    else{
      w1 <- n1/(n1+n2)
      pe <- p1*w1 + p2*(1 - w1)
      not0 <- pe>0
      ppe <- pe[not0]
      pp1 <- p1[not0]
      pp2 <- p2[not0]
      J <- length(ppe)
      chi.1 <- n1*sum( (pp1-ppe)^2/ppe)
      chi.2 <- n2*sum( (pp2-ppe)^2/ppe)
      chi.p <- chi.1 + chi.2  # Pearson's Chi-square
    }
    chis <- c(Pearson=chi.p, df=J-1, q0.05=qchisq(0.05, df=J-1, lower.tail=FALSE), 
              delta.h0=chi.p/qchisq(0.05, df=J-1, lower.tail=FALSE))
    list(meas=dd, chi.sq=chis, p.exp=pe)
}
