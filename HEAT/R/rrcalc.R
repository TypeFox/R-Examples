rrcalc <-
function (object, rrunit = 1)
{

    a <- object$best.fit
    beta_1<-a[1]
    se_1<-a[2]
    beta_2<-a[4]
    se_2<-a[5]

    fit_rr1 = exp(beta_1*rrunit)
    uci_rr1 = exp((beta_1 + 1.96 * se_1)*rrunit)
    lci_rr1 = exp((beta_1 - 1.96 * se_1)*rrunit)

    fit_rr2 = exp(beta_2*rrunit)
    uci_rr2 = exp((beta_2 + 1.96 * se_2)*rrunit)
    lci_rr2 = exp((beta_2 - 1.96 * se_2)*rrunit)

    RR <- rbind(c(fit_rr1, lci_rr1, uci_rr1),c(fit_rr2, lci_rr2, uci_rr2))
    colnames(RR) <- c("RR", "95% Lower CI.", "95% Upper CI")
    rownames(RR) <- c("<Threshold", ">=Threshold")
    return(RR)

}
