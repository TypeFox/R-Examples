Qt <- function(
   m,
   n,
   sigma2)
{

#NUMBER OF STUDIES
    t <- length(m)
#GLOBAL MEAN
    M <- sum(m*n)/sum(n)
    m.star <- m-M
    n.bar <- mean(n)

    Qd <- sum(m.star^2/sigma2)
    Qe <- sum(m.star^2*n/(n.bar*sigma2))

    bar.Qd <- 1/mean(sigma2)*sum(m.star^2)
    bar.Qe <- 1/(mean(sigma2)*n.bar)*sum(n*m.star^2)

    Md <- sum(m/sigma2)/sum(1/sigma2)
    Me <- sum(m*n/(n.bar*sigma2))/sum(n/(n.bar*sigma2))

    tilde.Qd <- sum((m-Md)^2/sigma2)
    tilde.Qe <- sum(n*(m-Me)^2/(n.bar*sigma2))

return(
       list(
           t=t,
           Qd=Qd,
           Qe=Qe,
           bar.Qd=bar.Qd,
           bar.Qe=bar.Qe,
           tilde.Qd=tilde.Qd,
           tilde.Qe=tilde.Qe
           )
       )
}
