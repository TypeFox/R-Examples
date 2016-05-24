hdr.den <- function(x=NULL, prob=c(50,95,99), den=NULL, h=hdrbw(BoxCox(x,lambda),mean(prob)), lambda=1, xlab=NULL, ylab="Density",...)
{
    if(missing(den))
        den <- tdensity(x,bw=h,lambda=lambda)
    else if(missing(x))
     x <- sample(den$x, 500, replace=TRUE, prob=den$y)
    hdr <- hdr(x,prob,den,h)
    maxden <- max(den$y)
    plot(den,type="l",xlab=xlab,ylab=ylab,...)
    cols <- rep(c(0, 4, 2, 3), 3)
    nregions <- nrow(hdr$hdr)
    for(i in 1:nregions)
    {
        lines(range(den$x),rep(hdr$falpha[i],2),col=5)
        for(j in 1:length(hdr$hdr[i,]))
            lines(rep(hdr$hdr[i,j],2),c((0.01+(i-1)*0.02)*maxden,hdr$falpha[i]),col=6)
    }
    for(i in 1:nrow(hdr$hdr))
        add.hdr(hdr$hdr[i,], (0.01+(i-1)*0.02)*maxden, 0.02*maxden, col=cols[i+1], horiz=TRUE)
    return(hdr)
}
