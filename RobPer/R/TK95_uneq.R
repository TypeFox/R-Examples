TK95_uneq <- function(tt, alpha=1.5) {
    if(!all(tt==sort(tt))) tt<-sort(tt)
    n<- length(tt)
    # resolution depends on minimal distance between three times
    aufl<-10*ceiling(1/ min(apply(matrix(tt[1:(2*(n%/%2))], ncol=2, byrow=TRUE),2, diff)))
    # length times ten
    N <- aufl*10
    # generate red noise
    rr_komplett<-TK95(N=N, alpha=alpha)
    # where to draw a tenth?
    stueck <- ceiling(runif(1,min=0, max=0.89)*length(rr_komplett))
    # draw a tenth
    ausgesucht <- rr_komplett[stueck:(stueck-1+aufl)]
    # distance between succesive times
    delta_gl <- (max(tt)-min(tt))/(aufl-1)
    # times to check in equally spaced sampling min(tt), min(tt)+delta,...,max(tt)-delta, max(tt)
    stellen <-(((tt-min(tt))+ delta_gl/2)%/%delta_gl)+1
    rr <- ausgesucht[stellen]
    return(rr)
}
