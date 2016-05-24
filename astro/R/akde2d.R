akde2d = function(x, y, n = 25, lims = c(range(x), range(y)), levels = seq(0,0.9,by=0.1), ...){
    adens = kde2d(x, y, n=n, lims=lims, ...)
    asort = sort(as.numeric(adens$z))
    asum=sum(asort)
    acum = cumsum(asort)
    alevs = {}
    for(i in 1:length(levels)){
        if(levels[i]==0){
            alevs = c(alevs, 0)
        }else if(levels[i]==1){
            alevs = c(alevs, max(adens$z))
        }else{
            alevs=c(alevs,asort[max(which(acum<asum*levels[i]))])
        }
    }
    alevs = t(as.matrix(alevs))
    colnames(alevs) = levels
    out = c(adens,l=list(alevs))
    return(out)
}

