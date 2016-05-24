`renyicomp` <-
function(x,y,factor,sites=Inf,scales=c(0,0.25,0.5,1,2,4,8,Inf),permutations=100,plotit=T,...) {
    groups <- table(y[,factor])
    if (sites == Inf) {sites <- min(groups)}
    m <- length(groups)
    n <- max(groups)
    s <- length(scales)
    levels <- names(groups)
    result <- array(NA,dim=c(m,s,6))
    dimnames(result) <- list(level=levels,scale=scales,c("mean","stdev","min","max","Qnt 0.025","Qnt 0.975"))
    names(dimnames(result)) <- c(factor,"scale","")
    for (i in 1:m) {
       if (groups[i] == sites) {result[i,,1] <- as.matrix(renyiresult(x,y,factor,levels[i],scales=scales))}
       if (groups[i] > sites) {result[i,,] <- renyiaccumresult(x,y,factor,levels[i],scales=scales,permutations=permutations)[sites,,]}
    }
    if (plotit==T) {renyiplot(result[,,1],...)}
    return(result)
}

