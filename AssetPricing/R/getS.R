getS <- function(dS,xc,t,jays) {
    if(is.list(dS)) {
        S <- lapply(dS[jays],function(f,x,t){f(x,t)},x=xc,t=t)
    } else {
        S <- lapply(jays,function(j,dS,x,t){dS(x,t,j)},dS=dS,x=xc,t=t)
    }
    matrix(unlist(S),nrow=length(xc),ncol=length(jays))
}
