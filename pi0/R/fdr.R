`fdr` <-
function(p,pi0=1){
    pi0=min(max(pi0,0),1)
    G=length(p)
    G0=(G*pi0)      ## or G0=round(G*pi0) ??
    ord=order(p,decreasing=TRUE)
    qval=p[ord]*G0/G:1  ## no need to deal with ties here because cummin will take care of that later
    ans=numeric(G)
    ans[ord]=cummin(qval)
    ans
}

