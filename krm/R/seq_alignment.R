alignment2count=function (alignment, level=20, weight=rep(1,nrow(alignment))){
    #t(sapply(1:ncol(alignment), function (i) table(factor(alignment[,i], levels=1:level)) ))
    n=nrow(alignment)
    T=ncol(alignment)
    out=matrix(0,T,level)
    for (t in 1:T)
        for (i in 1:n) {
            aa=alignment[i,t]
            if (aa<=level) out[t,aa]=out[t,aa] + weight[i]
        }
    out
}

alignment2trancount=function (alignment, weight=rep(1,nrow(alignment))){
    n=nrow(alignment)
    counts=matrix(0, ncol(alignment), 4)    
    for (i in 1:n) {
        x=alignment[i,,drop=F]
        if(x[1]!=21) counts[1,1]=counts[1,1]+weight[i]
        else counts[1,2]=counts[1,2]+weight[i]
        if (ncol(alignment)>1) {
            for (t in 2:ncol(alignment)) {
                if (x[t-1]!=21) 
                    if(x[t]!=21) counts[t,1]=counts[t,1]+weight[i]
                    else counts[t,2]=counts[t,2]+weight[i]
                else 
                    if(x[t]!=21) counts[t,3]=counts[t,3]+weight[i]
                    else counts[t,4]=counts[t,4]+weight[i]
            }
        }
    }
    counts
}

# remove gap from a seq
removeGap=function (seq) {
    tmp=aa2arabic(seq)
    tmp1=tmp[tmp!=21]
    concatList(aaList[tmp1])
}

calcPairwiseIdentity=function (alignment, dissimilarity, removeGap) {
    n=nrow(alignment)
    if (n>2) {
        out=matrix(0, nrow=n, ncol=n)
        for (i in 1:n) 
            for (j in 1:n)
                out[i,j]=calcPairwiseIdentity(alignment[c(i,j),], dissimilarity, removeGap)
        out 
    } else {
        if (removeGap) {
            tmpAlign=alignment[,!(alignment[1,]==21 & alignment[2,]==21)] # remove positions that are both gaps
        } else {
            tmpAlign=alignment
        }
        tmp=mean(tmpAlign[1,]==tmpAlign[2,])
        tmp=(tmp*100)
        if (dissimilarity) 100-tmp
        else tmp
    }            
}
