UQnormalize <-
function( raw.count, trt, ctrl, round=TRUE, add.one=TRUE, print=FALSE){
    cnames = colnames(raw.count)
    trtcol = sapply(trt, function(ch){ which(ch==cnames) })
    ctrlcol = sapply(ctrl, function(ch){ which(ch==cnames) })
    norm.count = UQnormalize.count( raw.count[, c(trtcol,ctrlcol)], round=round, add.one=add.one, print=print)
    K = length(trt)
    norm.dat = NULL
    for (k in 1:K){
        dat.k = data.frame(sgRNA=raw.count[,1], gene=raw.count[,2], trt=norm.count[,k] , ctrl=norm.count[,(K+k)])
        norm.dat=rbind(norm.dat, dat.k)
    }
    return(norm.dat)
    }
