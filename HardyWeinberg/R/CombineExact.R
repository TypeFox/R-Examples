CombineExact <-function (Xmat, alternative="less")
{
    m <- nrow(Xmat)
    pvecless <- HWExactMat(Xmat,pvaluetype="midp",alternative="less",verbose=FALSE)$pvalvec
#    pvecgrea <- HWExactMat(Xmat,pvaluetype="midp",alternative="greater")$pvalvec
#    pvecgrea <- 1 - pvecless
    mipless <- mipvalue(pvecless)
#    mipgrea <- mipvalue(pvecgrea)
    mipgrea <- 1 - mipless
    return(list(mipless=mipless,mipgrea=mipgrea))
}
