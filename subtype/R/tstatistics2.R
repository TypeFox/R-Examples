tstatistics2 <-
function (xdat, grp, logse = FALSE, paired = FALSE) 
{
    glab = unique(grp)
    if (length(glab) != 2 ) {
        #stop("Need exactly two different groups")
      ret = rep(0,dim(xdat)[1]);
      data.frame(tstat = ret);
    }
    if (!paired && length(glab)== 2) {
        ndx = grp == glab[1]
        n1 = table(ndx)["TRUE"]
        n2 = table(ndx)["FALSE"]
        if (sum(ndx)>1) {m1 = rowMeans(xdat[, ndx])}
      if (sum(ndx)==1) {m1 = xdat[, ndx]}
        if (sum(!ndx)>1) {m2 = rowMeans(xdat[, !ndx])}
      if (sum(!ndx)==1) {m2 = xdat[, !ndx]} 
        xdat[, ndx] = xdat[, ndx] - m1
        xdat[, !ndx] = xdat[, !ndx] - m2
        pvar = rowSums(xdat * xdat)/(n1 + n2 - 2)
        se = sqrt(pvar * (1/n1 + 1/n2))
    
    se[which(se==0)]<-10^-7

    mn = m1 - m2
        ret = mn/se
    }
    if (paired && length(glab)== 2) {
        n2 = ncol(xdat)
        n = floor(0.5 * n2)
        if (n2/2 != n) {
            stop("Paired t test requires an even number of observations")
        }
        ondx = seq(1, n2, by = 2)
        endx = seq(2, n2, by = 2)
        sodd = (-1)^as.numeric(grp[ondx] != glab[1])
        sevn = (-1)^as.numeric(grp[endx] != glab[2])
        if (!all(sodd == sevn)) {
            stop("paired observations must be consecutive")
        }
        d = (xdat[, ondx] - xdat[, endx]) %*% diag(sodd)
        m = rowMeans(d)
        d = d - m
        pvar = rowSums(d * d)/(n - 1)
        se = sqrt(pvar/n)
        mn = m
        ret = m/se
    }
    if (logse) {
        data.frame(tstat = ret, logse = log(se))
    }
    else {
        data.frame(tstat = ret)
    }
}
