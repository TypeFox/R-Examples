
 GSA.listsets=
function (GSA.obj, geneset.names = NULL, maxchar = 20, FDRcut = 0.2) 
{
if(is.null(geneset.names)){
  geneset.names=rep("xxxxxx", length(GSA.obj$GSA.scores))
}
negflag= !(GSA.obj$resp.type=="Multiclass" | GSA.obj$resp.type=="taCorr")
    r = GSA.obj$GSA.scores
    rstar = GSA.obj$GSA.scores.perm
    r[is.na(r)] = 0
    rstar[is.na(rstar)] = 0
    nperms = ncol(GSA.obj$GSA.scores.perm)
    np = length(r)
geneset.names = substring(geneset.names, 1, maxchar)
    pvalues.lo = GSA.obj$pvalues.lo
    pvalues.hi = GSA.obj$pvalues.hi
    m = sum(!is.na(pvalues.hi))

    make.monotone.increasing = function(x) {
        n = length(x)
        if (n==0) return(NULL)  ## added to prevent error when x is NULL
        if (n==1) x[1] <- x[1]  ## added to prevent error when length of x =1
        else {
   	      for (i in n:2) {
         	  if (x[i - 1] > x[i]) {
           	    x[i - 1] = x[i]
           	}
		      }
        }
        return(x)
    }



    oo = (1:length(r))[!is.na(pvalues.hi)]
    fdr.lo = fdr.hi = rep(NA, length(r))
    for (i in oo) {
    if(negflag){ fdr.lo[i] = round(m * pvalues.lo[i]/sum(pvalues.lo[!is.na(pvalues.lo)] <=
            pvalues.lo[i]), 4) }
        fdr.hi[i] = round(m * pvalues.hi[i]/sum(pvalues.hi[!is.na(pvalues.hi)] <= 
            pvalues.hi[i]), 4)
    }

 fdr.lo=pmin(fdr.lo,1)
 fdr.hi=pmax(fdr.hi,0)

    res1=NULL
    if(negflag){
        oo1 = (1:length(r))[r < 0]
        res1 = NULL
        for (i in oo1) {
                res1 = rbind(res1, c(i, geneset.names[i], round(GSA.obj$GSA.scores[i], 
                    4), pvalues.lo[i], fdr.lo[i]))
        }
        if (!is.null(res1)) o1 = order(res1[, 4],decreasing = FALSE)  
        else o1 = NULL

        res1 = res1[o1, , drop = F]
        res1[, 5] = make.monotone.increasing(as.numeric(res1[, 5]))
   }

    oo2 = (1:length(r))[r > 0]
    res2 = NULL
    for (i in oo2) {
            res2 = rbind(res2, c(i, geneset.names[i], round(GSA.obj$GSA.scores[i], 
                4), pvalues.hi[i], fdr.hi[i]))
    }
    if (!is.null(res2)) o2 = order(res2[, 4], decreasing = FALSE)  
    else o2 = NULL

    res2 = res2[o2, , drop = F]
    res2[, 5] = make.monotone.increasing(as.numeric(res2[, 5]))
    if (length(res1) == 0) {
        res1 = NULL
    }
    if (length(res2) == 0) {
        res2 = NULL
    }
    if ( (length(res1) > 0) & negflag) {
            dimnames(res1) = list(NULL, c("Gene_set", "Gene_set_name", 
                "Score", "p-value", "FDR"))
    }
    if (length(res2) > 0) {
            dimnames(res2) = list(NULL, c("Gene_set", "Gene_set_name", 
                "Score", "p-value", "FDR"))
    }
nsets.neg=NULL
   if(negflag){ res1 = res1[as.numeric(res1[, 5]) <= FDRcut,,drop=FALSE ]
        nsets.neg = nrow(res1)
      if (is.null(res1)) {
        nsets.neg = 0
    }}
    res2 = res2[as.numeric(res2[, 5]) <= FDRcut,,drop=FALSE ]
    nsets.pos = nrow(res2)
    if (is.null(res2)) {
        nsets.pos = 0
    }
    return(list(FDRcut = FDRcut, negative = res1, positive = res2, 
        nsets.neg = nsets.neg, nsets.pos = nsets.pos))
}

