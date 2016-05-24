seqgen.missing <- function(seqdata, p.cases=.1, p.left=.2, p.gaps=0, p.right=.3,
          mt.left="nr", mt.gaps="nr", mt.right="nr"){
          
    n <- nrow(seqdata)
    lgth <- max(seqlength(seqdata))
    nr.l <- attr(seqdata, mt.left)
    nr.g <- attr(seqdata, mt.gaps)
    nr.r <- attr(seqdata, mt.right)
    
    nm <- round(p.cases * n, 0)
    ## selecting cases
    idm <- sort(sample(1:n, nm))
    rdu.l <- runif(n,min=0,max=p.right)
    rdu.g <- runif(n,min=0,max=p.gaps)
    rdu.r <- runif(n,min=0,max=p.left)
    
    for (i in idm){
        # gaps
        gaps <- sample(1:lgth, round(rdu.g[i] * lgth, 0))
        seqdata[i,gaps] <- nr.g
        # left missing
        nl <- round(rdu.l[i] * lgth, 0)
        if (nl>0) seqdata[i,1:nl] <- nr.l
        # right missing
        nr <- round(rdu.r[i] * lgth, 0)
        if (nr>0) seqdata[i,(lgth-nr):lgth] <- nr.r
    }
    
	return(seqdata)
}


