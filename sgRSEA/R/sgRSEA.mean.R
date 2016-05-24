sgRSEA.mean <-
function(dat, multiplier=50, r.seed=NULL){

    v1 =as.character(dat[,1])
    v2 = as.character(dat[,2])
    v3 = as.numeric(dat[,3])
    v4 = as.numeric(dat[,4])
    cname=colnames(dat)
    dat = data.frame( v1,v2,v3, v4 )
    colnames(dat)=cname
    if (anyNA(dat)) stop("There is NA in dat.")
    if (min(dat[,3:4]) <1) stop("Minimum counts should be greater than or equal to 1.")

    dat = dat[ order(dat[,2]), ]

    set.seed(r.seed)
    Null.dat = dat[,3:4]

    p0 = pMME(Null.dat)
    Z0 = apply(Null.dat, MARGIN=1, Zstat, p.null=p0)

    ####    (1) FITTING

    Zvec <- Zfit(dat, p.null=p0)
    geneZ = data.frame(dat[,2], Zvec)
    Txx = T.mean( geneZ )
    datZ = outdatZ(dat, Zvec, Txx)
    datZ.s = datZ[ order(datZ$Z, decreasing=T),]

    ####    (2) Generating Null T values

    Null.Tss = Tmean.Null.Tss2(gnames=dat[,2], null.Z=Z0,  multiplier=multiplier)

    meanvec = tapply(Null.Tss[,1], Null.Tss[,2], mean)
    sdvec = tapply(Null.Tss[,1], Null.Tss[,2], sd)

    meanarr = cbind(meanvec,  as.numeric( names(meanvec) ) )
    sdarr = cbind(sdvec,  as.numeric( names(sdvec) ) )

    stdTmat = t( apply(Txx, MARGIN=1, stdofT, meanarr, sdarr) )
    stdNullTmat = t( apply(Null.Tss, MARGIN=1, stdofT, meanarr, sdarr) )

    ####    (3) Significance

    stdNullTmat = stdNullTmat[ order(stdNullTmat[,1], decreasing=T),]
    stdTmat = stdTmat[ order(stdTmat[,1], decreasing=T),]

    null.stdT = stdNullTmat[ ,1]
    obs.stdT = stdTmat[,1]

    Txx.pmat = PermPval(obs.stdT, null.stdT)
    Txx.pmat = cbind( m=stdTmat[,2], Txx.pmat)
    NES = Txx.pmat[,2]

    Txx.p.pos = Txx.pmat[,3]
    FDR.pos = p.adjust(Txx.p.pos, 'BH')
    rank.pos = rank( -Txx.pmat[,2])

    Txx.p.neg = Txx.pmat[,4]
    FDR.neg = p.adjust(Txx.p.neg, 'BH')
    rank.neg = rank( Txx.pmat[,2])

    Txx.pqrmat.pos = cbind(Txx.pmat[,1:2], p.value.pos=Txx.p.pos, FDR.pos,
    rank.pos)[order(NES, decreasing=T),]
    Txx.pqrmat.neg = cbind(Txx.pmat[,1:2], p.value.neg=Txx.p.neg, FDR.neg,
    rank.neg)[order(NES),]

    result = list(gene.pos=Txx.pqrmat.pos, gene.neg=Txx.pqrmat.neg, stdTmat=stdTmat,
    stdNullTmat=stdNullTmat,
                    Tmat=Txx, NullTmat=Null.Tss, sgRNA.stat=datZ.s)
    return(result)
    }
