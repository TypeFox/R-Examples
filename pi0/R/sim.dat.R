sim.dat=function(G=1e4,pi0=.75,gamma2=1,n1=5,n2=n1,
                errdist=rnorm,effdist=function(g,gamma2)rnorm(g,,sqrt(gamma2)),
                ErrArgs,EffArgs
){
    G1=round(G*(1-pi0))
    N=n1+n2
    dat=matrix({if(missing(ErrArgs)) errdist(G*N) else errdist(G*N,ErrArgs)},G,N)
    eff=if(missing(EffArgs)) effdist(G1,gamma2) else effdist(G1,gamma2,EffArgs)
    dat[1:G1,1:n1]=dat[1:G1,1:n1]+eff
    attr(dat,'G')=G
    attr(dat,'G1')=G1
    attr(dat,'n1')=n1
    attr(dat,'n2')=n2
    attr(dat,'gamma2')=if(is.null(attr(eff,'gamma2')))gamma2 else attr(eff,'gamma2')
    dat
}

