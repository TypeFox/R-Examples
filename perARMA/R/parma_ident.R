parma_ident <-
function(x,T,missval,datastr,...)
{
  parma_ident <- function(x,T,missval,datastr,outdir,details)
   {
       x=as.matrix(x)

      outdir=as.character(outdir)
      unlink(outdir, recursive=TRUE)
      sink("ident.txt")
      cat(paste('writing text output to',outdir,'\n'))

      BTAU=seq(0,16)
      BMETH=0
      MAXTAU=15
      orig_ind=seq(1,2*MAXTAU,2)
      adjust_ind=seq(2,2*MAXTAU,2)
      p_perpacf=8
      ALPHA=.01

     dirname=outdir

     dev.set(which=1)  
     pm<-permest(x,T,ALPHA,missval,datastr,pp=details)
     pmean=pm$pmean
     xd=pm$xd

     dev.set(which=1)  
     ps<-persigest(x,T,ALPHA,missval,datastr,pp=details)
     pstd=ps$pstd
     xn=ps$xn


     cat(paste('####################################','\n'))
     cat(paste('Bartlett test for equal sigmas T=',T,' pv=',ps$pspv,'\n'))
     cat(paste('####################################','\n'))
     cat(paste(' B(t,tau) = R(t+tau,t) from peracf.R','\n'))

     if (details)
     { PRTTAU=seq(0,MAXTAU)
       } else {
       PRTTAU=NULL}

    
     peracf_out<-peracf(x,T,seq(0,MAXTAU),missval,datastr,prttaus=PRTTAU,plottaus=BTAU)
     B=peracf_out$B
     B=as.matrix(B)

    cat(paste('####################################','\n'))
    cat(paste('B_k(tau) from Bcoeff','\n')) 

    Bc<-Bcoeff(t(x),T,BTAU,missval,datastr,printflg=details,meth=BMETH)
    pvecB=Bc$pv
    cat(paste('P-values for test of |B_k(tau)|^2=0 using neighbors of FFT bins; uses F distn','\n'))
    cat(paste('rows are index k=0:T-1 , cols are tau=0:maxtau','\n'))
    cat(paste('first line, k=0, is time average B_0(tau)=(1/T)sum_t B(t,tau);','\n',' first col is B_k(0), FT of sigma^2(t)','\n'))
    cat(paste('number of rows (k) fewer than T due to real','\n'))

    nc_pvecB=ncol(pvecB)
    nr_pvecB=nrow(pvecB)

    min1=matrix(0,nr_pvecB,1)
    bfpv_nztau=matrix(0,nr_pvecB,1)
    for (i in 1:nr_pvecB)
      { min1[i]=min(pvecB[i,2:nc_pvecB])}
        min1m=min1*(nc_pvecB-1)
        bfpv_nztau=min1m
        for(i in 1:nr_pvecB){
            if (bfpv_nztau[i]>=1) {bfpv_nztau[i]=1}
            }

    print(min1m)
    print(bfpv_nztau)
    cat(paste('last col is row pv  = min pv over non zero lags, bf corrected for', nc_pvecB-1,' lags','\n'))

    cat(paste('tau=','\n'))
    print(BTAU)
    cat(paste('\n'))
    print(c(pvecB,bfpv_nztau))

    bfpv_nonconst=min(pvecB[2:nr_pvecB,])*(nr_pvecB-1)
    cat(paste('col pv = min pv over non zero freqs, bf corrected for', nr_pvecB-1,'values of k','\n'))
     if (bfpv_nonconst>1) {bfpv_nonconst=1}
    cat(paste(bfpv_nonconst,'\n'))
    cat(paste('B_k(0)=0,k=1,2,... equiv to equal sigmas','\n'))
    cat(paste('\n'))
    ind_nztau=which(BTAU!=0);
     pvecB_nonconst_nztau=pvecB[2:nr_pvecB,ind_nztau]
    N_pvecB_nonconst_nztau=length(pvecB_nonconst_nztau[,])
    bfpv_nonconst_nztau=min(pvecB_nonconst_nztau[,]*N_pvecB_nonconst_nztau)
    cat(paste(bfpv_nonconst_nztau,'\n'))
     if (bfpv_nonconst_nztau>=1) {bfpv_nonconst_nztau=1}  
    cat(paste('min pv over', nr_pvecB-1,'non zero freqs and',length(ind_nztau),'nonzero lags, bf corrected: ',bfpv_nonconst_nztau, '\n'))
   
    cat(paste('####################################','\n'))
    cat(paste('Bcoeffa = phth2ab(B_k(tau))','\n'))
    cat(paste('\n'))
    Bca<-Bcoeffa(t(x),T,BTAU,missval,datastr,printflg=details,meth=BMETH)
    pvec=Bca$pvec
    cat(paste('rows are abcoeff index k=0:T-1 , cols are tau=0:maxtau','\n'))
    cat(paste('B(t,tau) = sum from k=1 to  floor(T/2) [a(1)(tau) + sum a(2k)(tau)*cos(2*pi*k*n/T)+a(2k+1)(tau)*sin(2*pi*k*n/T)]','\n'))
    cat(paste('first line, a1, is time average (1/T)sum_t B(t,tau)','\n'))
    cat(paste('next rows are a2,a3,a4,...','\n'))

    cat(paste('tau=','\n'))
    print(BTAU)
    cat(paste('\n'))
    print(pvec)

    nc_pvec=ncol(pvec)
    nr_pvec=nrow(pvec)
    bfpv_nonconst=min(pvec[2:nr_pvec,])*(T-1)
    cat(paste('min pv over non zero freqs, bf corrected','\n'))
    if (bfpv_nonconst>1) {bfpv_nonconst=1}
    cat(bfpv_nonconst)

    ind_nztau=which(BTAU !=0)
    pvec_nonconst_nztau=pvec[2:nr_pvec,ind_nztau] 
    N_pvec_nonconst_nztau=length(pvec_nonconst_nztau[,])
    bfpv_nonconst_nztau=min(pvec_nonconst_nztau[,]*N_pvec_nonconst_nztau)
      if (bfpv_nonconst_nztau>=1) {bfpv_nonconst_nztau=1}
    cat(paste('min pv over', nr_pvecB-1,'non zero freqs and',length(ind_nztau),'nonzero lags, bf corrected: ',bfpv_nonconst_nztau, '\n'))

    
    cat(paste('####################################','\n'))
    cat(paste('perpacf  pi(t,n)','\n'))
    cat(paste('following does not combine individual p-values but tests t-mean to be zero ','\n'))
    cat(paste('doing as chi2 or combining individual p-values may be better','\n'))

    perpacf_out<-perpacf(x,T,p_perpacf,missval)
    ppa=perpacf_out$ppa
    nsamp=perpacf_out$nsamp

    dev.set(which=1) 
    ppfplot(ppa,min(nsamp[,]),.05,datastr)

    cat(paste('####################################','\n'))
    cat(paste('ppfcoeffab','\n'))

    print(ppa)
    ppfcoeffab_out<-ppfcoeffab(ppa,nsamp,details,datastr)
    pkab=ppfcoeffab_out$pkab
    pkpv=ppfcoeffab_out$pkpv
    cat(paste('rows are abcoeff index k=0:T-1 , cols are n=0:p_perpacf','\n'))
    cat(paste('first line, k=0, is time average (1/T)sum_t pi(t,n)','\n'))
    print(seq(0,p_perpacf))
    cat(paste('\n'))
    ncpkpv=ncol(pkpv)
    nrpkpv=nrow(pkpv)

    min1=matrix(0,nrpkpv,1)
      for (i in 1:nrpkpv) { min1[i]=min(pkpv[i,2:ncpkpv])}
    min1m=min1*(ncpkpv-1)
      for(i in 1:nrpkpv)  { if (min1m[i]>1) {min1m[i]=1}
    }
     bfpv_nztau=min1m
     cat(paste('last col is row pv  = min pv over non zero lags, bf corrected for', ncpkpv-1,' lags','\n'))

     print(c(pkpv,bfpv_nztau))

     bfpv_nonconst=min(pkpv[2:nrpkpv,])*(T-1)
     cat(paste('min pv over non zero freqs, bf corrected','\n'))
     if (bfpv_nonconst>=1) {bfpv_nonconst=1}  
      print(bfpv_nonconst)

     pvec_nonconst_nztau=pkpv[2:nrpkpv,]
     N_pvec_nonconst_nztau=length(pvec_nonconst_nztau[,])
     bfpv_nonconst_nztau=min(pvec_nonconst_nztau[,]*N_pvec_nonconst_nztau)
        if (bfpv_nonconst_nztau>=1) {bfpv_nonconst_nztau=1}
     cat(paste('min pv over', T-1,'non zero freqs and',ncol(pvec),'nonzero lags, bf corrected: ',bfpv_nonconst_nztau, '\n'))
 
     dev.set(which=1) 
     acfpacf(t(x),16,p_perpacf,datastr)

      sink()
      dir.create(outdir)
      file.copy("ident.txt",outdir)
      unlink("ident.txt", recursive=TRUE)

    result = list(pmean=pmean,xd=xd,pstd=pstd,xn=xn)
   class(result) = "parma_ident"
    result
}

L<-modifyList(list(outdir='IDENT_OUT',details=1),list(x = x, T=T, missval=missval, datastr=datastr, ...))

 do.call(parma_ident,L)

}