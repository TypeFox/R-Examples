##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 
##
##     Written by Karim Rahim and Wesley Burr.
## 
##     This file is part of the multitaper package for R.
##     http://cran.r-project.org/web/packages/multitaper/index.html
## 
##     The multitaper package is free software: you can redistribute it and 
##     or modify it under the terms of the GNU General Public License as 
##     published by the Free Software Foundation, either version 2 of the 
##     License, or any later version.
##
##     The multitaper package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
##
##     If you wish to report bugs please contact the author:
## 
##     Karim Rahim
##     karim.rahim@gmail.com


.mweave <-  function (x,dw,swz,ndata,nord,ssqswz,dt_) {
    
    out <- .Fortran("mweave", as.double(x), as.double(dw),
                    as.double(swz), as.integer(ndata),
                    as.integer(nord), as.double(ssqswz),
                    cntr=double(1), as.double(dt_),
                    spz=double(1), varc=double(1),
                    PACKAGE='multitaper')
    return(list(cntr=out$cntr, spz=out$spz, varc=out$varc))
}

.HF4mp1 <- function(cft, swz, nord, ssqswz) {

    ## ######################################
    ## The notation and function names were chosen
    ## to map to original fortran (f77) code.
    ## Note to obtain:  swz <- apply(dw, 2, sum)
    ## swz is the zeroth frequency Fourier transform of the
    ## Slepian sequences. It is H_k(0) from P ercival and Walden (1993)
    ## pages 497--399.
    ## (just to define dw) dw <- dpssIN$v*sqrt(deltaT)    
    ## Vectorized from original F77 code
    ## Equation (13.5) of:
    ##   Thomson, D.J. Spectrum Estimation and Harmonic Analysis,
    ##   Proceedings of the IEEE, 1982.

    cmv <- (cft %*% swz) /ssqswz
    ssqave <-  (Mod(cmv)^2)*ssqswz
    swz <- as.matrix(swz)
    
    ssqres <- apply( Mod(cft - (cmv %*% t(swz)))^2,
                    1, sum)
    F_<- (nord-1)*ssqave/ssqres
    
    return(list(Ftest=F_,cmv=cmv))
}


.mw2wta <- function(sa, nfreq, nord,
                    var, dt_, ev, evp=(1-ev),
                    tol=.03, maxadaptiveiteration=100) {

    ## this is equation (5.3) and (5.4) form 
    ##   Thomson, D.J. Spectrum Estimation and Harmonic Analysis,
    ##   Proceedings of the IEEE, 1982.
    ## note that the weights are squared, they are |d_k(f)^2 from equation
    ## (5.4)

    out <- .Fortran("mw2wta", as.double(sa),
                    wt=matrix(as.double(0), nfreq, nord),
                    as.integer(nfreq), as.integer(nord),
                    s=double(nfreq), as.double(ev), as.double(evp),
                    dofs=double(nfreq), dofav=double(1),
                    as.double(var), as.double(dt_),
                    as.double(tol),
                    as.integer(maxadaptiveiteration),
                    mxiter=integer(1), aviter=double(1),
                    PACKAGE='multitaper')
    
    return(list(s=out$s, wt=out$wt, dofs=out$dofs, dofav=out$dofav,
                mxiter=out$mxiter, aviter=out$aviter))
}

.mw2jkw <-  function(sa, nfreq, nord, var, dt_, ev,
                     evp=(1-ev), tol=.03,
                     maxadaptiveiteration=100) {
    
    nordP2 <-  nord+2
    out <- .Fortran("mw2jkw", as.double(sa),
                    wt=matrix(as.double(0), nfreq, nord),
                    as.integer(nfreq), as.integer(nord),
                    s=double(nfreq), as.double(ev),
                    as.double(evp), dofs=double(nfreq),
                    dofav=double(1), as.double(var),
                    as.double(dt_), as.double(tol),
                    sjk=double(nordP2), varjk=double(nfreq),
                    bcjk=double(nfreq),
                    matrix(as.double(0), nord, nordP2),
                    double(nordP2), double(nord),
                    as.integer(maxadaptiveiteration),
                    mxiter=integer(1),
                    PACKAGE='multitaper')
                   
    
    return(list(s=out$s, wt=out$wt, dofs=out$dofs, dofav=out$dofav,
                mxiter=out$mxiter,
                varjk=out$varjk, bcjk=out$bcjk, sjk=out$sjk))
}

.qsF <- function(nFreqs,nFFT,k,cft,useAdapt,kadapt) {

    out <- .Fortran("quickSineF", as.integer(nFreqs),
                    as.integer(nFFT), as.integer(k),
                    cft=cft, as.logical(useAdapt),
                    kadapt=matrix(data=as.double(kadapt),nrow=nFreqs,ncol=1),
                    spec=matrix(data=double(nFreqs),nrow=nFreqs,ncol=1),
                    PACKAGE='multitaper')

    return(list(spec=out$spec))
}

.cF <- function(n,v) {

  out <- .Fortran("curbF",as.integer(n),as.double(v),PACKAGE='multitaper')
  opt <- out[[2]]
  return(list(opt=opt))
}

.nF <- function(n,i1,i2,s) {

   out <- .Fortran("northF",as.integer(n),as.integer(i1),
                   as.integer(i2),sx=matrix(data=as.double(s),nrow=n,ncol=1),
                   ds=double(1), dds=double(1),
                   PACKAGE='multitaper')

   return(list(ds=out$ds, dds=out$dds))
}

.adaptSine <- function(ntimes, k, nFreqs, sx, nFFT, cft, df, fact) {

    out <- .Fortran("adapt",as.integer(ntimes),as.integer(k),
                    as.integer(nFreqs),sx=matrix(data=as.double(sx),nrow=nFreqs,ncol=1),
                    as.integer(nFFT), cft=cft, as.double(df),
                    kopt=double(nFreqs),fact=as.double(fact),PACKAGE='multitaper')

    return(list(spec=out$sx,kadapt=out$kopt))

}

