glassopath=function (s, rholist=NULL, thr = 1e-04, maxit = 10000, approx = FALSE
, 
                penalize.diagonal = TRUE,  w.init = NULL, 
                wi.init = NULL, trace = 1) 
{
        n = nrow(s)
if(is.null(rholist)){
      rholist=seq(max(abs(s))/10,max(abs(s)),length=10)
  }
rholist=sort(rholist)
                ipen = 1 * (penalize.diagonal)
                ia = 1 * approx
                rho=xx=ww=matrix(0,n,n)
                nrho=length(rholist)
                beta=what=array(0,c(n,n,nrho))
                jerrs=rep(0,nrho)
                mode(rholist) = "double"
                mode(nrho) = "integer"
                mode(rho) = "double"
                mode(s) = "double"
                mode(ww) = "double"
                mode(xx) = "double"
                mode(n) = "integer"
                mode(maxit) = "integer"
                mode(ia) = "integer"
                mode(trace) = "integer"
                mode(ipen) = "integer"
                mode(thr) = "double"
                mode(beta) = "double"
                mode(what) = "double"
                mode(jerrs) = "integer"

        junk <- .Fortran("glassopath", beta=beta,what=what,jerrs=jerrs,rholist,nrho,n, s, rho, ia,   trace, ipen,
        thr, maxit = maxit, ww = ww, xx = xx, niter = integer(1), 
        del = double(1), ierr = integer(1), PACKAGE="glasso")

    xx = array(junk$beta,  c(n,n,nrho))
    what = array(junk$what,  c(n,n,nrho))
    return(list(w=what, wi = xx, 
        approx = approx, rholist=rholist, errflag = junk$jerrs))
}

