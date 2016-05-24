"spa.control"<-
  function (eps=1e-6,maxiter=20,gcv=c("lGCV","tGCV","fGCV","aGCV"),
            lqmax=0.2,lqmin=0.05,ldepth=10,ltmin=0.05,lgrid=NULL,
            lval=NULL,dissimilar=TRUE,pce=FALSE,adjust=0,warn=FALSE,...){
    if (eps <= 0) {
        warning("The value of epsilon supplied is <=0; the value 1e-6 was used instead")
        eps=1e-6
    }
    if (maxiter<=0) {
        warning("The value of maxiter supplied is <=0; the value 20 was used instead")
        maxiter=20
    }
    if(ldepth<3){
      warning("The value of ldepth supplied is <=0; the value 3 was used instead")
      ldepth=3
    }
    if(lqmax<lqmin){
      warning("The value of lqmin<lqmax, switching them")
      tmp=lqmax
      lqmax=lqmin
      lqmin=tmp
    }
    if(lqmax>1 || lqmax<0){
      warning("The value of lqmax must be between [0,1], default=0.95")
      lqmax=0.95
    }
    if(lqmin>1 || lqmin<0){
      warning("The value of lqmax must be between [0,1], default=0.05")
      lqmin=0.05
    }
    
    if(!is.null(lgrid)){
      if(length(lgrid)>1)
        lgrid=lgrid[1]
      if(lgrid<0)
        lgrid=-lgrid
      lgrid=ceiling(lgrid)
    }
    if(missing(gcv))
      gcv="aGCV"
    maxiter=round(maxiter,0)
    list(eps = eps, maxiter = maxiter, gcv=gcv,dissimilar = dissimilar,
         lqmax=lqmax,lqmin=lqmin,ldepth=ldepth,lgrid=lgrid,ltmin=ltmin,
         lval=lval,warn = warn,pce=pce,adjust = adjust)
}
