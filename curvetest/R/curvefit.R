curvefit <-
function(formula, data,  kernel="Quartic", alpha=0.5,  bw=NULL, 
      myx, bcorrect="simple", getit=F){    
     if(is.character(kernel)) {
      kernel.char=kernel
      kernel=get.weight.function(kernel)
     }else{kernel.char=attr(kernel, "name")}
    if(length(bcorrect)>1) bcorrect=bcorrect[1] 
    x<-model.frame(formula, data)[,2] 
    y<-model.frame(formula, data)[,1] 
    if(missing(myx)) myx=x else if(length(myx)==1) myx=seq(min(x), max(x), length=myx)
    if(!is.numeric(alpha)|is.null(alpha)) alpha=getoptimalalpha(formula, data)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data",  "kernel", "alpha", "myx", "plotit"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf[1:3], parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf)[,2]
    xlab<- deparse(formula[[3]])
    ylab<- deparse(formula[[2]])
    #
    o<-getWV(x=x, myx=myx, kernel=kernel.char, alpha=alpha, bw=bw, bc=bcorrect, getit=FALSE)
    x<-o$x;  myx<-o$myx; nn<-length(myx); n1<-length(x);  k0<-o$k0; 
    delta<-o$delta; delta2<-o$delta2;vv<-o$vv;TX<-o$TX; sx<-o$sx;
    LL<-o$LL; AA<-o$AA;
    ##############################0######    
    ey <-  LL%*% y
    bic<-log(sum(ey^2))*n1+log(n1)*delta/2 ##BIC statistic
    out<-list(formula=formula, data.origninal=data, data.model=data.frame(y=y, x=x), 
      xlab=xlab, ylab=ylab, fits=ey, 
      res=y-ey,  myx=myx, fitted=sx%*%y, BIC=bic,  k0=k0, vv=vv, 
      alpha=alpha, sx=sx,AA=AA, kernel=kernel, LL=LL, delta=delta,  delta2=delta2, bcorrect=bcorrect)
    class(out)<-"curvefit"    
    return(invisible(out))
 }
