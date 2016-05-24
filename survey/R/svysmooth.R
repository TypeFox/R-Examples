svysmooth<-function(formula,design,...) UseMethod("svysmooth", design)
svysmooth.default<-function(formula, design,method=c("locpoly","quantreg"),bandwidth=NULL,quantile,df=4,...){
  switch(match.arg(method),
         locpoly=svylocpoly(formula,design,bandwidth=bandwidth,...),
         quantreg=svyrqss(formula,design,quantile=quantile,df=df,...)
         )
}

fitted.rq<-function(object,...) object$x%*% object$coefficients/object$weights

svyrqss<-function(formula,design,quantile=0.5,df=4,...){
  require("quantreg") || stop("needs quantreg package")
  require("splines") || stop("needs splines package, which should be part of R")

  mf<-model.frame(formula, model.frame(design), na.action=na.omit)
  naa<-attr(mf,"na.action")

  tt<-attr(terms(formula),"term.labels")
  df<-rep(df, length=length(tt))
  quantile<-rep(quantile, length=length(tt))

  if (length(formula)==3){
    density<-FALSE
  } else {
    density<-TRUE
    stop("type='quantreg' does not do densities")
  }
  
  w<-weights(design,type="sampling")
  if (length(naa)) w<-w[-naa]
  environment(formula)<-environment()

  ll<-vector("list", length(tt))
  for(i in 1:length(tt)){
    termi<-as.name(tt[i])
    ff<-eval(bquote(update(formula,.~bs(.(termi),df=.(df[i])))))
    rqfit<-rq(ff, tau=quantile[i],weights=w,data=mf,...)
    xx<-mf[,i+1]
    oo<-order(xx)
    ll[[i]]<-list(x=xx[oo],y=fitted.rq(rqfit)[oo])
  }
  names(ll)<-attr(terms(formula),"term.labels")
  
  attr(ll,"call")<-sys.call(-2)
  attr(ll,"density")<-density
  if(density)
    attr(ll,"ylab")<-"Density"
  else
    attr(ll,"ylab")<-deparse(formula[[2]])

  class(ll)<-"svysmooth"
  
  ll

}
  
svylocpoly<-function(formula, design, ngrid=401, xlim=NULL,
                     ylim=NULL, bandwidth=NULL,...){
  require("KernSmooth") || stop("needs KernSmooth package")

  mf<-model.frame(formula,model.frame(design))
  mm<-model.matrix(terms(formula),mf)
  if(attr(terms(formula),"intercept"))
    mm<-mm[,-1,drop=FALSE]

  naa<-attr(mf,"na.action")


  if (length(formula)==3){
    Y<-model.response(mf)
    density<-FALSE
  } else density<-TRUE
  

  if (is.null(xlim)){
    xlim<-apply(mm,2,range)
  }
  if (!is.matrix(xlim))
    xlim<-matrix(xlim,nrow=2)


  if (is.null(bandwidth)){
    bandwidth<-numeric(ncol(mm))
    for(i in 1:ncol(mm)){
      bandwidth[i]<-if(density) dpik(mm[,i],gridsize=ngrid) else dpill(mm[,i],Y,gridsize=ngrid)
    }
  } else {
    bandwidth<-rep(bandwidth, length=ncol(mm))
  }
  
  w<-weights(design,type="sampling")
  if (length(naa)) w<-w[-naa]
  
  ll<-vector("list", ncol(mm))
  for(i in 1:NCOL(mm)){
    gx<-seq(min(xlim[,i]), max(xlim[,i]), length=ngrid)
    nx<-rowsum(c(rep(0,ngrid),w), c(1:ngrid, findInterval(mm[,i],gx)))
    if (density){
      ll[[i]]<-locpoly(rep(1,ngrid),nx*ngrid/(diff(xlim[,i])*sum(w)),
                           binned=TRUE, bandwidth=bandwidth[i], range.x=xlim[,i])
    }else{
      ny<-rowsum(c(rep(0,ngrid), Y*w), c(1:ngrid, findInterval(mm[,i],gx)))
      ll[[i]]<-locpoly(nx, ny, binned=TRUE, bandwidth=bandwidth[i], range.x=xlim[,i])
    }
    names(ll)<-attr(terms(formula),"term.labels")
  }
  attr(ll,"call")<-sys.call(-2)
  attr(ll,"density")<-density
  if(density)
    attr(ll,"ylab")<-"Density"
  else
    attr(ll,"ylab")<-deparse(formula[[2]])

  class(ll)<-"svysmooth"
  
  ll
 
}

print.svysmooth<-function(x,...){
  if(attr(x,"density"))
    cat("Density estimate: :")
  else
    cat("Scatterplot smoother :")
  print(attr(x,"call"))
  invisible(x)
}



plot.svysmooth<-function(x, which=NULL,type="l",xlabs=NULL,ylab=NULL,...){
  if (is.null(which))
    which<-seq(length=length(x))
  if (is.character(which))
    which<-match(which,names(x))
  
  if(is.null(xlabs)) xlabs<-names(x)[which]
  if(is.null(ylab)) ylab<-attr(x,"ylab")

  for(i in seq(length=length(which)))
    plot(x[[which[i]]], type=type, xlab=xlabs[i], ylab=ylab, ...)
  invisible(NULL)
}

lines.svysmooth<-function(x,which=NULL,...){
  for(i in names(x)) lines(x[[i]],...)
}
