
svykm<-function(formula, design, se=FALSE, ...) UseMethod("svykm",design)

svykm.survey.design<-function(formula, design,se=FALSE, ...){
  require("survival")
  if (!inherits(formula,"formula")) stop("need a formula")
  if (length(formula)!=3) stop("need a two-sided formula")
  mf<-model.frame(formula, model.frame(design), na.action=na.pass)
  mf<-na.omit(mf)
  drop<-attr(mf,"na.action")
  if (!is.null(drop)) 
    design<-design[-drop,] 
  y<-model.response(mf)
  if (!is.Surv(y) || attr(y,"type")!="right")
    stop("response must be a right-censored Surv object")      

  if (ncol(mf)==1) {
    if (se)
      s<-km.stderr(y,design)
    else 
      s<-svykm.fit(y,weights(design))
  } else {
    x<-mf[,-1]
    if (NCOL(x)>1)
      groups<-do.call(interaction,x)
    else
      groups<-as.factor(x)
    
    if (se){
      lhs<-formula
      lhs[[3]]<-1
      s<-lapply(levels(groups), function(g) svykm(lhs, subset(design,groups==g),se=TRUE))
    }else{
      s<-lapply(levels(groups), function(g) svykm.fit(y[groups==g],weights(design)[groups==g]))
    }
    names(s)<-levels(groups)
    class(s)<-"svykmlist"
  }
  call<-match.call()
  call[[1]]<-as.name(.Generic)
  attr(s,"call")<-call
  attr(s, "formula")<-formula
  attr(s, "na.action")<-drop
  return(s)
}


svykm.svyrep.design<-function(formula, design,se=FALSE, ...){
  require("survival")
  if (!inherits(formula,"formula")) stop("need a formula")
  if (length(formula)!=3) stop("need a two-sided formula")
  mf<-model.frame(formula, model.frame(design), na.action=na.pass)
  mf<-na.omit(mf)
  drop<-attr(mf,"na.action")
  if (!is.null(drop)) 
    design<-design[-drop,] 
  y<-model.response(mf)
  if (!is.Surv(y) || attr(y,"type")!="right")
    stop("response must be a right-censored Surv object")      

  if (ncol(mf)==1) {
    if (se)
      stop("SE not yet available")
    else 
      s<-svykm.fit(y,weights(design,"sampling"))
  } else {
    x<-mf[,-1]
    if (NCOL(x)>1)
      groups<-do.call(interaction,x)
    else
      groups<-as.factor(x)
    
    if (se){
      lhs<-formula
      lhs[[3]]<-1
      s<-lapply(levels(groups), function(g) svykm(lhs, subset(design,groups==g),se=TRUE))
    }else{
      s<-lapply(levels(groups), function(g) svykm.fit(y[groups==g],weights(design)[groups==g]))
    }
    names(s)<-levels(groups)
    class(s)<-"svykmlist"
  }
  call<-match.call()
  call[[1]]<-as.name(.Generic)
  attr(s,"call")<-call
  attr(s, "formula")<-formula
  attr(s, "na.action")<-drop
  return(s)
}


svykm.fit<-function(y,w){
  t<-y[,"time"]
  s<-y[,"status"]
  nn<-rowsum(cbind(s,1)*w,t)
  tt<-sort(unique(t))
  N<-c(sum(w),sum(w),sum(w)-cumsum(nn[-nrow(nn),2]))
  d<-c(0,nn[,1])
  surv<-pmax(0,cumprod(1-d/N))
  rval<-list(time=c(0,tt), surv=surv)
  class(rval)<-"svykm"
  rval
}

km.stderr<-function(survobj,design){
  time<-survobj[,'time']
  status<-survobj[,'status']
  ## Brute force and ignorance: compute Y and dN as totals, use delta-method
  keep<-which((status==1) & (weights(design)!=0))
  y<-outer(time,time[keep],">=")
  dN<-diag(status)[,keep]
  oo<-order(time[keep], -status[keep])
  okeep<-keep[oo]
  ntimes<-length(oo)
  ttime<-time[okeep]
  sstatus<-status[okeep]

  totals<-svytotal(cbind(dN[,oo],y[,oo]), design)
  rm(dN)
  y<-coef(totals)[-(1:ntimes)]
  dNbar<-coef(totals)[1:ntimes]
  
  h<-cumsum(dNbar/y)
  
  dVn<- vcov(totals)[(1:ntimes),(1:ntimes)]/outer(y,y)
  dVy <- vcov(totals)[-(1:ntimes),-(1:ntimes)]*outer(dNbar/y^2,dNbar/y^2)
  dCVny<- -vcov(totals)[(1:ntimes),-(1:ntimes)]*outer(1/y,dNbar/y^2)
  dV<-dVn+dVy+dCVny+t(dCVny)
  
  V<-numeric(ntimes)
  V[1]<-dV[1,1]
  for(i in 2:ntimes) V[i]<-V[i-1]+sum(dV[1:(i-1),i])+sum(dV[i,1:i])
  
  rval<-list(time=ttime,surv=exp(-h),varlog=V)
  class(rval)<-"svykm"
  rval
}


plot.svykm<-function(x,xlab="time",ylab="Proportion surviving",ylim=c(0,1),ci=NULL,lty=1,...){
  if (is.null(ci))
    ci<-!is.null(x$varlog)
  
  plot(x$time,x$surv,xlab=xlab,ylab=ylab, type="s",ylim=ylim,lty=lty,...)

  if (ci){
    if (is.null(x$varlog))
      warning("No standard errors available in object")
    else{
      lines(x$time,exp(log(x$surv)-1.96*sqrt(x$varlog)),lty=2,type="s",...)
      lines(x$time,pmin(1,exp(log(x$surv)+1.96*sqrt(x$varlog))),lty=2,type="s",...)
    }
  }
  invisible(x)
}


lines.svykm<-function(x,xlab="time",type="s",ci=FALSE,lty=1,...){
  lines(x$time,x$surv, type="s",lty=lty,...)
  if (ci){
    if (is.null(x$varlog))
      warning("no standard errors available in object")
    else {
      lines(x$time,exp(log(x$surv)-1.96*sqrt(x$varlog)),lty=2,type="s",...)
      lines(x$time,pmin(1,exp(log(x$surv)+1.96*sqrt(x$varlog))),lty=2,type="s",...)
    }
  }
  invisible(x)
}

plot.svykmlist<-function(x, pars=NULL, ci=FALSE,...){
  if (!is.null(pars)) pars<-as.data.frame(pars,stringsAsFactors=FALSE)

  if(is.null(pars))
    plot(x[[1]],ci=ci,...)
  else
    do.call(plot,c(list(x[[1]]),pars[1,,drop=FALSE],ci=ci,...))

  m<-length(x)
  if(m==1) return
  for(i in 2:m){
    if(is.null(pars))
      lines(x[[i]],ci=ci,...)
    else
      do.call(lines,c(list(x[[i]]),pars[i,,drop=FALSE],ci=ci,...))
  }
  invisible(x)
}

print.svykm<-function(x, digits=3,...,header=TRUE){
  if (header) {cat("Weighted survival curve: ")
               print(attr(x,"call"))}
  suppressWarnings({iq1<-min(which(x$surv<=0.75))
                    iq2<-min(which(x$surv<=0.5))
                    iq3<-min(which(x$surv<=0.25))})
  if (is.finite(iq1)) q1<-x$time[iq1] else q1<-Inf
  if (is.finite(iq2)) q2<-x$time[iq2] else q2<-Inf
  if (is.finite(iq3)) q3<-x$time[iq3] else q3<-Inf
 cat("Q1 =",round(q1,digits)," median =",round(q2,digits)," Q3 =",round(q3,digits),"\n")
 invisible(x)
}

print.svykmlist<-function(x, digits=3,...){
  cat("Weighted survival curves:\n")
  print(attr(x,"call"))
  for(i in 1:length(x)){
    cat(names(x)[i],": ")
    print(x[[i]],digits=digits,header=FALSE)
  }
  invisible(x)
}

quantile.svykm<-function(x, probs=c(0.75,0.5,0.25),ci=FALSE,level=0.95,...){
  
  iq<-sapply(probs, function(p) suppressWarnings(min(which(x$surv<=p))))
  qq<-sapply(iq, function(i) if (is.finite(i)) x$time[i] else Inf)
  names(qq)<-probs
  if (ci){
    if(is.null(x$varlog)){
      warning("no confidence interval available.")
    } else {
      halfalpha<-(1-level)/2
      z<-qnorm(halfalpha, lower.tail=FALSE)
      su<-exp(log(x$surv)+z*sqrt(x$varlog))
      iu<-sapply(probs, function(p) suppressWarnings(min(which(su<=p))))
      qu<-sapply(iu, function(i) if (is.finite(i)) x$time[i] else Inf)
      sl<-exp(log(x$surv)-z*sqrt(x$varlog))
      il<-sapply(probs, function(p) suppressWarnings(min(which(sl<=p))))
      ql<-sapply(il, function(i) if (is.finite(i)) x$time[i] else Inf)
      ci<-cbind(ql,qu)
      rownames(ci)<-probs
      colnames(ci)<-format(c(halfalpha,1-halfalpha),3)
      attr(qq,"ci")<-ci
    }
  }
  qq
}


confint.svykm<-function(object, parm, level=0.95,...){
  if (is.null(object$varlog)) stop("no standard errors in object")

  parm<-as.numeric(parm)
  idx<-sapply(parm, function(t) max(which(object$time<=t)))
  z<-qnorm((1-level)/2)
  ci<-exp(log(object$surv[idx])+outer(sqrt(object$varlog[idx]),c(z,-z)))
  ci[,2]<-pmin(ci[,2],1)
  rownames(ci)<-parm
  colnames(ci)<-format( c((1-level)/2, 1-(1-level)/2),3)
  ci
}


predict.svycoxph<-function(object, newdata, se=FALSE,
                           type=c("lp", "risk", "expected", "terms","curve"),
                           ...){
  
  type<-match.arg(type)
  if(type!="curve") return(NextMethod())
  
  design<-object$survey.design
  response<-object$y

  if (!is.null(attr(terms(object), "specials")$strata))
    stop("Stratified models are not supported yet")

  if (attr(response,"type")=="counting"){
    time<-object$y[,2]
    status<-object$y[,'status']
    entry<-object$y[,1]
  } else if (attr(response,'type')=="right"){
    time<-object$y[,"time"]
    status<-object$y[,"status"]
    entry<-rep(-Inf,length(time))
  } else stop("unsupported survival type")
  if(is.null(object$na.action)){
    design<-object$survey.design
  } else {
    design<-object$survey.design[-object$na.action,]
  }
  
  ff<-delete.response(terms(formula(object)))
  zmf<-model.frame(ff, newdata)
  z.pred<-model.matrix(ff, zmf)[,-1,drop=FALSE]
  
##
## The simple case first
##  
  risk<-getS3method("predict","coxph")(object,type="risk",se.fit=FALSE)
  if(se==FALSE){
    tt<-c(time,entry)
    ss<-c(status,rep(0,length(entry)))
    ee<-c(rep(1,length(status)),rep(-1,length(entry)))
    oo<-order(tt,-ee,-ss)
    dN<-ss[oo]
    w<-rep(weights(design),2)[oo]
    risks<-rep(risk,2)
    Y<-rev(cumsum(rev(risks[oo]*w*ee[oo])))
    keep<-dN>0
    s<-vector("list",nrow(z.pred))
    beta<-coef(object)
    h0<- cumsum( (w*dN/Y)[keep] )
    for(i in 1:nrow(z.pred)){
      zi<-z.pred[i,]-object$means
      s[[i]]<-list(time=time[oo][keep],
              surv=exp(-h0 * exp(sum(beta*zi))),
              call=sys.call())
    class(s[[i]])<-c("svykmcox","svykm")
    }
    names(s)<-rownames(newdata)
    return(s)
  }
##
## The hard case: curves with standard errors
##
  if(!inherits(design,"survey.design"))
    stop("replicate-weight designs not supported yet")
  
  keep<-which((status==1) & (weights(design)!=0))
  y<-outer(time,time[keep],">=")*risk*outer(entry,time[keep],"<=")
  dN<-diag(status)[,keep]
  oo<-order(time[keep], -status[keep])
  okeep<-keep[oo]
  ntimes<-length(oo)
  ttime<-time[okeep]
  sstatus<-status[okeep]
  totals<-svytotal(cbind(dN[,oo],y[,oo]), design)
  rm(dN)
  
  y<-coef(totals)[-(1:ntimes)]
  dNbar<-coef(totals)[1:ntimes]
  vtotals<-vcov(totals)
  rm(totals)
  
  h<-cumsum(dNbar/y)
  
  dVn<- vtotals[(1:ntimes),(1:ntimes)]/outer(y,y)
  dVy <- vtotals[-(1:ntimes),-(1:ntimes)]*outer(dNbar/y^2,dNbar/y^2)
  dCVny<- -vtotals[(1:ntimes),-(1:ntimes)]*outer(1/y,dNbar/y^2)
  dV<-dVn+dVy+dCVny+t(dCVny)
  
  det<-suppressWarnings(coxph.detail(object))
  ze<-sweep(as.matrix(det$means)[rep(1:length(det$time), det$nevent),,drop=FALSE],
            2, object$means)
  rm(det)
  
  dH<-dNbar/y
  h.ze<-dH*ze
  varbeta<-vcov(object)
  
  Vh<-numeric(ntimes)
  Vh[1]<-dV[1,1]
  for(i in 2:ntimes) Vh[i]<-Vh[i-1]+sum(dV[1:(i-1),i])+sum(dV[i,1:i])
  dVb<-numeric(ntimes)
  for(i in 1:ntimes) dVb[i]<-crossprod(h.ze[i,],varbeta%*%(h.ze[i,]))
  Vb<-cumsum(dVb)
  dCV<-matrix(nrow=ntimes,ncol=NCOL(ze))
  for(i in 1:ntimes) dCV[i,] <- -varbeta%*%(h.ze[i,])
  CV<-apply(dCV,2,cumsum)
  
  V0<-Vh+Vb
  s0<-exp(-h)
  
  s<-vector("list",nrow(z.pred))
  for(i in 1:nrow(z.pred)){
    zi<-z.pred[i,]-object$means
    riski<-exp(sum(zi*coef(object)))
    Vz<-crossprod(zi,varbeta%*%zi)*riski^2*h^2
    CVz<-colSums(t(dCV)*zi)*riski^2*h
    
    V<-V0*riski^2+Vz+CVz*2
    s[[i]]<-list(time=ttime,surv=exp(-h*riski), varlog=V)
    class(s[[i]])<-c("svykm.cox","svykm")
  }
  names(s)<-rownames(newdata)
  scall<-sys.call()
  scall[[1]]<-as.name(.Generic)
  attr(s,"call")<-scall
  class(s)<-c("svykmlist.cox","svykmlist")
  return(s)
}


