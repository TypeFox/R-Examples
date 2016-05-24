
imputationList<-function(datasets,...) UseMethod("imputationList")

imputationList.character<-function(datasets, dbtype, dbname, ...){
  if(dbtype=="ODBC"){
    connection<-RODBC::odbcConnect(dbname,...)
  } else {
    driver<-DBI::dbDriver(dbtype)
    connection<-DBI::dbConnect(driver,dbname,...)
  }
  
  rval<-list(imputations=datasets, db=list(connection=connection, dbname=dbname,dbtype=dbtype,...), call=sys.call(-1))
  class(rval)<-c("DBimputationList","imputationList")
  rval
}

imputationList.default<-function(datasets,...){
  rval<-list(imputations=datasets, call=sys.call(-1))
  class(rval)<-"imputationList"
  rval
}

print.imputationList<-function(x,...){
  cat("MI data with", length(x$imputations),"datasets\nCall: ")
  print(x$call)
}

with.DBimputationList<-function(data,expr,...){
  ee<-substitute(expr)
  imputations<-lapply(data$imputations,
                           function(thistable) getvars(all.vars(ee), data$db$connection, thistable,
                                                       db.only=FALSE,updates=data$updates))
  pf<-parent.frame()

  expr<-substitute(expr)
  results<-lapply(imputations, function(dataset) eval(expr, dataset, enclos=pf))
  

  if (all(sapply(results, inherits,  what="imputationResult"))){
    class(results)<-"imputationResultList"
    results$call<-sys.call(-1)
  } else {
    attr(results,"call")<-sys.call(-1)
  }
  
  results

}

close.DBimputationList<-function(con,...){
  dbcon<-con$db$connection
  if (is(dbcon,"DBIConnection"))
    DBI::dbDisconnect(dbcon)
  else
    close(dbcon)
  invisible(con)
}

open.DBimputationList<-function(con,...){
  if(con$db$dbtype=="ODBC"){
    oldenc<-attr(con$db$connection)
    con$db$connection<-RODBC::odbcReConnect(con$db$connection,...)
    attr(con$db$connection,"encoding")<-oldenc
  } else {
    dbdriver<-DBI::dbDriver(con$db$dbtype)
    con$db$connection<-DBI::dbConnect(dbdriver,dbname=con$db$dbname,...)
  }
  con
}

with.imputationList<-function(data, expr,fun,...){

  pf<-parent.frame()
  if (!is.null(match.call()$expr)){
    expr<-substitute(expr)
    results<-lapply(data$imputations, function(dataset) eval(expr, dataset, enclos=pf))
  } else {
    results<-lapply(data$imputations, fun,...)
  }

  if (all(sapply(results, inherits,  what="imputationResult"))){
    class(results)<-"imputationResultList"
    results$call<-sys.call(-1)
  } else {
    attr(results,"call")<-sys.call(-1)
  }
  
  results

  
}


MIcombine<-function(results, ...) UseMethod("MIcombine")


MIcombine.imputationResultList<-function(results, call=NULL,df.complete=Inf,...){

  vars<-suppressWarnings(lapply(results, vcov))
  thetas<-lapply(results, coef)
  rval<-MIcombine(thetas, vars,call=sys.call(-1),df.complete=df.complete)
  rval$call<-c(results$call, call)
  rval
}

MIcombine.default<-function(results, variances,call=sys.call(),df.complete=Inf,...){

  m<-length(results)
  oldcall<-attr(results,"call")
  if (missing(variances)){
    variances<-suppressWarnings(lapply(results, vcov))
    results<-lapply(results,coef)
  }
  ## yes, we do need a loop here.
  vbar<-variances[[1]]
  cbar<-results[[1]]
  for(i in 2:m){
    cbar<-cbar+results[[i]]
    vbar<-vbar+variances[[i]]
  }
  cbar<-cbar/m
  vbar<-vbar/m

  evar<-var(do.call("rbind",results))  
  r<- (1+1/m)*evar/vbar
  df<-(m-1)*(1+ 1/r)^2
  if (is.matrix(df))
    df<-diag(df)
  if(is.finite(df.complete)){
    dfobs<-((df.complete+1)/(df.complete+3))*df.complete*vbar/(vbar+evar)
    if(is.matrix(dfobs))
      dfobs<-diag(dfobs)
    df<-1/(1/dfobs+1/df)
  }
  if (is.matrix(r))
    r<-diag(r)
  rval<-list(coefficients=cbar,
             variance=vbar+evar*(m+1)/m,
             call=c(oldcall, call),
             nimp=m,
             df=df,
             missinfo=(r+2/(df+3))/(r+1)
           )
  class(rval)<-"MIresult"
  rval
}

vcov.MIresult<-function(object,...) object$variance

MIextract<-function(results, expr, fun){
  pf<-parent.frame()
  if (!is.null(match.call()$expr)){
    expr<-substitute(expr)
    lapply(results, function(result) eval(expr, result,pf))
  } else {
    lapply(results, fun)
  }
  
}

print.MIresult<-function(x,...){
  cat("Multiple imputation results:\n")
  lapply(x$call, function(a) {cat("      ");print(a)})
  out<-data.frame(results=coef(x), se=sqrt(diag(vcov(x))))
  print(out)
}

summary.MIresult<-function(object,...,alpha=0.05, logeffect=FALSE){
  cat("Multiple imputation results:\n")
  lapply(object$call, function(a) {cat("      ");print(a)})
  out<-data.frame(results=coef(object), se=sqrt(diag(vcov(object))))
  crit<-qt(alpha/2,object$df, lower.tail=FALSE)
  out$"(lower"<-out$results-crit*out$se
  out$"upper)"<-out$results+crit*out$se
  if (logeffect){
    out$results<-exp(out$results)
    out$se<-out$se*out$results
    out$"(lower"<-exp(out$"(lower")
    out$"upper)"<-exp(out$"upper)")
  }
  out$"missInfo" <-paste(round(100*object$missinfo), "%")
  print(out,...)

}


update.imputationList<-function(object,...){

  dots<-substitute(list(...))[-1]
  newnames<-names(dots)
  
  for(i in seq(along=object$imputations)){
    for(j in seq(along=dots)){
      object$imputations[[i]][,newnames[j]]<-eval(dots[[j]],object$imputations[[i]], parent.frame())
    }
  }
  
  object
}


dimnames.imputationList<-function(x) c(dimnames(x$imputations[[1]]), paste("imputation",seq(along=x$imputations)))
dim.imputationList<-function(x) c(dim(x$imputations[[1]]),length(x$imputations))


rbind.imputationList<-function(...){
  dots<-list(...)
  if(length(dots)==1)
    return(dots)
  if(!all(sapply(dots, inherits,"imputationList")))
    stop("arguments must all be imputationLists")

  ms<-sapply(dots, function(dot) length(dot$imputations))
  if (length(unique(ms))>1)
    stop("differing numbers of imputations")

  rval<-list(imputations=lapply(1:ms[1], function(i)
               do.call("rbind", lapply(dots, function(dot) dot$imputations[[i]]))),
             call=sys.call())
  class(rval)<-"imputationList"
  rval
}


cbind.imputationList<-function(...){
  dots<-list(...)
  if(length(dots)==1)
    return(dots)
  if(!all(sapply(dots, inherits,"imputationList")))
    stop("arguments must all be imputationLists")

  ms<-sapply(dots, function(dot) length(dot$imputations))
  if (length(unique(ms))>1)
    stop("differing numbers of imputations")

  rval<-list(imputations=lapply(1:ms[1], function(i)
               do.call("cbind", lapply(dots, function(dot) dot$imputations[[i]]))),
             call=sys.call())
  class(rval)<-"imputationList"
  rval
}
