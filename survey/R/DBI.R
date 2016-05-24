
svydesign.character<-function (ids, probs = NULL, strata = NULL, variables = NULL, 
                               fpc = NULL, data, nest = FALSE, check.strata = !nest,
                               weights = NULL,pps=FALSE,
                               dbtype="SQLite", dbname,
                               ...) 
{

  if (dbtype == "ODBC"){
    library(RODBC)
    if (dbname=="")
      dbconn<-odbcDriverConnect(dbname,...)
    else
      dbconn<-odbcConnect(dbname,...)
  } else {
    db<-dbDriver(dbtype)
    dbconn<- dbConnect(db, dbname,...)
  }
  design.vars<-c(all.vars(ids), all.vars(probs), all.vars(strata),
                 all.vars(fpc), all.vars(weights))
  design.query<-paste("select", paste(design.vars,collapse=","), "from", data)
  if (dbtype=="ODBC")
    design.data<-sqlQuery(dbconn, design.query)
  else
    design.data<-dbGetQuery(dbconn, design.query)
    
  rval<-svydesign(ids=ids, probs=probs, strata=strata, data=design.data,
                  fpc=fpc, variables=variables, nest=nest,check.strata=check.strata,
                  weights=weights)
  rval$db<-list(dbname=dbname, tablename=data, connection=dbconn, dbtype=dbtype)
  rval$variables<-NULL
  rval$call<-sys.call(-1)
  if (dbtype=="ODBC")
    class(rval)<-c("ODBCsvydesign",class(rval))
  else
    class(rval)<-c("DBIsvydesign",class(rval))
  rval
}

print.DBIsvydesign<-function(x,...){
  cat("DB-backed ")
  NextMethod()
   if (!checkConnection(x$db$connection, error=FALSE))
    cat("<DBI Connection closed>\n")
  invisible(x)
}

summary.DBIsvydesign<-function(object,...){
   class(object)<-c("summary.DBIsvydesign",class(object))
   object
}

print.summary.DBIsvydesign<-function(x,...){
   print.survey.design2(x,varnames=TRUE,design.summaries=TRUE,...)
   invisible(x)
}

close.DBIsvydesign<-function(con,...){
  dbDisconnect(con$db$connection,...)
  invisible(con)
}

open.DBIsvydesign<-function(con,...){
  db<-dbDriver(con$db$dbtype)
  con$db$connection<-dbConnect(db, dbname=con$db$dbname,...)
  con
}

svymean.DBIsvydesign<-function(x, design,...){
  design$variables<-getvars(x, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("svymean",design)
}


svytotal.DBIsvydesign<-function(x, design,na.rm=FALSE,...){
  design$variables<-getvars(x, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("svytotal",design)
}

svyquantile.DBIsvydesign<-function(x, design,quantiles,...){
  design$variables<-getvars(x, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("svyquantile",design)
}


dropFactor<-function(mf, w){
  if(!any(w==0)) return(mf)
  dropped<-w==0
  for(i in 1:ncol(mf)) {
    if (is.factor(mf[[i]])){
      fi<-mf[[i]]
      if (all(dropped[fi==levels(fi)[1]])){
        tt<-table(fi[!dropped])
        l<-min(which(tt>0))
        levs<-levels(fi)
        mf[[i]]<-relevel(mf[[i]],ref=levs[l])
      }
    }
  }
  mf
}

svyglm.DBIsvydesign<-function(formula, design,...){
  design$variables<-dropFactor(getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset),
                               weights(design))
  NextMethod("svyglm",design)
}



svyplot.DBIsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  design$variables[weights(design)==0,]<-NA
  NextMethod("svyplot",design)
}


svycoplot.DBIsvydesign<-function(formula,design, style=c("hexbin","transparent"),
                            basecol="black",alpha=c(0,0.8),hexscale=c("relative","absolute"),...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename, updates=design$updates, subset=design$subset)
  design$variables[weights(design)==0,]<-NA
  NextMethod("svycoplot",design)
}

svyboxplot.DBIsvydesign<-function(formula,design, ...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  design$variables[weights(design)==0,]<-NA
  class(design)<-setdiff(class(design),"DBIsvydesign")
  svyboxplot(formula,design,...)
}


svycdf.DBIsvydesign<-function(formula,design, na.rm=TRUE, ...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("svycdf",design)

}

svyolr.DBIsvydesign<-function(formula,design,...){
  design$variables<-dropFactor(getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset),
                               weights(design))
  NextMethod("svyolr",design)
}

svycoxph.DBIsvydesign<-function(formula,design,...){
  design$variables<-dropFactor(getvars(formula, design$db$connection, design$db$tablename,updates=design$updates),
                               weights(design))
  NextMethod("svycoxph",design)
}

svyvar.DBIsvydesign<-function(x,design,na.rm=FALSE,...){
  design$variables<-getvars(x, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("svyvar",design)
}



svykm.DBIsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("svykm",design)
}


svykappa.DBIsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("svykappa",design)
}


svysmooth.DBIsvydesign<-function(formula,design,method=c("locpoly","quantreg"),bandwidth,quantile,df,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("svysmooth",design)
}


svychisq.DBIsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("svychisq",design)
}

svyranktest.DBIsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("svyranktest",design)
}

svyratio.DBIsvydesign<-function(numerator, denominator, design,...){
  design$variables<-cbind(getvars(numerator,design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset),
                          getvars(denominator,design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset))
  NextMethod("svyratio",design)

}


svyby.DBIsvydesign<-function(formula, by, design,...){
  design$variables<-cbind(getvars(formula,design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset),
                          getvars(by,design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset))
  class(design)<-setdiff(class(design),"DBIsvydesign")
  svyby(formula,by,design,...)
}

svytable.DBIsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("svytable",design)
}


calibrate.DBIsvydesign<-function(design,formula,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates, subset=design$subset)
  NextMethod("calibrate",design)
}
postStratify.DBIsvydesign<-function(design, strata, population, partial = FALSE, ...) .NotYetImplemented()




subset.DBIsvydesign<-function (x, subset, ...) 
{
    e <- substitute(subset)
    x$variables<-getvars(make.formula(all.vars(e)), x$db$connection, x$db$tablename,updates=x$updates, subset=x$subset)
    r <- eval(e, x$variables, parent.frame())
    r <- r & !is.na(r)
    x <- x[r, ]
    x$call <- sys.call(-1)
    x
}





dim.DBIsvydesign<-function(x){
  w<-weights(x)
  nrow<-sum(w!=0)
   coln<-names(dbGetQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
   if (!is.null(x$updates)){
     update.names<-do.call(c, lapply(x$updates, names))
     ncol<-length(unique(c(coln,update.names)))
   } else ncol<-length(coln)
   c(nrow,ncol)
}


dim.DBIrepdesign<-function(x){
  if (is.null(x$subset))
    nrow <-nrow(x$repweights)
  else
    nrow<-length(x$subset)
  coln<-names(dbGetQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
  if (!is.null(x$updates)){
    update.names<-do.call(c, lapply(x$updates, names))
    ncol<-length(unique(c(coln,update.names)))
  } else ncol<-length(coln)
  c(nrow,ncol)
}

dimnames.DBIsvydesign<-function(x){
   rown<-rownames(x$cluster)[weights(x)!=0]
   coln<-names(dbGetQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
   if (!is.null(x$updates)){
     update.names<-do.call(c, lapply(x$updates, names))
     coln<-unique(c(coln,update.names))
   }
   list(rown,coln)
}


dimnames.DBIrepdesign<-function(x){
   if (is.null(x$subset))
     rown<-rownames(x$cluster)
   else 
     rown<-rownames(x$cluster)[x$subset]
   coln<-names(dbGetQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
   if (!is.null(x$updates)){
     update.names<-do.call(c, lapply(x$updates, names))
     coln<-unique(c(coln,update.names))
   }
   list(rown,coln)
}

"[.DBIsvydesign"<-function (x, i, ..., drop = TRUE) 
{
  if (!missing(i)) {
    if (is.logical(i)) 
      x$prob[!i] <- Inf
    else if (is.numeric(i) && length(i)) 
      x$prob[-i] <- Inf
    else {
      tmp <- x$prob[i, ]
      x$prob <- rep(Inf, length(x$prob))
      x$prob[i, ] <- tmp
    }
    index <- is.finite(x$prob)
    psu <- !duplicated(x$cluster[index, 1])
    tt <- table(x$strata[index, 1][psu])
    if (any(tt == 1)) {
      warning(sum(tt == 1), " strata have only one PSU in this subset.")
    }
    
  }
  else {
    if (!is.null(x$variables)) 
      x$variables <- x$variables[, ..1, drop = FALSE]
  }
  x
}


"[.DBIrepdesign"<-function (x, i, j, drop = FALSE) 
{
    if (!missing(i)) {
        pwt <- x$pweights
        if (is.data.frame(pwt)) 
            pwt <- pwt[[1]]
        x$pweights <- pwt[i]
        x$repweights <- x$repweights[i, , drop = FALSE]
        if (!is.null(x$selfrep)) 
            x$selfrep <- x$selfrep[i]
        if (is.null(x$subset))
          x$subset<-(1:nrow(x$variables))[i]
        else
          x$subset<-x$subset[i]
        if (!missing(j)) 
            x$variables <- x$variables[i, j, drop = FALSE]
        else x$variables <- x$variables[i, , drop = FALSE]
        x$degf <- NULL
        x$degf <- degf(x)
    }
    else {
        x$variables <- x$variables[, j, drop = FALSE]
    }
    x
}


svrepdesign.character<-function (variables=NULL,repweights=NULL, weights=NULL,
                                 data=NULL,type=c("BRR","Fay","JK1", "JKn","bootstrap","other"),
                                 combined.weights=TRUE, rho=NULL, bootstrap.average=NULL,
                                 scale=NULL,rscales=NULL,fpc=NULL, fpctype=c("fraction","correction"),
                                 mse=getOption("survey.replicates.mse"),dbtype="SQLite", dbname,
                                 ...) 
{

  if (dbtype == "ODBC"){
    library(RODBC)
    if (dbname=="")
      dbconn<-odbcDriverConnect(dbname,...)
    else
      dbconn<-odbcConnect(dbname,...)
  } else {
    db<-dbDriver(dbtype)
    dbconn<- dbConnect(db, dbname,...)
  }
  if (is.character(repweights)){
    allvars<-names(dbGetQuery(dbconn, paste("select * from",data,"limit 1"))) 
    design.vars<-c(all.vars(weights),grep(repweights,allvars,value=TRUE))
  } else {
    design.vars<-c(all.vars(weights),all.vars(repweights))
  }

  
  design.query<-paste("select", paste(design.vars,collapse=","), "from", data)
  if (dbtype=="ODBC")
    design.data<-sqlQuery(dbconn, design.query)
  else
    design.data<-dbGetQuery(dbconn, design.query)
    
  rval<-svrepdesign(variables=variables,repweights=repweights, weights=weights, type=type,
                    data=design.data,
                    combined.weights=combined.weights, rho=rho, bootstrap.average=NULL,
                    scale=scale,rscales=rscales,fpc=fpc, fpctype=c("fraction","correction"))
  
  rval$db<-list(dbname=dbname, tablename=data, connection=dbconn, dbtype=dbtype)
  rval$variables<-NULL
  rval$call<-sys.call(-1)
  if (dbtype=="ODBC")
    class(rval)<-c("ODBCrepdesign","ODBCsvydesign",class(rval))
  else
    class(rval)<-c("DBIrepdesign","DBIsvydesign",class(rval))
  rval
}

print.DBIrepdesign<-function(x,...){
  cat("DB-backed replicate weight design\n")
  print.svyrep.design(x,...)
  if (!checkConnection(x$db$connection, error=FALSE))
    cat("<DBI Connection closed>\n")
  invisible(x)
}

print.ODBCrepdesign<-function(x,...){
  cat("ODBC-backed replicate weight design\n")
  print.svyrep.design(x,...)
   if (!checkConnection(x$db$connection, error=FALSE))
    cat("<ODBC Connection closed>\n")
  invisible(x)
}

summary.DBIrepdesign<-function(object,...){
   summary.svyrep.design(object,...)
}

summary.ODBCrepdesign<-function(object,...){
   summary.svyrep.design(object,...)
}


