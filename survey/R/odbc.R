

print.ODBCsvydesign<-function(x,...){
  cat("ODBC-backed ")
  NextMethod()
   if (!checkConnection(x$db$connection, error=FALSE))
    cat("<ODBC Connection closed>\n")
  invisible(x)
}

summary.ODBCsvydesign<-function(object,...){
   class(object)<-c("summary.ODBCsvydesign",class(object))
   object
}

print.summary.ODBCsvydesign<-function(x,...){
   print.survey.design2(x,varnames=TRUE,design.summaries=TRUE,...)
   invisible(x)
}

close.ODBCsvydesign<-function(con,...){
  close(con$db$connection,...)
  invisible(con)
}

open.ODBCsvydesign<-function(con,...){
  oldenc<-attr(con$db$connection,"encoding")  ## bug in RODBC 1.2-3
  con$db$connection<-odbcReConnect(con$db$connection,...)
  attr(con$db$connection,"encoding")<-oldenc
  con
}

svymean.ODBCsvydesign<-function(x, design,...){
  design$variables<-getvars(x, design$db$connection, design$db$tablename, updates=design$updates)
  NextMethod("svymean",design)
}


svytotal.ODBCsvydesign<-function(x, design,na.rm=FALSE,...){
  design$variables<-getvars(x, design$db$connection, design$db$tablename, updates=design$updates)
  NextMethod("svytotal",design)
}

svyquantile.ODBCsvydesign<-function(x, design,quantiles,...){
  design$variables<-getvars(x, design$db$connection, design$db$tablename, updates=design$updates)
  NextMethod("svyquantile",design)
}

svyglm.ODBCsvydesign<-function(formula, design,...){
  design$variables<-dropFactor(getvars(formula, design$db$connection, design$db$tablename,updates=design$updates),
                               weights(design))
  NextMethod("svyglm",design)
}


svyranktest.ODBCsvydesign<-function(formula, design,...){
  design$variables<-dropFactor(getvars(formula, design$db$connection, design$db$tablename,updates=design$updates),
                               weights(design))
  NextMethod("svyranktest",design)
}


svyplot.ODBCsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename, updates=design$updates)
  design$variables[weights(design)==0,]<-NA
  NextMethod("svyplot",design)
}


svycoplot.ODBCsvydesign<-function(formula,design, style=c("hexbin","transparent"),
                            basecol="black",alpha=c(0,0.8),hexscale=c("relative","absolute"),...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  design$variables[weights(design)==0,]<-NA
  NextMethod("svycoplot",design)
}

svyboxplot.ODBCsvydesign<-function(formula,design, ...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  design$variables[weights(design)==0,]<-NA
  class(design)<-setdiff(class(design),"ODBCsvydesign")
  svyboxplot(formula,design,...)
}

svycdf.ODBCsvydesign<-function(formula,design, na.rm=TRUE, ...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svycdf",design)

}

svyolr.ODBCsvydesign<-function(formula,design,...){
  design$variables<-dropFactor(getvars(formula, design$db$connection, design$db$tablename,updates=design$updates),
                               weights(design))
  NextMethod("svyolr",design)
}

svycoxph.ODBCsvydesign<-function(formula,design,...){
  design$variables<-dropFactor(getvars(formula, design$db$connection, design$db$tablename,updates=design$updates),
                               weights(design))
  NextMethod("svycoxph",design)
}

svyvar.ODBCsvydesign<-function(x,design,na.rm=FALSE,...){
  design$variables<-getvars(x, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svyvar",design)
}



svykm.ODBCsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svykm",design)
}


svykappa.ODBCsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svykappa",design)
}


svysmooth.ODBCsvydesign<-function(formula,design,method=c("locpoly","quantreg"),bandwidth,quantile,df,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svysmooth",design)
}


svychisq.ODBCsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svychisq",design)
}

svyratio.ODBCsvydesign<-function(numerator, denominator, design,...){
  design$variables<-cbind(getvars(numerator,design$db$connection, design$db$tablename,updates=design$updates),
                          getvars(denominator,design$db$connection, design$db$tablename,updates=design$updates))
  NextMethod("svyratio",design)

}


svyby.ODBCsvydesign<-function(formula, by, design,...){
  design$variables<-cbind(getvars(formula,design$db$connection, design$db$tablename,updates=design$updates),
                          getvars(by,design$db$connection, design$db$tablename,updates=design$updates))
  class(design)<-setdiff(class(design),"ODBCsvydesign")
  svyby(formula,by,design,...)
}

svytable.ODBCsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svytable",design)
}

calibrate.ODBCsvydesign<-function(design,formula,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("calibrate",design)
}
postStratify.ODBCsvydesign<-function(design, strata, population, partial = FALSE, ...) .NotYetImplemented()

subset.ODBCsvydesign<-function (x, subset, ...) 
{
    e <- substitute(subset)
    x$variables<-getvars(make.formula(all.vars(e)), x$db$connection, x$db$tablename,updates=x$updates)
    r <- eval(e, x$variables, parent.frame())
    r <- r & !is.na(r)
    x <- x[r, ]
    x$call <- sys.call(-1)
    x
}


dim.ODBCsvydesign<-function(x){
  w<-weights(x)
  nrow<-sum(w!=0)
  coln<-names(sqlQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
  if (!is.null(x$updates)){
    update.names<-do.call(c, lapply(x$updates, names))
    ncol<-length(unique(c(coln,update.names)))
  } else ncol<-length(coln)
  c(nrow,ncol)
}


dimnames.ODBCsvydesign<-function(x){
  w<-weights(x)
  rown<-rownames(x$cluster)[w!=0]
   coln<-names(sqlQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
   if (!is.null(x$updates)){
     update.names<-do.call(c, lapply(x$updates, names))
     coln<-unique(c(coln,update.names))
   }
   list(rown,coln)
}


"[.ODBCsvydesign"<-function (x, i, ..., drop = TRUE) 
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

