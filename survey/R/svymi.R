
svydesign.imputationList<-function(ids, probs = NULL, strata = NULL, 
             variables = NULL, fpc = NULL, data, nest = FALSE, 
             check.strata = !nest,  weights = NULL, pps=FALSE,...){
    	designs <- lapply(data$imputations, function(d) svydesign(ids=ids, probs=probs,
              strata=strata,variables=variables,fpc=fpc,nest=nest,
              check.strata=check.strata, weights=weights,data=d,pps=pps,...))
    	rval <- list(designs=designs, call=sys.call(-1))
    	class(rval) <- "svyimputationList"
    	rval
    	}

svrepdesign.imputationList<-function(variables=NULL, repweights,weights,data,mse=getOption("survey.replicates.mse"),...){
  ## dispatch on data=
  if (!is.null(variables) && !inherits(variables,"imputationList"))
    stop("'variables' must also be an 'imputationList' (or NULL)")

  
  if(!is.null(variables)){
    if (inherits(repweights,"imputationList")){
      designs <- mapply(function(v,d,r) svrepdesign(variables=v, repweights=r, weights=weights,data=NULL,...),
                        variables$imputations,data$imputations, repweights$imputations,SIMPLIFY=FALSE)
    } else {
      designs <- mapply(function(d,v) svrepdesign(variables=v, repweights=repweights, weights=weights,data=d,...),
                      data$imputations,variables$imputations,SIMPLIFY=FALSE)
    }
  }else{
    if (inherits(repweights,"imputationList")){
      designs <- mapply(function(d,r) svrepdesign(repweights=r, weights=weights,data=NULL,...),
                        data$imputations, repweights$imputations,SIMPLIFY=FALSE)
    } else {
      designs <- lapply(data$imputations, function(d) svrepdesign( repweights=repweights, weights=weights,data=d,...))
    }
  }
  rval <- list(designs=designs, call=sys.call(-1))
  class(rval) <- "svyimputationList"
  rval
}

svydesign.DBimputationList<-function(ids, probs = NULL, strata = NULL, 
             variables = NULL, fpc = NULL, data, nest = FALSE, 
             check.strata = !nest,  weights = NULL, ...){
 
  design.vars<-c(all.vars(ids), all.vars(probs), all.vars(strata),all.vars(fpc), all.vars(weights))
  design.query<-paste("select", paste(design.vars,collapse=","), "from", data$imputations[1])
  if (data$db$dbtype=="ODBC")
    design.data<-sqlQuery(data$db$connection, design.query)
  else
    design.data<-dbGetQuery(data$db$connection, design.query)

  rval<-list()
  rval$design<-svydesign(ids=ids, probs=probs, strata=strata, data=design.data,
                  fpc=fpc, variables=variables, nest=nest,check.strata=check.strata,
                  weights=weights)
  class(rval$design)<-c(if(data$db$dbtype=="ODBC") "ODBCsvydesign" else "DBIsvydesign", class(rval$design))
  
  rval$design$updates<-data$updates
  rval$db<-data$db
  rval$imputations<-data$imputations
  rval$variables<-NULL
  rval$call<-sys.call(-1)
    class(rval)<-"svyDBimputationList"
  rval
}

print.svyDBimputationList<-function(x,...){
  cat("DB-backed Multiple (",length(x$imputations),") imputations: ",sep="")
  print(x$call)
}

print.svyimputationList<-function(x,...){
  cat("Multiple (",length(x$designs),") imputations: ",sep="")
  print(x$call)
}

dim.svyimputationList<-function(x){
  c(dim(x$designs[[1]]),length(x$designs))
}

dimnames.svyimputationList<-function(x){
   c(dimnames(x$designs[[1]]),list(paste("imputation",1:length(x$designs))))
}

subset.svyimputationList<-function(x, subset,...){
    n<-nrow(x$designs[[1]])
    e<-substitute(subset)
    r<-eval(e,x$designs[[1]]$variables, parent.frame())
    x$designs[[1]]<-x$designs[[1]][r,]
    same<-TRUE
    for(i in 2:length(x$designs)){
      r1<-eval(e,x$designs[[i]]$variables, parent.frame())
      x$designs[[i]]<-x$designs[[i]][r1,]
      r1<-r1 & !is.na(r1)
      if (any(r!=r1)) {
        same<-FALSE
      }
    }
    if (!same) warning('subset differed between imputations')
    
    x$call<-sys.call(-1)
    x
  }

subset.svyDBimputationList<-function(x, subset,...,all=FALSE){
    n<-nrow(x$designs[[1]])
    e<-substitute(subset)
    df<-getvars(all.vars(e), x$db$connection, x$imputations[1],
                db.only=FALSE, updates=x$design$updates)
    r<-eval(e,df, parent.frame())
    same<-TRUE
    for(i in 2:length(x$imputations)){
      df<-getvars(all.vars(e), x$db$connection, x$imputations[i],
                  db.only=FALSE, updates=x$design$updates)
      
      r1<-eval(e,df, parent.frame())
      r1<-r1 & !is.na(r1)
      if (any(r!=r1)) {
        same<-FALSE
        if (all) r <- r & r1 else r<- r | r1
      }
    }
    if (!same) warning('subset differed between imputations')
    x$design<-x$design[r,]
    x$call<-sys.call(-1)
    x
  }

with.svyimputationList<-function (data, expr, fun, ..., multicore=getOption("survey.multicore")) {
    pf <- parent.frame()
    if (multicore && !require("parallel",quietly=TRUE))
      multicore<-FALSE

    if (!is.null(match.call()$expr)) {
      expr <- substitute(expr)
      expr$design<-as.name(".design")
      if (multicore){
        results <- mclapply(data$designs,
                            function(.design) {
                            eval(expr, list(.design=.design),enclos=pf)
                          }
                            )
      } else{
        results <- lapply(data$designs,
                          function(.design) {
                            eval(expr, list(.design=.design),enclos=pf)
                          }
                          )
        
      }
    }
    else {
      results <- lapply(data$designs, fun, ...)
    }
    if (all(sapply(results, inherits, what = "imputationResult"))) {
      class(results) <- "imputationResultList"
      results$call <- sys.call(-1)
    }
    else {
      attr(results, "call") <- sys.call(-1)
    }
    results
  }


with.svyDBimputationList<-function (data, expr,  ..., multicore=getOption("survey.multicore")) {
    pf <- parent.frame()
    if (!is.null(match.call()$expr)) {
      expr <- substitute(expr)
      expr$design<-as.name(".design")
      if (multicore && !require("parallel")) multicore <-FALSE
      if (multicore){
        results<-mclapply(data$imputations,
                          function(tablename) {
                            close(data)
                            .design<-data$design
                            db<-data$db
                            db$tablename<-tablename
                            .design$db<-db
                            .design<-open(.design)
                            rval<-eval(expr, list(.design=.design),enclos=pf)
                            close(.design)
                            rval
                          }
                          )
      } else {
      results <- lapply(data$imputations,
                        function(tablename) {
                          .design<-data$design
                          db<-data$db
                          db$tablename<-tablename
                          .design$db<-db
                          eval(expr, list(.design=.design),enclos=pf)
                        }
                        )
    }
    }
    attr(results, "call") <- sys.call(-1)
    results
  }


update.svyDBimputationList<-function(object, ...){
  dots <- substitute(list(...))[-1]
  newnames <- names(dots)

  updates<-lapply(dots, function(dot){
    list(inputs=all.vars(dot),expression=dot)
  })

  if (is.null(object$design$updates))
    object$design$updates<-list(updates)
  else
    object$design$updates<-c(object$design$updates, list(updates))
  object
}

update.svyimputationList<-function(object, ...){
  dots <- substitute(list(...))[-1]
  newnames <- names(dots)
  for (i in seq(along = object$designs)) {
    for (j in seq(along = dots)) {
      object$designs[[i]]$variables[, newnames[j]] <- eval(dots[[j]], 
                           object$designs[[i]]$variables, parent.frame())
    }
  }
  object
}

close.svyDBimputationList<-function(con,...){
  dbcon<-con$db$connection
  if (is(dbcon,"DBIConnection"))
    dbDisconnect(dbcon)
  else
    close(dbcon)
  invisible(con)
}

open.svyDBimputationList<-function(con,...){
  if(con$db$dbtype=="ODBC"){
    oldenc<-attr(con$db$connection)
    con$db$connection<-odbcReConnect(con$db$connection,...)
    attr(con$db$connection,"encoding")<-oldenc
  } else {
    dbdriver<-dbDriver(con$db$dbtype)
    con$db$connection<-dbConnect(dbdriver,dbname=con$db$dbname,...)
  }
  con
}
