
##
## stored variable updates.
##

updatesInfilter<-function(varlist, updates){
  if (is.null(updates)) return(list(varlist=varlist))
  n<-length(updates)
  v<-vector("list",n)
  for(i in n:1){
    if (any(idx<-(varlist %in% names(updates[[i]])))){
      v[[i]]<-varlist[idx]
      ups<-match(v[[i]], names(updates[[i]]))
      varlist<-unique(c(varlist[!idx], do.call(c, lapply(updates[[i]][ups], "[[", "inputs"))))
    }
  }
 list(varlist=varlist, history=v)
}

updatesOutfilter<-function(df, varlist,history, updates){
  if (is.null(updates)) return(df)
  if (all(sapply(history,length)==0)) return(df)
  n<-length(updates)
  for(i in 1:n){
    if (mi<-length(history[[i]])){
      outputs<-vector("list", mi)
      for(j in 1:mi){
        idx.j<-match(history[[i]][j],names(updates[[i]]))
        outputs[[j]]<-eval(updates[[i]][[idx.j]]$expression, df)
      }
      names(outputs)<-history[[i]]
      if (any(mod<-names(df) %in% names(outputs))){
        df<-df[,!mod,drop=FALSE]
      }
      df<-cbind(df,outputs)
    }
  }
  df[, names(df) %in% varlist,drop=FALSE]
}



getvars<-function (vars, dbconnection, tables, db.only = TRUE, updates=NULL) 
{
  if(!length(vars)) return(NULL)
  var0<-vars
  infilter<-updatesInfilter(var0, updates)
  if (db.only) {
    in.db <- infilter$varlist
  }
  else {
    query <- sub("@tab@", tables, "select * from @tab@ limit 1")
    oneline <- DBI::dbGetQuery(dbconnection, query)
    in.db <- infilter$varlist[infilter$varlist %in% names(oneline)]
  }
  query <- paste("select", paste(in.db, collapse = ", "), "from", 
                 tables)
  
  if (is(dbconnection, "DBIConnection"))
    df <- DBI::dbGetQuery(dbconnection, query)
  else ##ODBC
    df<-RODBC::sqlQuery(dbconnection, query)
  
  df<-updatesOutfilter(df, var0, infilter$history, updates)
  
  is.string <- sapply(df, is.character)
  if (any(is.string)) {
    for (i in which(is.string)) df[[i]] <- as.factor(df[[i]])
  }
  df
}


update.DBimputationList<-function(object, ...){
  dots <- substitute(list(...))[-1]
  newnames <- names(dots)

  updates<-lapply(dots, function(dot){
    list(inputs=all.vars(dot),expression=dot)
  })

  if (is.null(object$updates))
    object$updates<-list(updates)
  else
    object$updates<-c(object$updates, list(updates))
  object
}



rbind.DBimputationList<-function(...) stop("Not implemented: use database facilities")
cbind.DBimputationList<-function(...) stop("Not implemented: use database facilities")
