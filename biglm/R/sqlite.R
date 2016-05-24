
setMethod("bigglm",
          c("ANY","DBIConnection"),
          function(formula, data, family = gaussian(),
                                  tablename, ..., chunksize=5000){
            terms<-terms(formula)
            modelvars<-all.vars(formula)	
            dots<-as.list(substitute(list(...)))[-1]
            dotvars<-unlist(lapply(dots,all.vars))
            vars<-unique(c(modelvars,dotvars))
            query<-paste("select ",paste(vars,collapse=", ")," from ",tablename)
            result<-dbSendQuery(data, query)
            got<-0
            on.exit(dbClearResult(result))
            chunk<-function(reset=FALSE){
              if(reset){
                if(got>0){
                  dbClearResult(result)
                  result<<-dbSendQuery(data,query)
                  got<<-0
                }
                return(TRUE)
              }
              rval<-fetch(result,n=chunksize)	
              got<<-got+nrow(rval)
              if (nrow(rval)==0) 
                return(NULL)
              return(rval)
            }
            rval<-bigglm(formula, data=chunk, family=family, ...)
            rval$call<-sys.call()
            rval$call[[1]]<-as.name("bigglm")
            rval
          }
          )
		
