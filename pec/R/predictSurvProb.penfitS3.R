penalizedS3 <- function(formula,data,...){
  # {{{ distangle the formula
  ff <- as.character(formula)
  response <- formula(paste(ff[[2]],"~1",sep=""))
  terms <- strsplit(ff[[3]],"\\+|\\-")[[1]]
  terms <- sapply(terms,function(tt){## remove whitespace
    gsub(" ","",tt)
  })
  strippedTerms <- strsplit(terms,"[()]")
  # }}}
  # {{{ extract the penalized and unpenalized parts
  penalTerms <- sapply(strippedTerms,function(x){length(x)==2 && x[[1]]=="pen"})
  unpenalVarnames <- unlist(strippedTerms[penalTerms==FALSE])
  if (length(unpenalVarnames)>0){
    unpenalized <- formula(paste("~",paste(unpenalVarnames,collapse="+")))
    response <- update.formula(response,unpenalized)
  }
  penalizedVarnames <- unlist(sapply(strippedTerms[penalTerms==TRUE],
                                     function(x){strsplit(x[[2]],",")}),use.names=FALSE)
  penalizedVarPositions <- unlist(lapply(penalizedVarnames,function(x){
    if (length(splitter <- strsplit(x,":")[[1]])>1)
      seq(as.numeric(splitter)[1],as.numeric(splitter)[2],1)
    else
      match(x,names(data),nomatch=0)
  }),use.names=FALSE)
  penalizedVarPositions <- unique(penalizedVarPositions)
  ## print(penalizedVarPositions)
  if (any(tested <- (penalizedVarPositions>NCOL(data))|penalizedVarPositions<0))
    stop("Cannot find variable(s): ",names(data[tested]))
  penalized <- data[,penalizedVarPositions]
  # }}}
  # {{{ call S4 method
  ## unpenalized terms are communicated via
  ## the left hand side of response
  fitS4 <- penalized(response=response,
                     penalized=penalized,
                     data=data,
                     ...)
  # }}}
  # {{{ deliver S3 object
  fit <- list(fitS4=fitS4,call=match.call())
  class(fit) <- "penfitS3"
  fit
  # }}}
}


penalizedOpt <- function(formula,data=data,...){
  ## require(prodlim)
  argList <- prodlim::SmartControl(call=list(...), keys=c("profL1","optL1","penalized"))
  # {{{ distangle the formula
  ff <- as.character(formula)
  response <- formula(paste(ff[[2]],"~1",sep=""))
  terms <- strsplit(ff[[3]],"\\+|\\-")[[1]]
  terms <- sapply(terms,function(tt){## remove whitespace
    gsub(" ","",tt)
  })
  strippedTerms <- strsplit(terms,"[()]")
  # }}}
  # {{{ extract the penalized and unpenalized parts
  penalTerms <- sapply(strippedTerms,function(x){length(x)==2 && x[[1]]=="pen"})
  unpenalVarnames <- unlist(strippedTerms[penalTerms==FALSE])
  if (length(unpenalVarnames)>0){
    unpenalized <- formula(paste("~",paste(unpenalVarnames,collapse="+")))
    response <- update.formula(response,unpenalized)
  }
  penalizedVarnames <- unlist(sapply(strippedTerms[penalTerms==TRUE],
                                     function(x){strsplit(x[[2]],",")}),use.names=FALSE)
  penalizedVarPositions <- unlist(lapply(penalizedVarnames,function(x){
    if (length(splitter <- strsplit(x,":")[[1]])>1)
      seq(as.numeric(splitter)[1],as.numeric(splitter)[2],1)
    else
      match(x,names(data),nomatch=0)
  }),use.names=FALSE)
  penalizedVarPositions <- unique(penalizedVarPositions)
  ## print(penalizedVarPositions)
  if (any(tested <- (penalizedVarPositions>NCOL(data))|penalizedVarPositions<0))
    stop("Cannot find variable(s): ",names(data[tested]))
  penalized <- data[,penalizedVarPositions]
  # }}}
  # {{{ determine optimal lambdas
  ## prof <- do.call("profL1",
  ## c(list(response=response,penalized=penalized,data=data),
  ## argList$profL1))
  pen <- do.call("optL1",
                 c(list(response=response,penalized=penalized,data=data),
                   argList$optL1))
  # }}}
  # {{{ call S4 method
  fitOpt <- do.call("penalized",c(list(response=response,penalized=penalized,data=data),
                                  list(lambda1=pen$lambda),
                                  argList$penalized))
  # }}}
  # {{{ deliver S3 object
  fit <- list(fitS4=fitOpt,call=match.call())
  class(fit) <- "penfitS3"
  fit
  # }}}
}

##' @export 
predictSurvProb.penfitS3 <- function(object,
                                     newdata,
                                     times,
                                     ...){
    penfit <- object$fit
    pCovaNames <- names(penfit@penalized)
    newPen <- newdata[,pCovaNames]
    ptemp <- predict(penfit,penalized=newPen,data=newdata)
    ## require(prodlim)
    pos <- prodlim::sindex(jump.times=ptemp@time,eval.times=times)
    ## Remark: currently it is possible, but theoretically
    ## not allowed to carry predictions forward beyond the
    ## last jump.time
    p <- cbind(1,ptemp@curves)[,c(pos+1)]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

