## setGeneric("formula", function(x, env = parent.frame(), ...) {
##           standardGeneric("formula")})
## don't know why behaviour of anova() and formula() are different?
## (used setGeneric for anova() without trouble, caused problems here)
## trying to avoid "creating a new generic" message on install?

setMethod("formula", "mle2",
          function(x, env = parent.frame(), ...) {
            as.formula(x@formula)
          })

## stdEr <- function(x, ...) {
##   UseMethod("stdEr")
##  }

setGeneric("stdEr", function(x, ...) { standardGeneric("stdEr")})

setMethod("stdEr","mle2",
          function(x, ...) {
            sqrt(diag(x@vcov)) ## why doesn't vcov(x) work here???
          })

## should this be object@fullcoef or object@coef??? or should
## it have an additional argument --- is that possible?
setMethod("coef", "mle2", function(object,exclude.fixed=FALSE) {
    if (!exclude.fixed) object@fullcoef else object@coef
})
## fullcoef <- function(object) object@fullcoef  ## this should be a method
setMethod("coef", "summary.mle2", function(object) { object@coef })
## hmmm.  Work on this. 'hessian' conflicts with numDeriv definition. Override?
## setMethod("Hessian", sig="mle2", function(object) { object@details$hessian })

setMethod("show", "mle2", function(object){
    cat("\nCall:\n")
    print(object@call.orig)
    cat("\nCoefficients:\n")
    print(coef(object))
    cat("\nLog-likelihood: ")
    cat(round(as.numeric(logLik(object)),2),"\n")
    if (object@optimizer=="optimx" && length(object@method)>1) {
      cat("Best method:",object@details$method.used,"\n")
    }
    if (object@details$convergence>0)
      cat("\nWarning: optimization did not converge (code ",
          object@details$convergence,": ",object@details$message,")\n",sep="")
  })

setMethod("show", "summary.mle2", function(object){
    cat("Maximum likelihood estimation\n\nCall:\n")
    print(object@call)
    cat("\nCoefficients:\n")
    printCoefmat(coef(object))
    cat("\n-2 log L:", object@m2logL, "\n")
})

setMethod("show", "profile.mle2", function(object){
    cat("Likelihood profile:\n\n")
    print(object@profile)
  })

setMethod("summary", "mle2", function(object, waldtest=TRUE, ...){
    cmat <- cbind(Estimate = object@coef,
                  `Std. Error` = sqrt(diag(object@vcov)))
    zval <- cmat[,"Estimate"]/cmat[,"Std. Error"]
    pval <- 2*pnorm(-abs(zval))
    coefmat <- cbind(cmat,"z value"=zval,"Pr(z)"=pval)
    m2logL <- 2*object@min
    new("summary.mle2", call=object@call.orig, coef=coefmat, m2logL= m2logL)
})



setMethod("logLik", "mle2",
function (object, ...)
{
    if(length(list(...)))
        warning("extra arguments discarded")
    val <- -object@min
    attr(val, "df") <- length(object@coef)
    attr(val, "nobs") <- attr(object,"nobs")
    class(val) <- "logLik"
    val
  })

setGeneric("deviance", function(object, ...) standardGeneric("deviance"))
setMethod("deviance", "mle2",
function (object, ...)
{
  -2*logLik(object)
})

setMethod("vcov", "mle2", function (object, ...) { object@vcov } )


setGeneric("anova", function(object, ...) standardGeneric("anova"))
setMethod("anova","mle2",
          function(object,...,width=getOption("width"), exdent=10) {
            mlist <- c(list(object),list(...))
            ## get names from previous call
            mnames <- sapply(sys.call(sys.parent())[-1],deparse)
            ltab <- as.matrix(do.call("rbind",
                                      lapply(mlist,
                                             function(x) {
                                               c("Tot Df"=length(x@coef),
                                                 Deviance=-2*logLik(x))
                                             })))
            terms=sapply(mlist,
              function(obj) {
                if (is.null(obj@formula) || obj@formula=="") {
                  mfun <- obj@call$minuslogl
                  mfun <- paste("[",if (is.name(mfun)) {
                    as.character(mfun)
                  } else { "..." },
                                "]",sep="")
                  paste(mfun,": ",paste(names(obj@coef),
                                        collapse="+"),sep="")
                } else {
                  as.character(obj@formula)
                }
              })
            mterms <- paste("Model ",
                            1:length(mnames),": ",mnames,", ",terms,sep="")
            mterms <- strwrapx(mterms,width=width,exdent=exdent,
                               wordsplit="[ \n\t]")
  ## trunc.term <- function(s,len) {
  ##     ## cat("***",nchar(s),length(grep("\\+",s)),"\n",sep=" ")    
  ##     if ((nchar(s)<len) || (length(grep("\\+",s))==0)) return(s)
  ##     ## cat("abc\n")
  ##     lens <- cumsum(sapply(strsplit(s,"\\+")[[1]],nchar)+1)
  ##     paste(substr(s,1,max(lens[lens<len])-1),"+...",sep="")
  ##   }
  ## WRAP here
  heading <- paste("Likelihood Ratio Tests",
                   paste(mterms,
                         collapse="\n"),
                   sep="\n")
  ltab <- cbind(ltab,Chisq=abs(c(NA,diff(ltab[,"Deviance"]))),
                Df=abs(c(NA,diff(ltab[,"Tot Df"]))))
  ltab <- cbind(ltab,"Pr(>Chisq)"=c(NA,pchisq(ltab[,"Chisq"][-1],
                       ltab[,"Df"][-1],lower.tail=FALSE)))
  rownames(ltab) <- 1:nrow(ltab)
  attr(ltab,"heading") <- heading
  class(ltab) <- "anova"
  ltab
})

## translate from profile to data frame, as either
## S3 or S4 method
as.data.frame.profile.mle2 <- function(x, row.names = NULL,
                                       optional = FALSE, ...) {
    m1 <- mapply(function(vals,parname) {
        ## need to unpack the vals data frame so that
        ## parameter names show up properly
        do.call("data.frame",
                c(list(param=rep(parname,nrow(vals))),
                  as.list(vals),focal=list(vals$par.vals[,parname])))
            },
                 x@profile,
                 as.list(names(x@profile)),
                 SIMPLIFY=FALSE)
    m2 <- do.call("rbind",m1)
    m2
}

setAs("profile.mle2","data.frame",
      function(from) {
          as.data.frame.profile.mle2(from)
          })


BIC.mle2 <- stats4:::BIC
