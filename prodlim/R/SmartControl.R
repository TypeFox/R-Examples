# {{{ SmartControl 
#' Function to facilitate the control of arguments passed to subroutines.
#' 
#' Many R functions need to pass several arguments to several different
#' subroutines. Such arguments can are given as part of the three magic dots
#' "...". The function SmartControl reads the dots together with a list of
#' default values and returns for each subroutine a list of arguments.
#' 
#' 
#' @param call A list of named arguments, as for example can be obtained via
#' \code{list(...)}.
#' @param keys A vector of names of subroutines.
#' @param ignore A list of names which are removed from the argument
#' \code{call} before processing.
#' @param defaults A named list of default argument lists for the subroutines.
#' @param forced A named list of forced arguments for the subroutines.
#' @param split Regular expression used for splitting keys from arguments.
#' Default is \code{"\."}.
#' @param ignore.case If \code{TRUE} then all matching and splitting is not
#' case sensitive.
#' @param replaceDefaults If \code{TRUE} default arguments are replaced by
#' given arguments. Can also be a named list with entries for each subroutine.
#' @param verbose If \code{TRUE} warning messages are given for arguments in
#' \code{call} that are not ignored via argument \code{ignore} and that do not
#' match any \code{key}.
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{plot.prodlim}}
#' @keywords Graphics
#' @examples
#' 
#' 
#' myPlot = function(...){
#'  ## set defaults
#'  plot.DefaultArgs=list(x=0,y=0,type="n")
#'  lines.DefaultArgs=list(x=1:10,lwd=3)
#'  ## apply smartcontrol
#'  x=SmartControl(call=list(...),
#'            defaults=list("plot"=plot.DefaultArgs, "lines"=lines.DefaultArgs),
#'         ignore.case=TRUE,keys=c("plot","axis2","lines"),
#'              forced=list("plot"=list(axes=FALSE),"axis2"=list(side=2)))
#'  ## call subroutines
#'  do.call("plot",x$plot)
#'  do.call("lines",x$lines)
#'  do.call("axis",x$axis2)
#' }
#' myPlot(plot.ylim=c(0,5),plot.xlim=c(0,20),lines.lty=3,axis2.At=c(0,3,4))
#'
#' @export
SmartControl <- function(call,
                         keys,
                         ignore,
                         defaults,
                         forced,
                         split,
                         ignore.case=TRUE,
                         replaceDefaults,
                         verbose=TRUE)
  # }}}
{
  if (missing(split)) split <- "\\."
  # {{{ set up argument list 
  SmartArgs <- as.list(call)
  SmartArgs <- SmartArgs[names(SmartArgs)!=""]
  if (ignore.case==TRUE){
    names(SmartArgs) <- tolower(names(SmartArgs))
  }
  # }}}
  # {{{remove ignorable arguments
  if (!missing(ignore) && is.character(ignore)){
    if (ignore.case==TRUE){
      ignore <- tolower(ignore)
    }
    SmartArgs <- SmartArgs[match(names(SmartArgs),
                                 ignore,
                                 nomatch=0)==0]
  }
  if (verbose==TRUE){
    allKeysRegexp <- paste("^",keys,split,sep="",collapse="|")
    notIgnored <- grep(allKeysRegexp,names(SmartArgs),value=FALSE,ignore.case=TRUE)
    Ignored <- names(SmartArgs)[-notIgnored]
    SmartArgs <- SmartArgs[notIgnored]
    if (length(Ignored)>0){
      paste(Ignored,collapse=", ")
      warning(paste("The following argument(s) are not smart and therefore ignored: ",paste(Ignored,collapse=", ")))
    }
  }
  # }}}
  # {{{ default arguments
  DefaultArgs <- vector(mode="list",length=length(keys))
  names(DefaultArgs) <- keys
  if (!missing(defaults)){
    whereDefault <- match(names(defaults),names(DefaultArgs),nomatch=0)
    if (all(whereDefault))
      DefaultArgs[whereDefault] <- defaults
    else
      stop("Could not find the following default arguments: ",paste(names(defaults[0==whereDefault]),","))
  }
  if (!missing(replaceDefaults)){
    if (length(replaceDefaults)==1){
      replaceDefaults <- rep(replaceDefaults,length(keys))
      names(replaceDefaults) <- keys
    }
    else {
      stopifnot(length(replaceDefaults)==length(keys))
      stopifnot(all(match(names(replaceDefaults),keys)))
      replaceDefaults <- replaceDefaults[keys]
    }
  }
  else{
    replaceDefaults <- rep(FALSE,length(keys))
    names(replaceDefaults) <- keys
  }
  # }}}
  # {{{ forced arguments
  keyForced <- vector(mode="list",length=length(keys))
  names(keyForced) <- keys
  if (!missing(forced)){
    whereDefault <- match(names(forced),names(keyForced),nomatch=0)
    if (all(whereDefault))
      keyForced[whereDefault] <- forced
    else stop("Not all forced arguments found.")
  }
  # }}}
  # {{{ loop over keys
  keyArgList <- lapply(keys,function(k){
    keyRegexp <- paste("^",k,split,sep="")
    foundArgs <- grep(keyRegexp,names(SmartArgs),value=TRUE,ignore.case=TRUE)
    if (length(foundArgs)>0){
      keyArgs <- SmartArgs[foundArgs]
      if (ignore.case)
        argNames <- sapply(strsplit(tolower(names(keyArgs)),tolower(keyRegexp)),function(x)x[[2]])
      else
        argNames <- sapply(strsplit(names(keyArgs),keyRegexp),function(x)x[[2]])
      
      keyArgs <- lapply(keyArgs,function(x){
        ## expressions for arrow labels in plot.Hist
        ## cannot be evaluated at this point
        ## if the expression is communicated
        ## more than one level higher
        maybeFail <- try(e <- eval(x),silent=TRUE)
        if (class(maybeFail)=="try-error")
          x
        else
          eval(x)
      })
      names(keyArgs) <- argNames
    }
    else{
      keyArgs <- NULL
    }
    # }}}
    # {{{ prepending the forced arguments-----------------
    if (length(keyForced[[k]])>0){
      keyArgs <- c(keyForced[[k]],keyArgs)
    }
    # }}}
    # {{{ appending default arguments
    if (length(DefaultArgs[[k]])>0 && replaceDefaults[k]==FALSE){
      keyArgs <- c(keyArgs,DefaultArgs[[k]])
    }
    # }}}
    # {{{ removing duplicates
    if (!is.null(names(keyArgs))){
      keyArgs[!duplicated(names(keyArgs))]
    }
  })
  
  names(keyArgList) <- keys
  keyArgList
  # }}}
}
