#' Converts an \code{evalGrid} object into a \code{data.frame}
#'
#'  Converts the results contained in the 
#'  object returned by \code{\link{evalGrids}} 
#'  into a \code{data.frame}. If the results can not
#'  be coerced automatically into a \code{data.frame}, the
#'  user can provide a function to pre-process the
#'  results (see \code{convert.result.fun}). Furthermore, 
#'  univariate functions to summarize the results over
#'  the replications can be specified via \code{summary.fun}.
#'
#'
#'@param x  an object returned by 
#'  \code{\link{evalGrids}}
#'@param \dots  only for S3 method consistency
#'@param convert.result.fun  a functions that converts the result
#'  object contained in \code{x} into a \code{data.frame}
#'@param summary.fun  univariate functions to summarize the results (numeric or logical) over
#'  the replications, e.g. mean, sd. Alternatively, \code{summary.fun} can be one
#'  function that may return a vector.
#'@param progress if \code{TRUE} a progress bar is shown in the console.
#'@return  a \code{data.frame} with the parameter constellations
#'  for the data generation and evaluation and the results
#'  (probably summarized).
#'@author  Marsel Scheer
#'@seealso  \code{\link{evalGrids}}
#'@examples
#'
#'genRegData <- function(){
#'  data.frame(
#'      x = 1:10,
#'      y = rnorm(10, mean=1:10))
#'}
#'
#'eg <- evalGrids(
#'  expandGrid(fun="genRegData"),
#'  expandGrid(proc="lm", formula=c("y ~ x", "y ~ x + I(x^2)")),
#'  replications=5)
#'
#'lm2df = function(lm.object) {
#'  ret = coef(summary.lm(lm.object))[, 1:2]
#'  data.frame(covariable = rownames(ret), ret, check.names=FALSE)
#'}
#'as.data.frame(eg, convert.result.fun=lm2df, progress=TRUE)
#'as.data.frame(eg, convert.result.fun=lm2df, summary.fun=c(mean, sd), progress=TRUE)
#'@import plyr
#'@export
as.data.frame.evalGrid <-
  function(x, ..., convert.result.fun = identity, summary.fun=NULL, progress=FALSE) {
    postFun = NULL

    ellipsis = list(...)
    if (is.element("post.proc", names(ellipsis))) stop("post.proc is deprecated. Please use summary.fun")    
    if (is.element("value.fun", names(ellipsis))) stop("value.fun is deprecated. Please use convert.result.fun")
    
    if (!is.null(summary.fun)){      
      if (length(summary.fun) == 1) {
        postFun = summary.fun
      } else {
        postFun = do.call(funstofun, as.list(match.call()$summary.fun[-1]))    
      }
    }

    .progress="none"
    if (progress)
      .progress="text"    

    with(x, {
      ddply(expandGrid(i=1:nrow(dataGrid), j=1:nrow(procGrid)), .(i, j), function(row){
        i=row[1,1]
        j=row[1,2]  

        ret = ldply(simulation[[i]], function(rep) convert.result.fun(rep$results[[j]]))
        if (nrow(ret) > 0)
          ret = cbind(replication=gl(length(simulation[[i]]), nrow(ret)/length(simulation[[i]])), ret)
        
        
        if (nrow(ret) > 0){
          if (!is.null(postFun)){            
            ret$replication=NULL
            idx = which(sapply(1:ncol(ret), function(i) all(is.numeric(ret[,i]) | is.logical(ret[,i]))))
            if (length(idx) == 0)
              stop("Only numeric or logical variables passed to in summary.fun. But the results does not seem to have numeric or logical variables.")
            mdf = melt(ret, measure.vars=idx)
            ret = cast(mdf, ... ~ variable, postFun)
          }
        } else {
          ret = data.frame(.evalGridComment="Results missing")
        }
        if (!is.null(x$summary.fun))
          ret$replication=NULL
        
        suppressWarnings(cbind(dataGrid[i,,drop=FALSE], procGrid[j,,drop=FALSE], ret))
      }, .progress=.progress)
    })
  }
