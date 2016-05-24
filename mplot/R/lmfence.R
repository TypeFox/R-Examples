#' The fence procedure for linear models
#'
#' This function implements the fence procedure to
#' find the best linear model.
#'
#' @param mf an object of class \code{\link[stats]{lm}}
#'   specifying the full model.
#' @param cstar the boundary of the fence, typically found
#'   through bootstrapping.
#' @param nvmax the maximum number of variables that will be
#'   be considered in the model.
#' @param adaptive logical. If \code{TRUE} the boundary of the fence is
#'   given by cstar.  Otherwise, it the original (non-adaptive) fence
#'   is performed where the boundary is cstar*hat(sigma)_{M,tildeM}.
#' @param trace logical. If \code{TRUE} the function prints out its
#'   progress as it iterates up through the dimensions.
#' @param force.in the names of variables that should be forced
#'   into all estimated models.
#' @param ... further arguments (currently unused)
#' @seealso \code{\link{af}}, \code{\link{glmfence}}
#' @references Jiming Jiang, Thuan Nguyen, J. Sunil Rao,
#'   A simplified adaptive fence procedure, Statistics &
#'   Probability Letters, Volume 79, Issue 5, 1 March 2009,
#'   Pages 625-629, http://dx.doi.org/10.1016/j.spl.2008.10.014.
#' @noRd
#' @family fence
#' @examples
#' n = 40 # sample size
#' beta = c(1,2,3,0,0)
#' K=length(beta)
#' set.seed(198)
#' X = cbind(1,matrix(rnorm(n*(K-1)),ncol=K-1))
#' e = rnorm(n)
#' y = X%*%beta + e
#' dat = data.frame(y,X[,-1])
#' # Non-adaptive approach (not recommended)
#' lm1 = lm(y~.,data=dat)
#' lmfence(lm1,cstar=log(n),adaptive=FALSE)

lmfence = function(mf, cstar,
                   nvmax,
                   adaptive=TRUE,
                   trace=TRUE,
                   force.in=NULL,...){
  method="ML"
    if(class(mf)!="lm"){
    stop("The argument to mf needs to be a lm object.")
  }
  if(attr(mf$terms,"intercept")==0){
    stop("Please allow for an intercept in your model.")
  }
  m = mextract(mf)
  kf = m$k
  fixed = m$fixed
  yname = m$yname
  data = m$X
  n = m$n
  wts = m$wts
  if(missing(nvmax)) nvmax=kf
  null.ff = stats::as.formula(paste(yname,"~1"))
  m0 = stats::lm(null.ff, data = data, weights=m$wts) # null model
  Qmf = Qm(mf, method=method) # Qm for the full model
  Qm0 = Qm(m0, method=method) # Qm for the null model
  ret = met = list()
  # Null model
  if(trace) cat(paste("Null model "))
  UB = Qmf + cstar*sigMM(k.mod = 1, method = method,
                         k.full = kf, adaptive = adaptive)
  if(Qm0<=UB){
    if(trace) txt.fn(Qm0,UB,m0)
    ret[[1]] = null.ff # record the result
    return(ret)
  } else if(trace) cat("(Not a candidate model) \n")

  if(cstar<5){ # avoids having to add variables to get the full model
    nvmax = kf
    prev.nvmax = nvmax
  } else if(nvmax<5){
    prev.nvmax = nvmax
    nvmax = nvmax+5
  } else prev.nvmax = nvmax
  prev.nvmax = min(prev.nvmax,kf)
  # look around for the best model at each model size
  while(prev.nvmax<=kf){
    prev.nvmax = nvmax
    # finds the best candidate for each model size
    rss = do.call(leaps::regsubsets,list(x=fixed,
                             data=data,
                             nbest = 5+kf,
                             nvmax = nvmax,
                             intercept=TRUE,
                             force.in=force.in,
                             really.big=TRUE,
                             weights = m$wts))
    rs = summary(rss)
    rs.which = data.frame(rs$which+0,row.names = NULL)
    rs.k = apply(rs.which,1,sum)
    rs.bic = split(rs$bic,f = rs.k)
    leaps.cands = lapply(split(rs.which,rs.k),FUN = function(x) x[1,])
    # best model of each size to test if it passes the fence
    leaps.cands = do.call(rbind,leaps.cands)
    lc.k = apply(leaps.cands,1,sum)
    # next best models of each size to test if also pass the fence
    other.cands = lapply(split(rs.which,rs.k),FUN = function(x) x[-1,])
    start = lc.k[1] #2+length(force.in)
    for(i in start:nvmax){
      if(trace) cat(paste("Model size:",i,""))
      UB = Qmf + cstar*sigMM(k.mod = i, method = method,
                             k.full = kf, adaptive = adaptive)
      mnames = colnames(leaps.cands)[which(leaps.cands[lc.k==i,]==1)]
      ff = stats::as.formula(paste(yname," ~ ",
                            paste(mnames[-1],collapse="+"),sep=""))
      em = stats::lm(formula=ff, data=data, weights=m$wts)
      hatQm = Qm(em,method=method)
      if(hatQm<=UB){
        if(trace){
          cat("\n Candidate model found via leaps. \n")
          cat("Exploring other options at this model size. ")
          txt.fn(hatQm,UB,em)
        }
        pos=1
        environment(ff) = globalenv()
        ret[[pos]] = ff # record the result
        met[[pos]] = hatQm #record its score
        lmf = other.cands[[paste(i)]]
        if(dim(lmf)[1]>0){
          for(j in 1:dim(lmf)[1]){
            mnames = colnames(lmf)[which(lmf[j,]==1)]
            ff = stats::as.formula(paste(yname," ~ ",
                                  paste(mnames[-1],collapse="+"),
                                  sep=""))
            em = stats::lm(ff, data = data, weights=m$wts)
            hatQm = Qm(em,method=method)
            if(hatQm<=UB){
              if(trace) txt.fn(hatQm,UB,em)
              pos = pos+1
              environment(ff) = globalenv()
              ret[[pos]] = ff
              met[[pos]] = hatQm
            } else break
          }
          return(ret[rank(unlist(met),ties.method="random")])
        } else return(ret)
      }
      if(trace) cat("(No candidate models found) \n")
    }
    if(trace) cat(" (No candidate models found: increasing nvmax) \n",cstar)
    nvmax = nvmax+5
  }
}
