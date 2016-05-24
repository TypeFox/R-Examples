#' The fence procedure for generalised linear models
#'
#' This function implements the fence procedure to
#' find the best generalised linear model.
#'
#' @param mf an object of class \code{\link[stats]{glm}}
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
#' @param ... further arguments (currently unused)
#' @seealso \code{\link{af}}, \code{\link{lmfence}}
#' @references Jiming Jiang, Thuan Nguyen, J. Sunil Rao,
#'   A simplified adaptive fence procedure, Statistics &
#'   Probability Letters, Volume 79, Issue 5, 1 March 2009,
#'   Pages 625-629, http://dx.doi.org/10.1016/j.spl.2008.10.014.
#' @noRd
#' @family fence

glmfence = function(mf,
                    cstar,
                    nvmax,
                    adaptive=TRUE,
                    trace=TRUE,...){
  method="ML"
  if(any(class(mf)=="glm")!=TRUE){
    stop("The argument to mf needs to be a glm object.")
  }
  if(attr(mf$terms,"intercept")==0){
    stop("Please allow for an intercept in your model.")
  }
  m = mextract(mf)
  kf = m$k
  fixed = m$fixed
  family = m$family
  yname = m$yname
  Xy = m$X
  n = m$n
  wts = m$wts
  if(missing(nvmax)) nvmax=kf
  null.ff = stats::as.formula(paste(yname,"~1")) # null formula
  m0 = stats::glm(null.ff, data = Xy, family=family, weights = wts) # null model
  Qmf = Qm(mf, method=method) # Qm for the full model
  Qm0 = Qm(m0, method=method) # Qm for the null model
  ret = met = list()

  # Null model
  if(trace) cat(paste("Null model "))
  UB = Qmf + cstar*sigMM(k.mod=1, method, k.full=kf,adaptive=adaptive)
  if(Qm0<=UB){
    if(trace) txt.fn(Qm0,UB,m0)
    ret[[1]] = null.ff
    return(ret)
  } else if(trace) cat("(Not a candidate model) \n")

  if(cstar<5){ # avoids having to add variables to get the full model
    nvmax = kf
    prev.nvmax = nvmax
  } else if(nvmax<5){
    prev.nvmax = nvmax
    nvmax = nvmax+5
  } else prev.nvmax = nvmax
  # look around for the best model at each model size
  while(prev.nvmax<=kf){
    prev.nvmax = nvmax
    bg = bestglm::bestglm(Xy=Xy, family=family,
                 IC = "BIC",
                 TopModels = 5*kf,
                 nvmax = nvmax, weights=wts)
    lc = bg$Subsets[,1:kf]+0 # 'leaps' candidates
    for(i in 2:nvmax){
      if(trace) cat(paste("Model size:",i,""))
      UB = Qmf + cstar*sigMM(k.mod=i, method, k.full=kf, adaptive=adaptive)
      mnames = colnames(lc)[which(lc[i,]==1)]
      ff = stats::as.formula(paste(yname," ~ ",paste(mnames[-1],collapse="+"),sep=""))
      em = stats::glm(formula=ff, data=Xy, family=family, weights=wts)
      hatQm = Qm(em,method=method)
      if(hatQm<=UB){
        if(trace){
          cat("\n")
          cat("Candidate model found via bestglm. \n")
          cat("Exploring other options at this model size. ")
          txt.fn(hatQm,UB,em)
        }
        pos=1
        environment(ff) = globalenv()
        ret[[pos]] = ff # record the result
        met[[pos]] = hatQm #record its score
        # look for others at this model size:
        lfm = bg$BestModels[,1:(kf-1)]+0
        lfm.sum = apply(lfm,1,sum)
        lfm = lfm[lfm.sum==i-1,]
        # remove already estimated model from lfm:
        check.fn = function(x) !all(x==lc[i,-1])
        lfm = lfm[apply(lfm,1,check.fn),]
        if(dim(lfm)[1]>0){
          for(j in 1:dim(lfm)[1]){
            mnames = colnames(lfm)[which(lfm[j,]==1)]
            ff = stats::as.formula(paste(yname," ~ ",
                                  paste(mnames,collapse="+"),
                                  sep=""))
            em = stats::glm(ff, data=Xy, family=family, weights=wts)
            hatQm = Qm(em,method=method)
            if(hatQm<=UB){
              if(trace) txt.fn(hatQm,UB,em)
              pos = pos+1
              environment(ff) = globalenv()
              ret[[pos]] = ff
              met[[pos]] = hatQm
            }
          }
          return(ret[rank(unlist(met),ties.method="random")])
        } else return(ret)
      }
      if(trace) cat("(No candidate models found) \n")
    }
    if(trace) cat(" (No candidate models found: increasing nvmax) \n", cstar)
    nvmax = nvmax+5
  }
}
