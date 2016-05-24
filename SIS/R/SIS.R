SIS <- function(x, y, family = c("gaussian","binomial","poisson","cox"), penalty=c("SCAD","MCP","lasso"), concavity.parameter = switch(penalty, SCAD=3.7, 3), tune = c("bic","ebic","aic","cv"), nfolds = 10, 
                type.measure = c("deviance","class","auc","mse","mae"), gamma.ebic = 1, nsis = NULL, iter = TRUE, iter.max = ifelse(greedy==FALSE, 10, floor(nrow(x)/log(nrow(x)))), 
                varISIS = c("vanilla","aggr","cons"), perm = FALSE, q = 1, greedy = FALSE, greedy.size = 1, seed = 0, standardize = TRUE){
  
    this.call=match.call()
    family = match.arg(family)
    penalty = match.arg(penalty)
    tune = match.arg(tune)
    type.measure = match.arg(type.measure)
    varISIS = match.arg(varISIS)
                
    if(is.null(x)||is.null(y)) stop("The data is missing!")  
    if(class(concavity.parameter) != "numeric") stop("concavity.parameter must be numeric!")   
    if(class(nfolds) != "numeric") stop("nfolds must be numeric!")
    if(class(seed) != "numeric") stop("seed must be numeric!")
       
    if(family == "cox" && penalty%in%c("SCAD","MCP"))
       stop("Cox model currently not implemented with selected penalty")
            
    if(type.measure%in%c("class","auc") && family%in%c("gaussian","poisson","cox"))
       stop("'class' and 'auc' type measures are only available for logistic regression")       
       
    if(type.measure%in%c("class","auc","mse","mae") && penalty%in%c("SCAD","MCP"))
       stop("Only 'deviance' is available as type.measure for non-convex penalties")              
    
    fit = switch(family,
          "gaussian"=sisglm(x, y, "gaussian", penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize),
          "binomial"=sisglm(x, y, "binomial", penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize),
          "poisson"=sisglm(x, y, "poisson", penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize),
          "cox"=sisglm(x, y, "cox", penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize)
          )
    fit$call=this.call
    class(fit)=c(class(fit),"SIS")
    return(fit)
}

sisglm <- function(x, y, family, penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize, s1=NULL, s2=NULL, split.tries=0){
  
  storage.mode(x) = "numeric"
  n = dim(x)[1]; p = dim(x)[2]
  models = vector("list")
  if(is.null(nsis)==TRUE) nsis = calculate.nsis(family=family, varISIS=varISIS, n=n, p=p)
  if(is.null(s1)==TRUE){
     set.seed(seed)
     split.sample = sample(1:n); s1 = split.sample[1:ceiling(n/2)]; s2 = setdiff(split.sample,s1)
  }
  if(standardize == TRUE){ 
     old.x = x
     x = scale(x)
  }
  iterind = 0 
  
  if(iter == TRUE){
     ix0 = sort(obtain.ix0(x=x, y=y, s1=s1, s2=s2, family=family, nsis=nsis, iter=iter, varISIS=varISIS, perm=perm, q=q, greedy=greedy, greedy.size=greedy.size, iterind=iterind))
     repeat{
        iterind = iterind + 1
        cat("Iter ", iterind, ": ", ix0,"\n")
        coef.beta = obtain.beta(x, y, ix0, family, penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic)
        ix1 = sort(ix0[which(abs(coef.beta)>1e-10)])
        if(length(ix1) == 0){
           split.tries = split.tries + 1
           split.sample = sample(1:n); s1 = split.sample[1:ceiling(n/2)]; s2 = setdiff(split.sample,s1)
           cat("Sample splitting attempt: ", split.tries, "\n")
           if(split.tries == 20){ 
              varISIS = "vanilla"; perm = TRUE; greedy = FALSE; tune = "cv"
              cat("No variables remaining after ", split.tries," sample splitting attempts! \n")
              cat("Trying a more conservative variable screening approach with a data-driven threshold for marginal screening! \n")
              return(sisglm(x, y, family, penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize, s1=NULL, s2=NULL))
           }else 
              return(sisglm(x, y, family, penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize, s1, s2, split.tries))
        }      
        cat("Iter ", iterind, ": ix1:", ix1,"\n")
        if(length(ix1) >= nsis || iterind >= iter.max){
           ix0 = ix1
           if(length(ix1) >= nsis)
              cat("Maximum number of variables selected \n")
           if(iterind >= iter.max)
              cat("Maximum number of iterations reached \n")              
           break
        }

        models[[iterind]] = ix1; flag.models = 0
        if(iterind > 1){
           for(j in 1:(iterind-1)){
               if(identical(models[[j]],ix1) == TRUE) flag.models = 1
           }
        }
        if(flag.models==1){
           ix0 = ix1
           cat("Model already selected \n")
           break   
        }     
        
        candind = setdiff(1:p, ix1)
        pleft = nsis - length(ix1)
        newix = sort(obtain.newix(x=x, y=y, candind=candind, ix1=ix1, s1=s1, s2=s2, family=family, pleft=pleft, varISIS=varISIS, perm=perm, q=q, greedy=greedy, greedy.size=greedy.size, iterind=iterind))
        cat("Iter ", iterind, ": newix:", newix,"\n")
        ix1 = sort(c(ix1, newix))
        if(setequal(ix1,ix0)){
           ix0 = sort(setdiff(ix1,newix))
           cat("Model already selected \n")
           break
        }
        ix0 = ix1
     }                        # end repeat
 
  }else{                      # end if(iter==TRUE)
     ix0 = sort(obtain.ix0(x=x, y=y, s1=s1, s2=s2, family=family, nsis=nsis, iter=iter, varISIS=varISIS, perm=perm, q=q, greedy=greedy, greedy.size=greedy.size, iterind=iterind))     
     coef.beta = obtain.beta(x, y, ix0, family, penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic)
     ix0 = sort(ix0[which(abs(coef.beta)>1e-10)])   
  }    
   
  # Unstandardize for final fit      
  if(standardize == TRUE) x.final = old.x[,ix0]
  else x.final = x[,ix0]
  final = final.fit(as.matrix(x.final), y, family, penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic)
  if(tune == "cv") final.pre = as.vector(coef(final$fit, s="lambda.min"))
  else final.pre = as.vector(coef(final$fit)[,final$ind])
  if(family== "cox") names(final.pre) = paste("V", ix0, sep="") 
  else names(final.pre) = c("(Intercept)", paste("V", ix0, sep=""))
  return(list(ix=ix0, coef.est=final.pre, fit=final$fit, path.index=final$ind))
}
