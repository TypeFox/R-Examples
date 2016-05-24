#' Wild Cluster Bootstrapped p-Values For For Regression With Instrumental Variables
#'
#' This software estimates p-values using wild cluster bootstrapped t-statistics for instrumental variables regression models (Cameron, Gelbach, and Miller 2008). Residuals are repeatedly re-sampled by cluster to form a pseudo-dependent variable, a model is estimated for each re-sampled data set, and inference is based on the sampling distribution of the pivotal (t) statistic. Users may choose whether to impose the null hypothesis for independent variables; the null is never imposed for the intercept or any model that includes factor variables. Confidence intervals are only reported when the null hypothesis is \emph{not} imposed.
#'
#' @param mod A linear (identity link) model estimated using \code{ivreg}.
#' @param dat The data set used to estimate \code{mod}.
#' @param cluster A formula of the clustering variable.
#' @param ci.level What confidence level should CIs reflect? (Note: only reported when \code{impose.null == FALSE}).
#' @param impose.null Should we impose the null Ho?
#' @param boot.reps The number of bootstrap samples to draw.
#' @param report Should a table of results be printed to the console?
#' @param prog.bar Show a progress bar of the bootstrap (= TRUE) or not (= FALSE).
#'
#' @return A list with the elements
#' \item{p.values}{A matrix of the estimated p-values.}
#' \item{ci}{A matrix of confidence intervals (if null not imposed).}
#' @author Justin Esarey
#' @note Code to estimate clustered standard errors by Mahmood Arai: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/. Cluster SE degrees of freedom correction = (M/(M-1)) with M = the number of clusters.
#' @examples
#' \dontrun{
#' 
#' # example one: predict cigarette consumption
#' data("CigarettesSW", package = "AER") 
#' CigarettesSW$rprice <- with(CigarettesSW, price/cpi)
#' CigarettesSW$rincome <- with(CigarettesSW, income/population/cpi)
#' CigarettesSW$tdiff <- with(CigarettesSW, (taxs - tax)/cpi)
#' fm <- ivreg(log(packs) ~ log(rprice) + log(rincome) | 
#'     log(rincome) + tdiff + I(tax/cpi), data = CigarettesSW)
#' 
#' # compute cluster-adjusted p-values
#' cluster.wd.c <- cluster.wild.ivreg(fm, dat=CigarettesSW, cluster = ~state, report = T)
#' 
#' 
#' 
#' # example two: pooled IV analysis of employment
#' require(plm)
#' require(AER)
#' data(EmplUK)
#' EmplUK$lag.wage <- lag(EmplUK$wage)
#' emp.iv <- ivreg(emp ~ wage + log(capital+1) | output + lag.wage + log(capital+1), data = EmplUK)
#' 
#' # compute cluster-adjusted p-values
#' cluster.wd.e <- cluster.wild.ivreg(mod=emp.iv, dat=EmplUK, cluster = ~firm)
#' 
#' }
#' @rdname cluster.wild.ivreg
#' @import AER
#' @import Formula
#' @importFrom lmtest coeftest
#' @importFrom sandwich estfun
#' @importFrom sandwich sandwich
#' @references Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. 2008. "Bootstrap-Based Improvements for Inference with Clustered Errors." \emph{The Review of Economics and Statistics} 90(3): 414-427.
#' @import stats
#' @importFrom utils write.table
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @export
#' 
#

cluster.wild.ivreg<-function(mod, dat, cluster, ci.level = 0.95, impose.null = TRUE, boot.reps = 1000, report = TRUE, prog.bar = TRUE){
  
  form <- mod$formula                                            # what is the formula of this model?
  variables <- all.vars(form)                                    # what variables are in this model?
  clust.name <- all.vars(cluster)                                # what is the name of the clustering variable?
  used.idx <- which(rownames(dat) %in% rownames(mod$model))      # what were the actively used observations in the model?
  dat <- dat[used.idx,]                                          # keep only active observations (drop the missing)
  clust <- as.vector(unlist(dat[[clust.name]]))                  # store cluster index in convenient vector
  G<-length(unique(clust))                                       # how many clusters are in this model?
  "%w/o%" <- function(x, y) x[!x %in% y]                         # a little function to create a without function (see ?match)
  ind.variables <- names(mod$coefficients) %w/o% "(Intercept)"   # what independent variables are in this data set?
  dv <- variables %w/o% all.vars(update(form, 1 ~ .))            # what is the dependent variable?
  ind.variables.data <- all.vars(update(form, 1 ~ .))            # gives the names of IVs in the data set (before function transforms)
  ind.variables.names <- names(coefficients(mod))                # printed names of all fitted coefficients (incld. intercept)
  
  # in case dv is wrapped in a function, need to set it to its functional value
  # so that residuals can be added w/o incident
  dat$dv.new <- mod$y                                            # add transformed DV into data set
  form.new <- update(as.Formula(form), dv.new ~ .)               # formula subbing in transformed DV
  

  # check to see whether any IVs are factors
  fac <- c()
  for(i in 1:length(ind.variables.data)){
    fac[i] <- is.factor(dat[,ind.variables.data[i]])
  }
  fac <- max(fac)
  
  # do not impose the null for factor variables
  if(fac == 1 & impose.null == TRUE){
    cat("\n","\n", "Note: null not imposed (factor variables are present).", "\n")
    impose.null<-FALSE
  }
    
  
  # load in a function to create clustered standard errors
  # by Mahmood Arai: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/
  cl   <- function(dat, fm, cluster){
    #require(sandwich, quietly = TRUE)
    #require(lmtest, quietly = TRUE)
    M <- length(unique(cluster))
    N <- length(cluster)
    K <- fm$rank
    dfc <- (M/(M-1))
    uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
    vcovCL <- dfc*sandwich(fm, meat.=crossprod(uj)/N)
    coeftest(fm, vcovCL) }
  
  
  se.clust <- cl(dat, mod, clust)[ind.variables.names,2]               # retrieve the clustered SEs
  beta.mod <- coefficients(mod)[ind.variables.names]                   # retrieve the estimated coefficients
  w <- beta.mod / se.clust                                             # calculate the wald test statistics
 
  
  # if the null is to be imposed, execute the following code
  if(impose.null==TRUE){
    
    p.store <- c()                                                            # store wild boostrapped p-values
    w.store <- matrix(data=NA, nrow=boot.reps, ncol=length(ind.variables))    # store bootstrapped test statistics
    
    if(names(mod$coefficients)[1] == "(Intercept)"){offset <- 1}else{offset <- 0}
    
    if(prog.bar==TRUE){cat("\n")}
    for(j in 1:length(ind.variables)){
      
      if(prog.bar==TRUE){cat("Independent variable being bootstrapped: ", ind.variables[j], "\n")}
      
      # run model imposing the null hypothesis
      mod.form <- as.Formula(mod$form)
      form.null <- update(mod.form, as.formula(paste("~ . -", ind.variables[j], "| . -", ind.variables[j] ) ))
      mod.call <- mod$call
      mod.call[[2]] <- form.null              # change formula to impose the null hypothesis
      mod.call[[3]] <- quote(dat)             # use data set with dropped NA observations
      mod.null <- eval(mod.call)
      null.resid <- residuals(mod.null)
      
      boot.dat <- dat           # copy the data set into a bootstrap resampling dataset
      wald.store <- c()         # create a container for storing the test statistics    
      
      if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
      for(i in 1:boot.reps){
        
        if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
        
        weight <- c(1, -1)[rbinom(G, size=1, prob=0.5)                       # assign wild bootstrap weights
                           + 1][match(clust, unique(clust))]      
        pseudo.resid <- null.resid*weight                                    # create pseudo-residuals using weights
        pseudo.dv <- predict(mod.null)+ pseudo.resid                         # create pseudo-observations using pseudo-residuals
        boot.dat[,"dv.new"] <- pseudo.dv                                     # create a bootstrap replicate data set (note use of dv.new)
        
        boot.mod.call <- mod$call                                            # retrieve the original model call
        boot.mod.call[[2]] <- quote(form.new)                                # use formula with transformed DV
        boot.mod.call[[3]] <- quote(boot.dat)                                # use bootstrap data set
        boot.mod <- suppressWarnings(tryCatch(eval(boot.mod.call),           # attempt to re-run the model
                            error = function(e){return(NA)})) 
        
        se.boot <- cl(boot.dat, boot.mod, clust)[offset+j,2]                  # retrieve the bootstrap clustered SE
        beta.boot <- coefficients(boot.mod)[offset+j]                         # store the bootstrap beta coefficient
        wald.store[i] <- beta.boot / se.boot                                  # store the bootstrap test statistic
        
      }
      if(prog.bar==TRUE){close(pb)}
      
      p.store[j] <- 1 - ( sum( abs(w[offset + j]) > abs(wald.store) ) / boot.reps )    # calculate the wild bootstrap p-value
      w.store[,j] <- wald.store
       
    }
    
    # calculate t-stat for intercept, if present, w/o imposing the null
    if(names(mod$coefficients)[1] == "(Intercept)"){
      
      if(prog.bar==TRUE){cat("Independent variable being bootstrapped:  Intercept (null not imposed)", "\n")}
  
      # don't impose the null for the constant (but still call it null.resid)
      null.resid <- residuals(mod)
      
      boot.dat <- dat           # copy the data set into a bootstrap resampling dataset
      wald.store <- c()         # create a container for storing the test statistics    
      
      if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
      for(i in 1:boot.reps){
        
        if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
        
        weight <- c(1, -1)[rbinom(G, size=1, prob=0.5)                 # assign wild bootstrap weights
                           + 1][match(clust, unique(clust))]    
        pseudo.resid <- null.resid*weight                              # create pseudo-residuals using weights
        pseudo.dv <- predict(mod)+ pseudo.resid                        # create pseudo-observations using pseudo-residuals
        boot.dat[,"dv.new"] <- pseudo.dv                               # create a bootstrap replicate data set (note use of dv.new)
        
        boot.mod.call <- mod$call                                      # retrieve the original model call
        boot.mod.call[[2]] <- quote(form.new)                          # use formula with transformed DV
        boot.mod.call[[3]] <- quote(boot.dat)                          # use bootstrap data set
        boot.mod <- suppressWarnings(tryCatch(eval(boot.mod.call),     # attempt to re-run the model
                            error = function(e){return(NA)})) 
        
        se.boot <- cl(boot.dat, boot.mod, clust)[1,2]                  # retrieve the bootstrap clustered SE
        beta.boot <- coefficients(boot.mod)[1]                         # store the bootstrap beta coefficient
        wald.store[i] <- (beta.boot - beta.mod[1]) / se.boot           # store the bootstrap test statistic
        
      }
      if(prog.bar==TRUE){close(pb)}
      
      p.store <- c( 1 - ( sum( abs(w[1]) > abs(wald.store) ) / boot.reps ), p.store)    # calculate the wild bootstrap p-value
      w.store <- cbind(wald.store, w.store)
      
    
    }

    ci.lo = NULL
    ci.hi = NULL
    print.ci = NULL
    out.ci = NULL
    
  # if the null is NOT to be imposed...
  }else{
    
    if(prog.bar==TRUE){cat("Wild Cluster bootstrapping w/o imposing null...", "\n")}

    boot.dat <- dat                                                                 # copy the data set into a bootstrap resampling dataset
    w.store <- matrix(data=NA, nrow=boot.reps, ncol=length(ind.variables.names))    # store bootstrapped test statistics
    
    resid <- residuals(mod)                                                         # get the residuals for the model
    
    if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
    for(i in 1:boot.reps){
      
      if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
      
      weight <- c(1, -1)[rbinom(G, size=1, prob=0.5)                         # assign wild bootstrap weights
                         + 1][match(clust, unique(clust))]                
      pseudo.resid <- resid*weight                                           # create pseudo-residuals using weights
      pseudo.dv <- predict(mod)+ pseudo.resid                                # create pseudo-observations using pseudo-residuals
      boot.dat[,"dv.new"] <- pseudo.dv                                       # create a bootstrap replicate data set (note use of dv.new)
      
      boot.mod.call <- mod$call                                              # retrieve the original model call
      boot.mod.call[[2]] <- quote(form.new)                                  # use formula with transformed DV
      boot.mod.call[[3]] <- quote(boot.dat)                                  # use bootstrap data set
      boot.mod <- suppressWarnings(tryCatch(eval(boot.mod.call),             # attempt to re-run the model
                          error = function(e){return(NA)})) 
      
      se.boot <- cl(boot.dat, boot.mod, clust)[,2]                           # retrieve the bootstrap clustered SE
      beta.boot <- coefficients(boot.mod)                                    # store the bootstrap beta coefficient
      w.store[i,] <- (beta.boot-beta.mod) / se.boot                          # store the bootstrap test statistic
      
    }
    if(prog.bar==TRUE){close(pb)}
    
    comp.fun<-function(vec2, vec1){as.numeric(vec1>vec2)}                              # a simple function comparing v1 to v2
    p.store.s <- t(apply(X = abs(w.store), FUN=comp.fun, MARGIN = 1, vec1 = abs(w)))   # compare the BS test stats to orig. result
    p.store <- 1 - ( colSums(p.store.s) / dim(w.store)[1] )                            # calculate the cluster bootstrap p-value
    

    # compute critical t-statistics for CIs
    crit.t <- apply(X=abs(w.store), MARGIN=2, FUN=quantile, probs=ci.level )
    ci.lo <- beta.mod - crit.t*se.clust
    ci.hi <- beta.mod + crit.t*se.clust

    

    print.ci <- cbind(ind.variables.names, ci.lo, ci.hi)
    print.ci <- rbind(c("variable name", "CI lower", "CI higher"), print.ci)
    
    out.ci <- cbind(ci.lo, ci.hi)
    rownames(out.ci) <- ind.variables.names
    colnames(out.ci) <- c("CI lower", "CI higher")
        
  }
  
  out <- matrix(p.store, ncol=1)
  colnames(out) <- c("wild cluster BS p-value")
  rownames(out) <- ind.variables.names
  out.p <- cbind(ind.variables.names, round(out, 3))
  out.p <- rbind(c("variable name", "wild cluster BS p-value"), out.p)
  
  printmat <- function(m){
    write.table(format(m, justify="right"), row.names=F, col.names=F, quote=F, sep = "   ")
  }
  
  if(report==T){
    
    cat("\n", "\n", "Wild Cluster Bootstrapped p-values: ", "\n", "\n")
    printmat(out.p)
    if(is.null(print.ci) == FALSE){
      cat("\n", "Confidence Intervals (derived from bootstrapped t-statistics): ", "\n", "\n")
      printmat(print.ci)
    }

   }
  
  out.list<-list()
  out.list[["p.values"]]<-out
  out.list[["ci"]] <- out.ci
  return(invisible(out.list))
  
  
   
}

