#' Wild Cluster Bootstrapped p-Values For PLM
#'
#' This software estimates p-values using wild cluster bootstrapped t-statistics for fixed effects panel linear models (Cameron, Gelbach, and Miller 2008). Residuals are repeatedly re-sampled by cluster to form a pseudo-dependent variable, a model is estimated for each re-sampled data set, and inference is based on the sampling distribution of the pivotal (t) statistic. The null is never imposed for PLM models.
#'
#' @param mod A "within" model estimated using \code{plm}.
#' @param dat The data set used to estimate \code{mod}.
#' @param cluster A formula of the clustering variable.
#' @param ci.level What confidence level should CIs reflect? (Note: only reported when \code{impose.null == FALSE}).
#' @param boot.reps The number of bootstrap samples to draw.
#' @param report Should a table of results be printed to the console?
#' @param prog.bar Show a progress bar of the bootstrap (= TRUE) or not (= FALSE).
#'
#' @return A list with the elements
#' \item{p.values}{A matrix of the estimated p-values.}
#' \item{ci}{A matrix of confidence intervals (if null not imposed).}
#' @author Justin Esarey
#' @examples
#' \dontrun{
#' 
#' # predict employment levels, cluster on group
#' require(plm)
#' data(EmplUK)
#'
#' emp.1 <- plm(emp ~ wage + log(capital+1), data = EmplUK, model = "within",
#'              index=c("firm", "year"))
#' cluster.wild.plm(mod=emp.1, dat=EmplUK, cluster="group", ci.level = 0.95,
#'         boot.reps = 1000, report = TRUE, prog.bar = TRUE)
#' 
#' # cluster on time
#' cluster.wild.plm(mod=emp.1, dat=EmplUK, cluster="time", ci.level = 0.95, 
#'             boot.reps = 1000, report = TRUE, prog.bar = TRUE)
#' 
#' }
#' @rdname cluster.wild.plm
#' @import Formula
#' @import plm
#' @references Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. 2008. "Bootstrap-Based Improvements for Inference with Clustered Errors." \emph{The Review of Economics and Statistics} 90(3): 414-427.
#' @import stats
#' @importFrom utils write.table
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @export
#' 

cluster.wild.plm<-function(mod, dat, cluster, ci.level = 0.95, boot.reps = 1000, report = TRUE, prog.bar = TRUE){
  
  if( min( class(dat) != "pdata.frame" ) ){                           # if data not pdata.frame
    dat <- pdata.frame(dat, index=colnames(dat)[1:2])                 # convert it
  }
  
  if(cluster=="group"){                                               
    
    clust <- attr(mod$mod, "index")[,1]                               # which clusters were used?                               
    clust.name <- colnames(attr(mod$mod, "index"))[1]                 # what is the name of the cluster?
    clust.full <- attr(dat, "index")[,1]                              # all clusters in the dataset
    used.idx <- as.numeric(rownames(mod$model))                       # what were the actively used observations in the model?
    G<-length(unique(clust))                                          # how many clusters are in this model?
    
  }
  if(cluster=="time"){
    
    cat("\n", "\n", "Note: specify lags as variables when clustering on time", "\n")
    clust <- attr(mod$mod, "index")[,2]                               # which clusters were used?                               
    clust.name <- colnames(attr(mod$mod, "index"))[2]                 # what is the name of the cluster?
    clust.full <- attr(dat, "index")[,2]                              # all clusters in the dataset
    used.idx <- as.numeric(rownames(mod$model))                       # what were the actively used observations in the model?
    G<-length(unique(clust))                                          # how many clusters are in this model?
    
  }
  if(cluster != "group" & cluster != "time"){
    stop("invalid clustering variable; see help file")
  }
  
  form <- mod$formula                                            # what is the formula of this model?  
  variables <- all.vars(form)                                    # what variables are in this model?
  ind.variables <- names(coefficients(mod))                      # what independent variables are in this model?
  ind.variables.data <- all.vars(update(form, 1 ~ .))            # gives the names of IVs in the data set (before function transforms)
  
  se.clust <- sqrt(diag(vcovHC(mod, cluster=cluster)))           # retrieve the clustered SEs
  beta.mod <- coefficients(mod)[ind.variables]                   # retrieve the estimated coefficients
  w <- beta.mod / se.clust                                       # calculate the t-test statistic  
  
  # in case dv is wrapped in a function, need to set it to its functional value
  # so that residuals can be added w/o incident
  dat$dv.new[used.idx] <- mod$model[,1]                     # add transformed DV into data set
  form.new <- update(form, dv.new ~ .)                      # formula subbing in transformed DV
  
  # check to see whether any IVs are factors
  fac <- c()
  for(i in 1:length(ind.variables.data)){
    fac[i] <- is.factor(dat[,ind.variables.data[i]])
  }
  fac <- max(fac)
  
  if(prog.bar==TRUE){cat("Wild Cluster bootstrapping w/o imposing null...", "\n")}

  boot.dat <- dat                                              # copy the data set into a bootstrap resampling dataset
  w.store <- matrix(data=NA, nrow=boot.reps, 
                ncol=length(ind.variables))                    # store bootstrapped test statistics
  
  resid <- residuals(mod)                                      # get the residuals for the model
  
  if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
  for(i in 1:boot.reps){
    
    if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
    
    weight <- c(1, -1)[rbinom(G, size=1, prob=0.5)                   # assign wild bootstrap weights
                       + 1][match(clust, unique(clust))] 
    pseudo.resid <- resid*weight                                     # create pseudo-residuals using weights
    pseudo.dv <- predict(mod)+ pseudo.resid                          # create pseudo-observations using pseudo-residuals
    boot.dat[used.idx,"dv.new"] <- pseudo.dv                         # create a bootstrap replicate data set (note dv.new)
    
    boot.mod.call <- mod$call                                            # retrieve the original model call
    boot.mod.call[[2]] <- quote(form.new)                                # use formula with transformed DV
    boot.mod.call[[3]] <- quote(boot.dat)                                # use bootstrap data set
    boot.mod <- suppressWarnings(tryCatch(eval(boot.mod.call),           # attempt to re-run the model
                        error = function(e){return(NA)})) 
    
    se.boot <- sqrt(diag(vcovHC(boot.mod, cluster=cluster)))                   # retrieve the bootstrap clustered SE
    beta.boot <- coefficients(boot.mod)[ind.variables]                         # store the bootstrap beta coefficient
    w.store[i,] <- (beta.boot-beta.mod) / se.boot                              # store the bootstrap test statistic
    
    
  }
  if(prog.bar==TRUE){close(pb)}
  
  comp.fun<-function(vec2, vec1){as.numeric(vec1>vec2)}                            # a simple function comparing v1 to v2
  p.store.s <- t(apply(X = abs(w.store), FUN=comp.fun, MARGIN = 1, vec1 = abs(w))) # compare the BS test stats to orig. result
  p.store <- 1 - ( colSums(p.store.s) / dim(w.store)[1] )                          # calculate the cluster bootstrap p-value
  

  # compute critical t-statistics for CIs
  crit.t <- apply(X=abs(w.store), MARGIN=2, FUN=quantile, probs=ci.level )
  ci.lo <- beta.mod - crit.t*se.clust
  ci.hi <- beta.mod + crit.t*se.clust

  

  print.ci <- cbind(ind.variables, ci.lo, ci.hi)
  print.ci <- rbind(c("variable name", "CI lower", "CI higher"), print.ci)
  
  out.ci <- cbind(ci.lo, ci.hi)
  rownames(out.ci) <- ind.variables
  colnames(out.ci) <- c("CI lower", "CI higher")
    
  
  out <- matrix(p.store, ncol=1)
  colnames(out) <- c("wild cluster BS p-value")
  rownames(out) <- ind.variables
  out.p <- cbind(ind.variables, round(out, 3))
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

