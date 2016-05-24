#' Pairs Cluster Bootstrapped p-Values For PLM
#'
#' This software estimates p-values using pairs cluster bootstrapped t-statistics for fixed effects panel linear models (Cameron, Gelbach, and Miller 2008). The data set is repeatedly re-sampled by cluster, a model is estimated, and inference is based on the sampling distribution of the pivotal (t) statistic. 
#'
#' @param mod A "within" model estimated using \code{plm}.
#' @param dat The data set used to estimate \code{mod}.
#' @param cluster Clustering dimension ("group", the default, or "time").
#' @param ci.level What confidence level should CIs reflect?
#' @param boot.reps The number of bootstrap samples to draw.
#' @param cluster.se Use clustered standard errors (= TRUE) or ordinary SEs (= FALSE) for bootstrap replicates.
#' @param report Should a table of results be printed to the console?
#' @param prog.bar Show a progress bar of the bootstrap (= TRUE) or not (= FALSE).
#'
#' @return A list with the elements
#' \item{p.values}{A matrix of the estimated p-values.}
#' \item{ci}{A matrix of confidence intervals.}
#' @author Justin Esarey
#' @examples
#' \dontrun{
#' 
#' # predict employment levels, cluster on group
#' require(plm)
#' data(EmplUK)
#'
#' emp.1 <- plm(emp ~ wage + log(capital+1), data = EmplUK, 
#'              model = "within", index=c("firm", "year"))
#' cluster.bs.plm(mod=emp.1, dat=EmplUK, cluster="group", ci.level = 0.95, 
#'           boot.reps = 1000, cluster.se = TRUE, report = TRUE, 
#'           prog.bar = TRUE)
#' 
#' # cluster on time
#' 
#' cluster.bs.plm(mod=emp.1, dat=EmplUK, cluster="time", ci.level = 0.95, 
#'             boot.reps = 1000, cluster.se = TRUE, report = TRUE, 
#'             prog.bar = TRUE)
#' 
#' }
#' @rdname cluster.bs.plm
#' @import Formula
#' @import plm
#' @import stats
#' @importFrom utils write.table
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom lmtest coeftest
#' @importFrom sandwich estfun
#' @importFrom sandwich sandwich
#' @references Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. 2008. "Bootstrap-Based Improvements for Inference with Clustered Errors." \emph{The Review of Economics and Statistics} 90(3): 414-427.
#' @export

cluster.bs.plm<-function(mod, dat, cluster="group", ci.level = 0.95, boot.reps = 1000, cluster.se = TRUE, report = TRUE, prog.bar = TRUE){
  
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
  
  form <- mod$formula                                                 # what is the formula of this model?  
  variables <- all.vars(form)                                         # what variables are in this model?
  ind.variables <- names(coefficients(mod))                           # what independent variables are in this model?
  
   if(cluster.se == T){
     
     se.clust <- sqrt(diag(vcovHC(mod, cluster=cluster)))           # retrieve the clustered SEs
     beta.mod <- coefficients(mod)[ind.variables]                   # retrieve the estimated coefficients
     w <- beta.mod / se.clust                                       # calculate the t-test statistic
     
   }else{

     se.beta <- summary(mod)$coefficients[ind.variables,2]          # retrieve the vanilla SEs
     beta.mod <- coefficients(mod)[ind.variables]                   # retrieve the estimated coefficients
     w <- beta.mod / se.beta                                        # calculate the t-test statistic

   }
  
  w.store <- matrix(data=NA, nrow=boot.reps, ncol=length(ind.variables))    # store bootstrapped test statistics
  
  if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
  for(i in 1:boot.reps){
    
    if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
    
    boot.sel <- sample(1:G, size=G, replace=T)                            # randomly select clusters
    
    # pick the observations corresponding to the randomly selected clusters
    boot.ind <- c()                                                       # where the selected obs will be stored
    boot.clust <- c()                                                     # create new cluster index for the bootstrap data
    
    for(k in 1:G){

      obs.sel <- which(as.character(clust.full) 
                   == as.character(unique(clust)[boot.sel[k]] ))          # which observations are in the sampled cluster?
      boot.ind <- c(boot.ind, obs.sel)                                    # append the selected obs index to existing index
      
      boot.clust <- c(boot.clust, rep(k, length(obs.sel)))                # store the new bootstrap cluster index
      
      
    }
    
    boot.dat <- as.data.frame(dat[boot.ind,])                             # create the bootstrapped data
    if(cluster=="group"){
      boot.dat$group.idx <- boot.clust
      boot.dat[[clust.name]] <- boot.clust
      boot.dat$time.idx <- attr(dat, "index")[boot.ind,2]
    }else{
      boot.dat$group.idx <- attr(dat, "index")[boot.ind,1]
      boot.dat$time.idx <- boot.clust
      boot.dat[[clust.name]] <- boot.clust
    }
  
    boot.dat.p <- pdata.frame(boot.dat, index=c("group.idx", "time.idx"))
    boot.call <- mod$call
    boot.call[[3]] <- quote(boot.dat.p)

    # run a model on the bootstrap replicate data
    boot.mod <- suppressWarnings(tryCatch( eval(boot.call), 
                error = function(e){return(NA)}))                                  

    fail <- max(is.na(boot.mod))                                     # determine whether the process created an error
    
    if(fail==0){                                                     # proceed if the model was not in error

      if(cluster.se == T){
        
        se.boot <- tryCatch(sqrt(diag(vcovHC(mod, cluster=cluster))),
                   error = function(e){return(NA)}, 
                   warning = function(w){return(NA)})                              # retrieve the bootstrap clustered SE
        beta.boot <- tryCatch(coefficients(boot.mod)[ind.variables],
                     error = function(e){return(NA)}, 
                     warning = function(w){return(NA)})                            # store the bootstrap beta coefficient
        w.store[i,] <- (beta.boot-beta.mod) / se.boot                              # store the bootstrap test statistic
        
      }else{
        
        se.boot <- tryCatch(summary(boot.mod)$coefficients[ind.variables,2],
                   error = function(e){return(NA)}, 
                   warning = function(w){return(NA)})                               # retrieve the bootstrap vanilla SE
        beta.boot <- tryCatch(coefficients(boot.mod)[ind.variables],
                     error = function(e){return(NA)}, 
                     warning = function(w){return(NA)})                             # retrieve the bootstrap beta coefficient
        w.store[i,] <- (beta.boot-beta.mod) / se.boot                               # calculate the t-test statistic
                
      }
    
    }else{
      w.store[i,] <- NA                                                  # if model didn't converge, store NA as a result 
    }
  
  }
  if(prog.bar==TRUE){close(pb)}
  
  num.fail <- length(attr(na.omit(w.store), "na.action"))         # count the number of times something went wrong
  w.store <- na.omit(w.store)                                     # drop the erroneous bootstrap replicates
  
  
  comp.fun<-function(vec2, vec1){as.numeric(vec1>vec2)}                              # a simple function comparing v1 to v2
  p.store.s <- t(apply(X = abs(w.store), FUN=comp.fun, MARGIN = 1, vec1 = abs(w)))   # compare the BS test stats to orig. result
  p.store <- 1 - ( colSums(p.store.s) / dim(w.store)[1] )                            # calculate the cluster bootstrap p-value

  # compute critical t-statistics for CIs
  crit.t <- apply(X=abs(w.store), MARGIN=2, FUN=quantile, probs=ci.level )
  if(cluster.se == TRUE){
    ci.lo <- beta.mod - crit.t*se.clust
    ci.hi <- beta.mod + crit.t*se.clust
  }else{
    ci.lo <- beta.mod - crit.t*se.beta
    ci.hi <- beta.mod + crit.t*se.beta
  }
  
  print.ci <- cbind(ind.variables, ci.lo, ci.hi)
  print.ci <- rbind(c("variable name", "CI lower", "CI higher"), print.ci)
  
  out.ci <- cbind(ci.lo, ci.hi)
  rownames(out.ci) <- ind.variables
  colnames(out.ci) <- c("CI lower", "CI higher")
  
  out <- matrix(p.store, ncol=1)
  colnames(out) <- c("clustered bootstrap p-value")
  rownames(out) <- ind.variables
  out.p <- cbind(ind.variables, out)
  out.p <- rbind(c("variable name", "cluster bootstrap p-value"), out.p)
  

  printmat <- function(m){
    write.table(format(m, justify="right"), row.names=F, col.names=F, quote=F, sep = "   ")
  }
  
  if(report==T){
    
    if(num.fail!=0){
    cat("\n", "\n", "\n", "****", "Warning: ", num.fail, " out of ", boot.reps, "bootstrap replicate models failed to estimate.", "****", "\n")
    }
    
    cat("\n", "Cluster Bootstrap p-values: ", "\n", "\n")
    printmat(out.p)

    cat("\n", "Confidence Intervals (derived from bootstrapped t-statistics): ", "\n", "\n")
    printmat(print.ci)
    
    
  }
  
  out.list<-list()
  out.list[["p.values"]]<-out
  out.list[["ci"]] <- out.ci
  return(invisible(out.list))
  
}

