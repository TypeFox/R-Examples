#' Cluster-Adjusted Confidence Intervals And p-Values For mlogit
#'
#' Computes p-values and confidence intervals for multinomial logit models based on cluster-specific model estimation (Ibragimov and Muller 2010). A separate model is estimated in each cluster, and then p-values and confidence intervals are computed based on a t/normal distribution of the cluster-specific estimates.
#'
#' @param mod A model estimated using \code{mlogit}.
#' @param dat The data set used to estimate \code{mod}.
#' @param cluster A formula of the clustering variable.
#' @param ci.level What confidence level should CIs reflect?
#' @param report Should a table of results be printed to the console?
#' @param truncate Should outlying cluster-specific beta estimates be excluded?
#'
#' @return A list with the elements
#' \item{p.values}{A matrix of the estimated p-values.}
#' \item{ci}{A matrix of confidence intervals.}
#' @author Justin Esarey
#' @note Confidence intervals are centered on the cluster averaged estimate, which can diverge from original model estimates if clusters have different numbers of observations. Consequently, confidence intervals may not be centered on original model estimates. Any cluster for which all coefficients cannot be estimated will be automatically dropped from the analysis. If truncate = TRUE, any cluster for which any coefficient is more than 6 times the interquartile range from the cross-cluster mean will also be dropped as an outlier.
#' @examples
#' \dontrun{
#' 
#' # example: predict type of heating system installed in house
#' require(mlogit)
#' data("Heating", package = "mlogit")
#' H <- Heating
#' H.ml <- mlogit.data(H, shape="wide", choice="depvar", varying=c(3:12))
#' m <- mlogit(depvar~ic+oc, H.ml)
#' 
#' # compute cluster-adjusted p-values
#' cluster.im.h <- cluster.im.mlogit(m, H.ml, ~ region)
#' 
#' }
#' @rdname cluster.im.mlogit
#' @importFrom mlogit mlogit mlogit.data hmftest mFormula is.mFormula mlogit.optim cov.mlogit cor.mlogit rpar scoretest med rg stdev qrpar prpar drpar
#' @references Ibragimov, Rustam, and Ulrich K. Muller. 2010. "t-Statistic Based Correlation and Heterogeneity Robust Inference." \emph{Journal of Business & Economic Statistics} 28(4): 453-468. 
#' @import stats
#' @importFrom utils write.table
#' @export

cluster.im.mlogit<-function(mod, dat, cluster, ci.level = 0.95, report = TRUE, truncate = FALSE){
  
  form <- mod$formula                                                    # what is the formula of this model?  
  variables <- all.vars(form)                                            # what variables are in this model?
  used.idx <- which(rownames(dat) %in% rownames(mod$mod))                # what observations are used?
  dat <- dat[used.idx,]                                                  # keep only active observations
  ind.variables <- names(coefficients(mod))                              # what independent variables are in this model?
  "%w/o%" <- function(x, y) x[!x %in% y]                                 # create a without function (see ?match)
  dv <- variables %w/o% all.vars(update(form, 1 ~ .))                    # what is the dependent variable?
  
  # obtain the clustering variable
  clust.name <- all.vars(cluster)                                        # name of the cluster variable
  dat.rs <- subset(dat, select = clust.name )                            # select cluster variable from data set
  dat.rs$id.zz <- attr(dat, 'index')[used.idx,1]                                 # choice index
  dat.rs$ti.zz <- attr(dat, 'index')[used.idx,2]                                 # alternative index
  clust <- reshape(dat.rs, timevar="ti.zz",                              # reshape long to wide, store as clust
          idvar=c("id.zz", clust.name), direction="wide")[[clust.name]]
  if(sum(is.na(clust)>0)){stop("missing cluster indices")}               # check for missing cluster indices
  G<-length(unique(clust))                                               # how many clusters are in this model?
  
  
  b.clust<-matrix(data=NA, nrow=G, ncol=length(coefficients(mod)))            # store bootstrapped beta values
  
  dropped.nc <- 0                                                             # store number of non-converged clusters
  for(i in 1:G){
    
    clust.ind <- which(dat[[clust.name]] == unlist(unique(clust))[i])         # select obs in cluster i
    
    clust.dat <- dat[clust.ind,]                                              # create the cluster i data set
    clust.mod.call <- mod$call                                                # get original model call
    clust.mod.call[[3]] <- quote(clust.dat)                                   # modify call for cluster data set
    
    clust.mod <- suppressWarnings(tryCatch(eval(clust.mod.call),              # estimate model on cluster dataset
                                  error = function(e){return(NULL)}))    
    
    fail <- is.null(clust.mod)                                                # detect model error
    

    if(fail==F){                                                              # if model ran...

      b.clust[i,] <- coefficients(clust.mod)                                  # store the cluster i beta coefficient

    }
  }
  
  # purge non-converged clusters
  b.clust <- na.omit(b.clust)
  dropped.nc <- length(attr(b.clust, "na.action"))                            # record number of models dropped
  
  G.t <- dim(b.clust)[1]
  if(G.t == 0){stop("all clusters were dropped (see help file).")}
  
  # remove clusters with outlying betas
  dropped <- 0                                                                 # store number of outlying estimates
  if(truncate==TRUE){                                                          # if drop outlying betas...

    IQR <- apply(FUN=quantile, MARGIN=2, X=b.clust, probs=c(0.25, 0.75))       # calculate inter-quartile range
    
    b.clust.save <- b.clust                                                    # non-outlying beta identifier
    for(i in 1:dim(b.clust)[2]){
      b.clust.save[,i] <- ifelse( abs(b.clust[,i]) >                           # determine which estimates to save
            (abs(mean(b.clust[,i])) + 6*abs(IQR[2,i] - IQR[1,i])), 0, 1)
    }
    
    save.clust <- apply(X=b.clust.save, MARGIN=1, FUN=min)                     # which rows had outlying betas
    dropped <- dim(b.clust)[1] - sum(save.clust)                               # how many clusters dropped?
    
    b.clust.adj <- cbind(b.clust, save.clust)                                  # append the outlier ID to beta matrix
    
    b.clust <- subset(b.clust,                                                 # drop the outlying betas
              subset=(save.clust==1), select=1:dim(b.clust)[2])
    
  }
  
  G <- dim(b.clust)[1]
  if(G == 0){stop("all clusters were dropped (see help file).")}
  
  b.hat <- colMeans(b.clust)                                # calculate the avg beta across clusters
  b.dev <- sweep(b.clust, MARGIN = 2, STATS = b.hat)        # sweep out the avg betas
  #s.hat <- sqrt( (1 / (G-1)) * colSums(b.dev^2) )          # manually calculate the SE of beta estimates across clusters (deprecated)
  vcv.hat <- cov(b.dev)                                     # calculate VCV matrix
  s.hat <- sqrt(diag(vcv.hat))                              # calculate standard error
  
  t.hat <- sqrt(G) * (b.hat / s.hat)                        # calculate t-statistic
  
  # compute p-val based on # of clusters
  p.out <- 2*pmin( pt(t.hat, df = G-1, lower.tail = TRUE), pt(t.hat, df = G-1, lower.tail = FALSE) )
  
  # compute CIs
  ci.lo <- b.hat - qt((1-ci.level)/2, df=(G-1), lower.tail=FALSE)*(s.hat/sqrt(G))
  ci.hi <- b.hat + qt((1-ci.level)/2, df=(G-1), lower.tail=FALSE)*(s.hat/sqrt(G))
    
  out <- matrix(p.out, ncol=1)
  out.p <- cbind( names(coefficients(summary(mod))), round(out,3) )
  
  out.ci <- cbind(ci.lo, ci.hi)
  rownames(out.ci) <- names(coefficients(summary(mod)))
  colnames(out.ci) <- c("CI lower", "CI higher")
  
  print.ci <- cbind(names(coefficients(summary(mod))), ci.lo, ci.hi)
  print.ci <- rbind(c("variable name", "CI lower", "CI higher"), print.ci)
  

  printmat <- function(m){
    write.table(format(m, justify="right"), row.names=F, col.names=F, quote=F, sep = "   ")
  }
    
  if(report==T){
    
    cat("\n", "Cluster-Adjusted p-values: ", "\n", "\n")
    printmat(out.p)
    
    cat("\n", "Confidence Intervals (centered on cluster-averaged results):", "\n", "\n")
    printmat(print.ci)

    if(dropped.nc > 0){cat("\n", "Note:", dropped.nc, "clusters were dropped due to missing coefficients (see help file).", "\n", "\n")}
    if(dropped > 0){cat("\n", "Note:", dropped, "clusters were dropped as outliers (see help file).", "\n", "\n")}
    
  }
  

  out.list<-list()
  out.list[["p.values"]]<-out
  out.list[["ci"]]<-out.ci
  return(invisible(out.list))
  
}
