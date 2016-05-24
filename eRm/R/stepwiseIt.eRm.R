#function for stepwise item elimination

stepwiseIt.eRm <- function(object, criterion = list("itemfit"), alpha = 0.05, verbose = TRUE,
                       maxstep = NA)
{
# object of class dRm
# criterion: either list("itemfit") or list("LRtest", splitcr) od list("Waldtest", splitcr)

  #-------- sanity checks ---------
  dummy <- match.arg(criterion[[1]], c("itemfit","LRtest","Waldtest"))
  if (!is.list(criterion)) stop("Criterion must be provided as list!")
  if (!any(class(object) == "dRm")) stop("Stepwise elimination implemented for dichotomous Rasch models only!")
  #------- end sanity checks ------

  X.new <- object$X
  K <- dim(X.new)[2]
  if (is.na(maxstep)) maxstep <- K

  if (length(criterion) == 2) {
      splitcr <- criterion[[2]]
  } else {
      splitcr <- "median"
  }


  #---------------- start elimination  ----------------
  i <- 0
  it.el <- rep(NA, K)                              #initialize outputs
  el.names <- rep(NA, K)
  wald.mat <- matrix(NA, ncol = 2, nrow = K)
  itemfit.mat <- matrix(NA, ncol = 3, nrow = K)
  LR.mat <- matrix(NA, ncol = 3, nrow = K)
    
  repeat
  {
    if((dim(X.new)[2]) == 2) {
      warning("Only 2 items left: No Rasch homogeneous itemset found!", call. = FALSE)
      break
    }
    if (i == maxstep) {
      warning("Maximum number of steps reached!", call. = FALSE)
      break
    }

    i <- i + 1
    res <- RM(X.new)                                     #fit Rasch

    #---------------- itemfit criterion ------------
    if (criterion[[1]] == "itemfit") {
      pres <- person.parameter(res)                        #person parameters
      it.res <- itemfit(pres)                              #compute itemfit
      pvalvec <- 1-pchisq(it.res$i.fit, it.res$i.df)       #vector with pvalues
      pvalsig <- which(pvalvec < alpha)                    #significant p-values

      if (length(pvalsig) > 0) {
        it.el[i] <- which(it.res$i.fit == max(it.res$i.fit))[1]
        ie <- it.el[i]
        itemfit.mat[i,] <- c(it.res$i.fit[ie], it.res$i.df[ie], pvalvec[ie])
        if (verbose) cat("Eliminated item - Step ",i,": ",colnames(X.new)[it.el[i]],"\n", sep = "")
        el.names[i] <- colnames(X.new)[it.el[i]]
        X.new <- X.new[,-it.el[i]]
      } else break
    }
    #-------------- end itemfit criterion -----------
      
    #------------------ Waldtest criterion ----------
    if (criterion[[1]] == "Waldtest")
    {
      wald.res <- Waldtest(res, splitcr = splitcr)        #compute Waldtest
      zvalvec <- abs(wald.res$coef.table[,1])             #absolute z-values
      pvalvec <- wald.res$coef.table[,2]                  #vector with pvalues
      pvalsig <- which(pvalvec < alpha)                   #significant p-values 
      if (length(pvalsig) > 0) {
        elpos <- which(zvalvec == max(zvalvec))[1]       #exclude maximum z-value Waldtest
        wald.mat[i,] <- wald.res$coef.table[elpos,]
        if (length(wald.res$it.ex) > 0) elpos <- elpos + sum(wald.res$it.ex <= elpos)  #if items couldn't computed in Waldtest
        it.el[i] <- elpos
        el.names[i] <- colnames(X.new)[it.el[i]]
        if (verbose) cat("Eliminated item - Step ",i,": ",el.names[i],"\n", sep = "")
        X.new <- X.new[,-it.el[i]]
      } else break
    }
      
    #-------------- LRtest criterion ----------------
    if (criterion[[1]] == "LRtest")                           #uses Waldtest but stops when LRtest is sig.
    {
      lr.res <- LRtest(res, splitcr = splitcr)
      if(lr.res$pvalue < alpha) {
        wald.res <- Waldtest(res, splitcr = splitcr)        #compute Waldtest
        zvalvec <- abs(wald.res$coef.table[,1])             #absolute z-values
        elpos <- which(zvalvec == max(zvalvec))[1]       #exclude maximum z-value Waldtest
        if (length(wald.res$it.ex) > 0) elpos <- elpos + sum(wald.res$it.ex <= elpos)  #if items couldn't computed in Waldtest
        it.el[i] <- elpos
        LR.mat[i,] <- c(lr.res$LR, lr.res$df, lr.res$pvalue)
        el.names[i] <- colnames(X.new)[it.el[i]]
        if (verbose) cat("Eliminated item - Step ",i,": ",el.names[i],"\n", sep = "")
        X.new <- X.new[,-it.el[i]]
      } else break
    }
    #----------- end LRtest criterion ---------
  }
 #--------------------- end stepwise------------------
 
  #labeling
  el.names <- el.names[!is.na(el.names)]
  if (all(is.na(el.names))) {
    warning("No items eliminated! Each of them fits the Rasch model!", call. = FALSE)
    itemfit.mat <- NULL
    LR.mat <- NULL
    wald.mat <- NULL
    criterion[[1]] <- "none"
  }

  if (criterion[[1]] == "itemfit")
  {
   itemfit.mat <- rbind(itemfit.mat[!is.na(rowSums(itemfit.mat)),])
   rownames(itemfit.mat) <- paste("Step ",1:length(el.names),": ",el.names,sep = "")
   colnames(itemfit.mat) <- c("Chisq", "df","p-value")
  } else {
    itemfit.mat <- NULL
  }
  if (criterion[[1]] == "Waldtest")
  {
    wald.mat <- rbind(wald.mat[!is.na(rowSums(wald.mat)),])
    rownames(wald.mat) <- paste("Step ",1:length(el.names),": ",el.names,sep = "")
    colnames(wald.mat) <- c("z-statistic", "p-value")
  } else {
    wald.mat <- NULL
  }
  if (criterion[[1]] == "LRtest")
  {
    if (i == maxstep) {
      LR.mat <- rbind(LR.mat[!is.na(rowSums(LR.mat)),])
      rownames(LR.mat) <- paste("Step ",1:length(el.names),": ",el.names,sep = "")
    } else {
      LR.mat <- rbind(LR.mat[!is.na(rowSums(LR.mat)),], c(lr.res$LR, lr.res$df, lr.res$pvalue))
      rownames(LR.mat) <- c(paste("Step ",1:length(el.names),": ",el.names,sep = ""), paste("Step ", i,": None", sep = ""))
    }
    colnames(LR.mat) <- c("LR-value", "Chisq df", "p-value")
  } else {
    LR.mat <- NULL
  }

  result <- list(X = X.new, fit = res, it.elim = el.names, res.wald = wald.mat, res.itemfit = itemfit.mat,
                res.LR = LR.mat, nsteps = i-1)
  class(result) <- "step"
  result
}
