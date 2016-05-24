#' Permutation-based variable importance metric for high dimensional datasets appropriate for time to event outcomes, 
#' in the presence of imperfect self-reports or laboratory-based diagnostic tests.
#' 
#' @description Let N and P denote the number of subjects and number of variables in the dataset, respectively. 
#'              Let N** denote the total number of visits, summed over all subjects in the study 
#'              [i.e. N** denotes the number of diagnostic test results available for all subjects in the study].
#'               This algorithm builds a user-defined number of survival trees, using bootstrapped datasets. 
#'               Using the out of bag (OOB) data in each tree, a permutation-based measure of 
#'               variable importance for each of the P variables is obtained.  
#' 
#' @param data name of the data frame that includes the variables subject, testtimes, result
#' @param subject vector of subject IDs of length N**x1.
#' @param testtimes vector of visit or test times of length N**x1.
#' @param result vector of binary diagnostic test results (0 = negative for event of interest; 1 = positive 
#'        for event of interest) of length N**x1.
#' @param sensitivity the sensitivity of the diagnostic test
#' @param specificity the specificity of the diagnostic test
#' @param Xmat a N x P matrix of covariates.
#' @param root.size minimum number of subjects in a terminal node.
#' @param ntree number of survival trees.
#' @param p number of covariate selected at each node to split the tree.
#' @param node For parallel computation, specify the number of nodes.
#' 
#' @return a vector of the ensembled variable importance for modified random survival forest (icRSF).
#' @export
#' @examples
#' library(parallel)
#' data(Xmat)
#' data(pheno)
#' vimp <- icrsf(data=pheno, subject=ID, testtimes=time, result=result, sensitivity=1, 
#'              specificity=1, Xmat=Xmat, root.size=30, ntree=1, p=sqrt(ncol(Xmat)), node=1)
#'
#' @useDynLib icRSF
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim rbinom rexp
#' @import parallel icensmis 
#' 
icrsf <- function(data, subject, testtimes, result, sensitivity, specificity, Xmat, 
                  root.size, ntree, p, node){
  #library(parallel)
  # getchisq function
  getchisq <- function(Dmat, x) {
    J <- dim(Dmat)[2]-1
    parmi <- c(0, log(-log((J:1)/(J+1))))
    parmD <- c(parmi[1], diff(c(0, parmi[-1])))
    parmD0 <- parmD[-1] 
      
    ## log-likelihood function
    if(length(unique(x))==1){
      chisq <- 0
      parm <- rep(0, J+1)
    }else{
      q <- try(optim(parmD, loglikCD, gradlikCD, lower = c(rep(-Inf, 2), rep(1e-8, J-1)), Dmat = Dmat, x = x, method="L-BFGS-B"))
      q0 <- try(optim(parmD0, loglikCD0, gradlikCD0, lower = c(-Inf, rep(1e-8, J-1)), Dmat = Dmat, method="L-BFGS-B"))
      if (class(q)=="try-error" | class(q0)=="try-error") return(rep(0, J+2))
      chisq <- 2*(q0$value - q$value)
      parm <- q$par
    }
    return(c(chisq, parm))
  }
  
  # split function 
  splitpoint <- function(Dmat, x) {
    ux <- sort(unique(x))
    xlist <- (ux[-1] + ux[-length(ux)])/2
    csqs <- sapply(xlist, function(i) getchisq(Dmat, x>=i))
    maxid <- which.max(csqs[1, , drop=F])
    return(c(xlist[maxid], csqs[, maxid]))
  }
  # get parameter function 
  getparm <- function(id, Xmat, treemat, hdidx = 1) {
    stopifnot(is.numeric(id), length(id) == 1)
    vec <- treemat[hdidx, ]
    obs <- Xmat[id, ]
    dir <- (obs[vec[3]] >= vec[4]) + 1
    if (is.na(vec[dir])) {
      fparm <- vec[-(1:5)] ### changes here
      fparm[1] <- fparm[1] + (dir-1)*vec[5]*obs[vec[3]] ### changes here
      return(list(fparm=fparm, id=hdidx)) ### changes here
    } else {
      getparm(id, Xmat, treemat, vec[dir])
    }
  }
  
  # calculate varimp function
  varimp1 <- function(Dm, Xmat, treemat, OOBid){
    #   id <- eval(substitute(subject), data, parent.frame())
    #   time <- eval(substitute(testtimes), data, parent.frame())
    #   result <- eval(substitute(result), data, parent.frame())  
    #   ord <- order(id, time)
    #   if (is.unsorted(ord)) {
    #     id <- id[ord]
    #     time <- time[ord]
    #     result <- result[ord]
    #     data <- data[ord, ]
    #   }  
    #   utime <- sort(unique(time))
    #   
    #   #############################################################################
    #   # Check variables and input parameters check
    #   ############################################################################# 
    #   stopifnot(is.numeric(sensitivity), is.numeric(specificity), is.character(weight))
    #   stopifnot(length(sensitivity) == 1, sensitivity >= 0, sensitivity <= 1, 
    #             length(specificity) == 1, specificity >= 0, specificity <= 1,
    #             length(weight) == 1, all(time > 0), all(is.finite(time)))
    #   if (!weight %in% c("Y", "N")) stop("weight must be 'Y' or 'N'")
    #   if (!all(result %in% c(0, 1))) stop("result must be 0 or 1")
    #   if (any(tapply(time, id, anyDuplicated)))
    #     stop("existing duplicated visit times for some subjects")
    #   #############################################################################
    #   # Compute D matrix
    #   #############################################################################  
    #   
    #   timen0 <- (time != 0)
    #   Dm <- icensmis:::dmat(id[timen0], time[timen0], result[timen0], sensitivity,
    #                         specificity, 1)
    #   
    nOOB <- length(OOBid)
    ncov <- ncol(Xmat)
    
    newpred <- function(inc){
      m <- rep(NA, inc)
      m[1] <- 1
      #    m[1, -1] <- parm
      return(list(mat=m, nxt=2, inc=inc))
    }
    insertpred <- function(pred, newval){
      newidx <- pred$nxt
      if (pred$nxt == length(pred) + 1) {
        pred$mat <- c(pred$mat, rep(NA, pred$inc))
      }
      pred$mat[newidx] <- newval
      #  pred$mat[newidx, -1] <- parm
      pred$nxt <- pred$nxt + 1
      return(pred)
    }
    
    getpred <- function(x, treemat, hdidx = 1, pmat) {
      vec <- treemat[hdidx, ]
      dir <- (x[vec[3]] >= vec[4]) + 1
      if(is.na(hdidx) == FALSE){
        if(hdidx == 1){
          pmat <- newpred(3)
        } else
        {
          pmat <- insertpred(pmat, hdidx)
        }
        hdidx <- vec[dir]
        pmat <- getpred(x, treemat, hdidx, pmat)
      }
      return(pmat)
    }
    #############################################################################
    # Retreive terminal parameters and compute log likelihood
    #############################################################################  
    OOBloglik <- function(Dmat, Xmat, OOBid, covx, treemat) {
      OOBDm <- Dmat[OOBid, ]
      OOBxm <- Xmat[OOBid, ]
      if(covx > 0) {
        OOBxm[, covx] <- sample(OOBxm[, covx], nOOB, replace=FALSE)
      }
      p <- lapply(1:nOOB, function(x) getparm(x, OOBxm, treemat, 1))
      parms <- t(sapply(p, function(x) cumsum(x$fparm)))
      l <- rowSums(OOBDm*cbind(1, exp(-exp(parms))))
      freq <- unlist(sapply(1:nOOB, function(x) getpred(OOBxm[x, ], treemat, 1, 1)[[1]]))
      # get rid of the NA's
      freq <- freq[!is.na(freq)]
      # this is to calculate the number of OOB data pass through each row of the tree
      f <- sapply(1:dim(treemat)[1], function(x) sum(freq==x)/nOOB)
      # this is to get the id of the parent node of the terminal node    
      tid <- sapply(p, function(x) x$id)
      loglik <- ifelse(l>0, log(l), NA)
      
      return(loglik)
    }
    cov <- unique(treemat[, 3])
    cov <- cov[!is.na(cov)]
    varimp <- rep(0, ncov)
    lik.org <- OOBloglik(Dm, Xmat, OOBid, 0, treemat)
    
    ## This is to get the loglik after permuation of each of the selected covariates in the tree
    permutelik <- lapply(cov, function(x) OOBloglik(Dm, Xmat, OOBid, x, treemat))  
    varimp[cov] <- sapply(permutelik, function(x) sum(lik.org - x, na.rm=T))
    return(varimp)
  }

  
  #############################################################################
  # Function to build trees
  ############################################################################# 
  
  treebuilder <- function(Dm, Xmat, root.size, p){
    
    getchisq <- function(Dmat, x) {
      J <- dim(Dmat)[2]-1
      parmi <- c(0, log(-log((J:1)/(J+1))))
      parmD <- c(parmi[1], diff(c(0, parmi[-1])))
      parmD0 <- parmD[-1] 
      
      ## log-likelihood function
      if(length(unique(x))==1){
        chisq <- 0
        parm <- rep(0, J+1)
      }else{
        q <- try(optim(parmD, loglikCD, gradlikCD, lower = c(rep(-Inf, 2), rep(1e-8, J-1)), Dmat = Dmat, x = x, method="L-BFGS-B"))
        q0 <- try(optim(parmD0, loglikCD0, gradlikCD0, lower = c(-Inf, rep(1e-8, J-1)), Dmat = Dmat, method="L-BFGS-B"))
        if (class(q)=="try-error" | class(q0)=="try-error") return(rep(0, J+2))
        chisq <- 2*(q0$value - q$value)
        parm <- q$par
      }
      return(c(chisq, parm))
    }
    
    
    #   id <- eval(substitute(subject), data, parent.frame())
    #   time <- eval(substitute(testtimes), data, parent.frame())
    #   result <- eval(substitute(result), data, parent.frame())  
    #   ord <- order(id, time)
    #   if (is.unsorted(ord)) {
    #     id <- id[ord]
    #     time <- time[ord]
    #     result <- result[ord]
    #     data <- data[ord, ]
    #   }  
    #   utime <- sort(unique(time))
    #   
    #   #############################################################################
    #   # Check variables and input parameters check
    #   ############################################################################# 
    #   stopifnot(is.numeric(sensitivity), is.numeric(specificity),
    #             is.numeric(root.size))
    #   stopifnot(length(sensitivity) == 1, sensitivity >= 0, sensitivity <= 1, 
    #             length(specificity) == 1, specificity >= 0, specificity <= 1,
    #             length(root.size) == 1, root.size >= 1, root.size <= nrow(Xmat),
    #             all(time > 0), all(is.finite(time)))
    #   if (!all(result %in% c(0, 1))) stop("result must be 0 or 1")
    #   if (any(tapply(time, id, anyDuplicated)))
    #     stop("existing duplicated visit times for some subjects")
    #############################################################################
    # Compute D matrix
    #############################################################################  
    #   
    #   timen0 <- (time != 0)
    #   Dm <- icensmis:::dmat(id[timen0], time[timen0], result[timen0], sensitivity,
    #                         specificity, 1)
    J <- ncol(Dm) - 1
    nsub <- nrow(Dm)
    
    parmi <- c(0, log(-log((J:1)/(J+1))))
    parmD <- c(parmi[1], diff(c(0, parmi[-1])))
    parmD0 <- parmD[-1] 
    
    #  simdata <- simoutcome$data
    #  id <- simoutcome$idx
    nsub <- nrow(Dm)
    ncov <- ncol(Xmat)
    row.names(Dm) <- 1:nsub
    J <- ncol(Dm) - 1
    parmi <- c(0, log(-log((J:1)/(J+1))))
    parmD <- c(parmi[1], diff(c(0, parmi[-1])))
    parmD0 <- parmD[-1] 
    ns <- p
    
    bootid <- sample(1:nsub, nsub, replace = T)
    #bootdm <- Dmat[bootid, ]
    #bootxm <- Xmat[bootid, ]
    
    OOBid <- setdiff(1:nsub, unique(bootid))
    #  OOB <- simdata[setdiff(1:nsub, unique(bootid)), ]
    #  nOOB <- nrow(OOB)
    ################# TREE BUILDING ####################
    
    ## function to create a new tree
    newtree <- function(firstval, firstc, inc, parm) {
      m <- matrix(rep(NA, inc*(J+5)), nrow = inc, ncol = J + 5)
      m[1, 3] <- firstval
      m[1, 4] <- firstc
      m[1, -(1:4)] <- parm
      return(list(mat=m, nxt=2, inc=inc))
    }
    
    ## function to insert value to the tree
    insert <- function(hdidx, dir, tr, newval, newc, parm) {
      newidx <- tr$nxt
      if (tr$nxt == nrow(tr$mat) + 1) {
        tr$mat <- rbind(tr$mat, matrix(rep(NA, tr$inc*(J+5)), nrow=tr$inc, ncol=J+5))
      }
      tr$mat[newidx, 3] <- newval
      tr$mat[newidx, 4] <- newc
      tr$mat[newidx, -(1:4)] <- parm
      tr$mat[hdidx, dir] <- newidx
      tr$nxt <- tr$nxt + 1
      return(tr)
    }
    
    training <- function(bootid, Dmat, Xmat, hdidx, dir, tree, root.size) {
      ## stopping rule
      if (length(unique(row.names(Dmat[bootid, ])))>=root.size) {
        #print(paste("unique observations:", length(unique(row.names(bootdata)))))
        #Dmat <- bootdata[, 1:(J+1)]
        #print(paste("Dmat:", dim(Dmat)))
        #Xmat <- bootdata[, -(1:(J+1))]
        #       set.seed(2)
        selid <- sample(1:ncov, ns, replace=FALSE)
        #       print(paste("firstselect:", selid))
        bootdm <- Dmat[bootid, ]
        bootxm <- Xmat[bootid, ]
        ind <- apply(bootxm[, selid], 2, function(x) length(unique(x))>1)
        selid <- selid[ind]
        X <- bootxm[, selid]
        bsplit <- bestsplitC(bootdm, X, getchisq)
        #print(paste("number of candidate covariates:", length(selid)))
        if(length(selid) > 0){
          bvar <- selid[bsplit[[1]]]
          #print(paste("bestsplit variable:", bvar))
          bcp <- bsplit[[2]][1]
          bootid.L <- bootid[X[, bsplit[[1]]] <= bcp]
          bootid.R <- bootid[X[, bsplit[[1]]] > bcp]
          parm <- bsplit[[2]][-(1:2)]
          if (length(bootid) == nsub) {
            tree <- newtree(bvar, bcp, 3, parm)
          } else {
            tree <- insert(hdidx, dir, tree, bvar, bcp, parm)
          }
          a <- tree$nxt - 1
          tree <- training(bootid.L,Dmat, Xmat, a, 1, tree, root.size)
          tree <- training(bootid.R, Dmat, Xmat, a, 2, tree, root.size)
        } 
      }
      return(tree)
    }
    tr <- training(bootid, Dm, Xmat, 1, 1, 1, root.size)
    return(list(tree=tr[[1]], OOBid=OOBid))    
  }
  
  ## Main body
  id <- eval(substitute(subject), data, parent.frame())
  time <- eval(substitute(testtimes), data, parent.frame())
  result <- eval(substitute(result), data, parent.frame())  
  ord <- order(id, time)
  if (is.unsorted(ord)) {
    id <- id[ord]
    time <- time[ord]
    result <- result[ord]
    data <- data[ord, ]
  }  
  utime <- sort(unique(time))
  
  #############################################################################
  # Check variables and input parameters check
  ############################################################################# 
  stopifnot(is.numeric(sensitivity), is.numeric(specificity),
            is.numeric(root.size),  is.numeric(node), 
            is.numeric(node), is.numeric(p))
  stopifnot(length(sensitivity) == 1, sensitivity >= 0, sensitivity <= 1, 
            length(specificity) == 1, specificity >= 0, specificity <= 1,
            length(root.size) == 1, root.size >= 1, root.size <= nrow(Xmat),
            length(node) ==1, node >=1,
            length(p) ==1, p >=1, p<=ncol(Xmat),
            all(time > 0), all(is.finite(time)))
  if (!all(result %in% c(0, 1))) stop("result must be 0 or 1")

  if (any(tapply(time, id, anyDuplicated)))
    stop("existing duplicated visit times for some subjects")
  #############################################################################
  # Compute D matrix
  #############################################################################  
  
  timen0 <- (time != 0)
  Dm <- dmat(id[timen0], time[timen0], result[timen0], sensitivity,
                        specificity, 1)
  #sidx <- dat$idx
  trees <- mclapply(1:ntree, function(x) treebuilder(Dm, Xmat, root.size, p),
                    mc.cores=node)
  treemats <- mclapply(trees, function(x) x[[1]], mc.cores=node)  
  OOBids <- mclapply(trees, function(x) x[[2]], mc.cores=node)

    vimp <- mclapply(1:ntree, function(x) varimp1(Dm, Xmat, treemats[[x]], OOBid=OOBids[[x]]))

  vimp.adj <- mclapply(vimp, function(x) ifelse(x>0, x, 0))
  vimp1 <- do.call("rbind", vimp.adj)
  vimp2 <- colSums(vimp1, na.rm=T)
  return(ensemble.vimp=vimp2)
  }
  
