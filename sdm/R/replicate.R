# Author: Babak Naimi, naimi.b@gmail.com
# Date :  April 2015
# Version 2.0
# Licence GPL v3


.replicateMethods=new('.ReplicateMethods')

.replicateMethods$Methods['cross_validation'] <- new('.methodTemplate',
                                                     name='cross_validation',
                                                     aliases=c('cv','cross-validation'),
                                                     arguments=c(x='numeric',replicates='numeric',nfolds='numeric',stratify='logical',family='character'),
                                                     Help='K-fold cross validation partitioning...\n The procedure can be repeated r times (i.e., r separate k-fold cross-validation; r=replicates)',
                                                     Function=function(x,...) {
                                                       o <- list()
                                                       r <- list(...)[['replicates']]
                                                       k <- list(...)[['nfolds']]
                                                       family <- list(...)[['family']]
                                                       st <- list(...)[['stratify']]
                                                       if (is.null(r)) r <- 1
                                                       if (is.null(st)) st <- TRUE
                                                       if (is.null(family)) {
                                                         if (all(sort(unique(x)) == c(0,1))) family <- 'binomial'
                                                       }
                                                       
                                                       if (family %in% c('bernouli','binomial')) {
                                                         
                                                         np <- length(x[x == 1])
                                                         na <- length(x) - np
                                                         min.n <- min(c(np,na))
                                                         
                                                         if (min.n == 0) stop('The distribution is not binomial!')
                                                         if (is.null(k)) {
                                                           if (min.n > 5) k <- 5
                                                           else k <- min.n
                                                         } else if (k > min.n) {
                                                           k <- min.n
                                                           warning(paste('number of folds cannot be greater than the number of either presence or absence records. nfolds is changed to:',min.n))
                                                         }
                                                         if (!st) {
                                                           o[[1]] <- rep(seq_len(k), length.out=length(x))
                                                           o[[2]] <- replicate(r, sample(1:length(x)))
                                                           
                                                           if (!all(unlist(lapply(unique(o[[1]]),function(i) {
                                                             unlist(lapply(1:r,function(j) length(unique(x[o[[2]][o[[1]] == i,j]])) == 2))
                                                           })))) {
                                                             st <- TRUE
                                                             warning('To avoid error in partitioning, stratify is turned to TRUE')
                                                           }
                                                         }
                                                         if (st) {
                                                           o[[1]] <- c(rep(seq_len(k), length.out=np),rep(seq_len(k), length.out=na))
                                                           o[[2]] <- rbind(replicate(r, sample(which(x == 1))),replicate(r, sample(which(x == 0))))
                                                         }
                                                       } else {
                                                         if (k > length(x)) {
                                                           k <- length(x)
                                                           warning('nfolds cannot be greater than the number of observations! It is changed to the number of observations (became "leave one out"/Jackknife procedure)')
                                                         }
                                                         o[[1]] <- rep(seq_len(k), length.out=length(x))
                                                         o[[2]] <- replicate(r, sample(1:length(x)))
                                                       }
                                                       oo <- list()
                                                       for (k in unique(o[[1]])) {
                                                         w <- which(o[[1]] != k)
                                                         for (i in 1:ncol(o[[2]])) {
                                                           oo[[length(oo)+1]] <- list(method='cross_validation',train=o[[2]][w,i],test=o[[2]][-w,i])
                                                         }
                                                       }
                                                       oo
                                                     })


.replicateMethods$Methods['bootstrap'] <- new('.methodTemplate',
                                              name='bootstrap',
                                              aliases=c('boot','bootstrapping','bootstraping'),
                                              arguments=c(x='numeric',replicates='numeric'),
                                              Help='Partition data into train and test using bootstrapping.',
                                              Function=function(x,...) {
                                                o <- list()
                                                r <- list(...)[['replicates']]
                                                family <- list(...)[['family']]
                                                if (is.null(r)) r <- 1
                                                if (is.null(family)) {
                                                  if (all(sort(unique(x)) == c(0,1))) family <- 'binomial'
                                                }
                                                .sampleboot <- function(x,family) {
                                                  .s <- function(x) {
                                                    o <- list()
                                                    n <- length(x)
                                                    s <- sample(c(1:n),n,replace=TRUE)
                                                    w <- c(1:n) %in% s
                                                    if (length(which(w)) == n) {
                                                      # if all records are selected for train (i.e., there is no records remained for test), it is repeated to avoid error (up to 100 times)
                                                      for (i in 1:100) {
                                                        s <- sample(1:n,n,replace=TRUE)
                                                        w <- c(1:n) %in% s
                                                        if (length(which(w)) < n) break
                                                      }
                                                    }
                                                    o[[1]] <- ifelse(w,1,2)
                                                    o[[2]] <- s
                                                    o
                                                  }
                                                  o <- .s(x)
                                                  if (family %in% c('bernouli','binomial')) {
                                                    w1 <- length(unique(x[o[[2]][o[[1]] == 1]])) == 2
                                                    w2 <- length(unique(x[o[[2]][o[[1]] == 2]])) == 2
                                                    if (!c(w1 & w2)) {
                                                      for (i in 1:100) {
                                                        o <- .s(x)
                                                        w1 <- length(unique(x[o[[2]][o[[1]] == 1]])) == 2
                                                        w2 <- length(unique(x[o[[2]][o[[1]] == 2]])) == 2
                                                        if (c(w1 & w2)) break
                                                      }
                                                    }
                                                  }
                                                  o
                                                }
                                                w <- replicate(r,.sampleboot(x,family=family),simplify=FALSE)
                                                o[[1]] <- list()
                                                o[[2]] <- matrix(ncol=r,nrow=length(x))
                                                for (i in 1:r) {
                                                  o[[1]][[i]] <- w[[i]][[1]]
                                                  o[[2]][,i] <- w[[i]][[2]]
                                                }
                                                oo <- list()
                                                for (i in 1:ncol(o[[2]])) {
                                                  folds <- o[[1]][[i]]
                                                  w <- which(folds == 1)
                                                  oo[[i]] <- list(method='bootstrap',train=o[[2]][w,i],test=o[[2]][-w,i])
                                                }
                                                oo
                                              })
#-----
.replicateMethods$Methods['subsampling'] <- new('.methodTemplate',
                                                name='subsampling',
                                                aliases=c('sub','subsample'),
                                                arguments=c(x='numeric',replicates='numeric',test.percent='numeric',stratify='logical',family='character'),
                                                Help='Subsampling procedure randomly draw (without replacement) a proportion of data as test dataset, and keep the remained records for training',
                                                Function=function(x,...) {
                                                  o <- list()
                                                  r <- list(...)[['replicates']]
                                                  testp <- list(...)[['test.percent']]
                                                  st <- list(...)[['stratify']]
                                                  family <- list(...)[['family']]
                                                  if (is.null(r)) r <- 1
                                                  if (is.null(st)) st <- TRUE
                                                  if (is.null(family)) {
                                                    if (all(sort(unique(x)) == c(0,1))) family <- 'binomial'
                                                  }
                                                  if (is.null(testp)) testp <- 0.3
                                                  else if (testp > 100) {
                                                    warning('test percentage should be less than 100, it is changed to the default, 30!')
                                                    testp <- 0.3
                                                  }
                                                  
                                                  if (testp > 1) testp <- testp / 100
                                                  if (family %in% c('bernouli','binomial')) {
                                                    np <- length(x[x == 1])
                                                    na <- length(x) - np
                                                    npt <- ceiling(np * testp)
                                                    nat <- ceiling(na * testp)
                                                    if (np == npt) {
                                                      if (np > 1) {
                                                        npt <- npt - 1
                                                        st <- TRUE
                                                      } else stop('Subsampling is not possible with only 1 presence record...!')
                                                    }
                                                    if (na == nat) {
                                                      if (na > 1) {
                                                        nat <- nat - 1
                                                        st <- TRUE
                                                      } else stop('Subsampling is not possible with only 1 absence record...!')
                                                    }
                                                    if (!st) {
                                                      test.nr <- ceiling(length(x) * testp)
                                                      o[[1]] <- rep(1:2,times=c((length(x)-test.nr),test.nr))
                                                      o[[2]] <- replicate(r, sample(1:length(x)))
                                                      u <- unlist(lapply(unique(o[[1]]),function(i) {
                                                        unlist(lapply(1:r,function(j) length(unique(x[o[[2]][o[[1]] == i,j]])) == 2))
                                                      }))
                                                      if (!all(u)) {
                                                        st <- TRUE
                                                        warning('To avoid error in partitioning, stratified partitioning is used!')
                                                      }
                                                    }
                                                    
                                                    if (st) {
                                                      o[[1]] <- c(rep(1:2,times=c(np-npt,npt)),rep(1:2,times=c(na-nat,nat)))
                                                      o[[2]] <- rbind(replicate(r, sample(which(x == 1))),replicate(r, sample(which(x == 0))))
                                                    } 
                                                  } else {
                                                    test.nr <- ceiling(length(x) * testp)
                                                    o[[1]] <- rep(1:2,times=c((length(x)-test.nr),test.nr))
                                                    o[[2]] <- replicate(r, sample(1:length(x)))
                                                  }
                                                  oo <- list()
                                                  w <- which(o[[1]] == 1)
                                                  for (i in 1:ncol(o[[2]])) {
                                                    oo[[i]] <- list(method='subsampling',train=o[[2]][w,i],test=o[[2]][-w,i])
                                                  }
                                                  oo
                                                })





# .replicateMethods <- new('.ReplicateMethods')
# .replicateMethods$Methods['cross_validation'] <- new('.methodTemplate',
#                                 name='cross_validation',
#                                 aliases=c('cv','cross-validation'),
#                                 arguments=c(x='numeric',replicates='numeric',nfolds='numeric',stratify='logical',family='character'),
#                                 Help='K-fold cross validation partitioning...\n The procedure can be repeated r times (i.e., r separate k-fold cross-validation; r=replicates)',
#                                 Function=function(x,...) {
#                                   o <- list()
#                                   r <- list(...)[['replicates']]
#                                   k <- list(...)[['nfolds']]
#                                   family <- list(...)[['family']]
#                                   st <- list(...)[['stratify']]
#                                   if (is.null(r)) r <- 1
#                                   if (is.null(st)) st <- TRUE
#                                   if (is.null(family)) {
#                                     if (all(sort(unique(x)) == c(0,1))) family <- 'binomial'
#                                   }
#                                   
#                                   if (family %in% c('bernouli','binomial')) {
#                                     
#                                     np <- length(x[x == 1])
#                                     na <- length(x) - np
#                                     min.n <- min(c(np,na))
#                                     
#                                     if (min.n == 0) stop('The distribution is not binomial!')
#                                     if (is.null(k)) {
#                                       if (min.n > 5) k <- 5
#                                       else k <- min.n
#                                     } else if (k > min.n) {
#                                       k <- min.n
#                                       warning(paste('number of folds cannot be greater than the number of either presence or absence records. nfolds is changed to:',min.n))
#                                     }
#                                     if (!st) {
#                                       o[[1]] <- rep(seq_len(k), length.out=length(x))
#                                       o[[2]] <- replicate(r, sample(1:length(x)))
#                                       
#                                       if (!all(unlist(lapply(unique(o[[1]]),function(i) {
#                                         unlist(lapply(1:r,function(j) length(unique(x[o[[2]][o[[1]] == i,j]])) == 2))
#                                       })))) {
#                                         st <- TRUE
#                                         warning('To avoid error in partitioning, stratify is turned to TRUE')
#                                       }
#                                     }
#                                     if (st) {
#                                       o[[1]] <- c(rep(seq_len(k), length.out=np),rep(seq_len(k), length.out=na))
#                                       o[[2]] <- rbind(replicate(r, sample(which(x == 1))),replicate(r, sample(which(x == 0))))
#                                     }
#                                   } else {
#                                     if (k > length(x)) {
#                                       k <- length(x)
#                                       warning('nfolds cannot be greater than the number of observations! It is changed to the number of observations (became "leave one out"/Jackknife procedure)')
#                                     }
#                                     o[[1]] <- rep(seq_len(k), length.out=length(x))
#                                     o[[2]] <- replicate(r, sample(1:length(x)))
#                                   }
#                                   o
#                                 })
# 
# 
# .replicateMethods$Methods['bootstrap'] <- new('.methodTemplate',
#                                     name='bootstrap',
#                                     aliases=c('boot','bootstrapping','bootstraping'),
#                                     arguments=c(x='numeric',replicates='numeric'),
#                                     Help='Partition data into train and test using bootstrapping.',
#                                     Function=function(x,...) {
#                                       o <- list()
#                                       r <- list(...)[['replicates']]
#                                       family <- list(...)[['family']]
#                                       if (is.null(r)) r <- 1
#                                       if (is.null(family)) {
#                                         if (all(sort(unique(x)) == c(0,1))) family <- 'binomial'
#                                       }
#                                       .sampleboot <- function(x,family) {
#                                         .s <- function(x) {
#                                           o <- list()
#                                           n <- length(x)
#                                           s <- sample(c(1:n),n,replace=TRUE)
#                                           w <- c(1:n) %in% s
#                                           if (length(which(w)) == n) {
#                                             # if all records are selected for train (i.e., there is no records remained for test), it is repeated to avoid error (up to 100 times)
#                                             for (i in 1:100) {
#                                               s <- sample(1:n,n,replace=TRUE)
#                                               w <- c(1:n) %in% s
#                                               if (length(which(w)) < n) break
#                                             }
#                                           }
#                                           o[[1]] <- ifelse(w,1,2)
#                                           o[[2]] <- s
#                                           o
#                                         }
#                                         o <- .s(x)
#                                         if (family %in% c('bernouli','binomial')) {
#                                           w1 <- length(unique(x[o[[2]][o[[1]] == 1]])) == 2
#                                           w2 <- length(unique(x[o[[2]][o[[1]] == 2]])) == 2
#                                           if (!c(w1 & w2)) {
#                                             for (i in 1:100) {
#                                               o <- .s(x)
#                                               w1 <- length(unique(x[o[[2]][o[[1]] == 1]])) == 2
#                                               w2 <- length(unique(x[o[[2]][o[[1]] == 2]])) == 2
#                                               if (c(w1 & w2)) break
#                                             }
#                                           }
#                                         }
#                                         o
#                                       }
#                                       
#                                       w <- replicate(r,.sampleboot(x,family=family),simplify=FALSE)
#                                       o[[1]] <- list()
#                                       o[[2]] <- matrix(ncol=r,nrow=length(x))
#                                       for (i in 1:r) {
#                                         o[[1]][[i]] <- w[[i]][[1]]
#                                         o[[2]][,i] <- w[[i]][[2]]
#                                       }
#                                       o
#                                     })
# #-----
# .replicateMethods$Methods['subsampling'] <- new('.methodTemplate',
#                                 name='subsampling',
#                                 aliases=c('sub','subsample'),
#                                 arguments=c(x='numeric',replicates='numeric',test.percent='numeric',stratify='logical',family='character'),
#                                 Help='Subsampling procedure randomly draw (without replacement) a proportion of data as test dataset, and keep the remained records for training',
#                                 Function=function(x,...) {
#                                   o <- list()
#                                   r <- list(...)[['replicates']]
#                                   testp <- list(...)[['test.percent']]
#                                   st <- list(...)[['stratify']]
#                                   family <- list(...)[['family']]
#                                   if (is.null(r)) r <- 1
#                                   if (is.null(st)) st <- TRUE
#                                   if (is.null(family)) {
#                                     if (all(sort(unique(x)) == c(0,1))) family <- 'binomial'
#                                   }
#                                   if (is.null(testp)) testp <- 0.3
#                                   else if (testp > 100) {
#                                     warning('test percentage should be less than 100, it is changed to the default, 30!')
#                                     testp <- 0.3
#                                   }
#                                   
#                                   if (testp > 1) testp <- testp / 100
#                                   if (family %in% c('bernouli','binomial')) {
#                                     np <- length(x[x == 1])
#                                     na <- length(x) - np
#                                     npt <- ceiling(np * testp)
#                                     nat <- ceiling(na * testp)
#                                     if (np == npt) {
#                                       if (np > 1) {
#                                         npt <- npt - 1
#                                         st <- TRUE
#                                       } else stop('Subsampling is not possible with only 1 presence record...!')
#                                     }
#                                     if (na == nat) {
#                                       if (na > 1) {
#                                         nat <- nat - 1
#                                         st <- TRUE
#                                       } else stop('Subsampling is not possible with only 1 absence record...!')
#                                     }
#                                     if (!st) {
#                                       test.nr <- ceiling(length(x) * testp)
#                                       o[[1]] <- rep(1:2,times=c((length(x)-test.nr),test.nr))
#                                       o[[2]] <- replicate(r, sample(1:length(x)))
#                                       u <- unlist(lapply(unique(o[[1]]),function(i) {
#                                         unlist(lapply(1:r,function(j) length(unique(x[o[[2]][o[[1]] == i,j]])) == 2))
#                                       }))
#                                       if (!all(u)) {
#                                         st <- TRUE
#                                         warning('To avoid error in partitioning, stratified partitioning is used!')
#                                       }
#                                     }
#                                     
#                                     if (st) {
#                                       o[[1]] <- c(rep(1:2,times=c(np-npt,npt)),rep(1:2,times=c(na-nat,nat)))
#                                       o[[2]] <- rbind(replicate(r, sample(which(x == 1))),replicate(r, sample(which(x == 0))))
#                                     } 
#                                   } else {
#                                     test.nr <- ceiling(length(x) * testp)
#                                     o[[1]] <- rep(1:2,times=c((length(x)-test.nr),test.nr))
#                                     o[[2]] <- replicate(r, sample(1:length(x)))
#                                   } 
#                                   o
#                                 })
# 
