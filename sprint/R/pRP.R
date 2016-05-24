##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008-2011 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.org/licenses/>.  #
#                                                                        #
##########################################################################

`pRS` <- function(data, cl, num.perm, logged=TRUE, na.rm=FALSE, gene.names=NULL,
                  plot=FALSE, rand=NULL)
  {
    pRP(data, cl, num.perm, logged, na.rm, gene.names, plot, rand,
        sum=TRUE)
  }

`pRP` <- function(data, cl, num.perm, logged=TRUE, na.rm=FALSE,
                  gene.names=NULL, plot=FALSE, rand=NULL,
                  sum=FALSE)
  {
    init.rng(seed=rand)
    nsamples <- length(cl)
    if ( nsamples != ncol(data) )
      stop("Number of classes should match the number of samples (columns) in the data")

    ngene <- nrow(data)

    xy <- xyCall(cl)
    x <- xy$x
    y <- xy$y

    if ( is.null(gene.names) )
      gene.names <- rownames(data)

    data <- as.matrix(data)

    mode(data) <- "numeric"
    if ( any(is.na(data)) ) {
      if ( !na.rm ) {
        stop("Cannot deal with NA values in data.  Try with na.rm=TRUE")
      }
      data <- na.replace(data)
    }
    if ( any(is.nan(data)) ) {
      stop("Data contain NaN values, please remove them")
    }
    data1 <- data[,x]                   # class 1 samples
    data2 <- data[,y]                   # class 2 samples
    nclass1 <- ncol(data1)
    nclass2 <- ncol(data2)

    data1.mean <- apply(data1, 1, mean)
    data2.mean <- apply(data2, 1, mean)

    if ( nclass2 == 0 ) {
      fold.change <- data1.mean
    } else {
      if ( logged ) {
        fold.change <- data1.mean - data2.mean
      } else {
        fold.change <- data1.mean/data2.mean
      }
    }

    RP.experiment.upin2 <- rankprod(data1, data2, logged,
                                    rev.sorting=FALSE, sum)
    rank.experiment.upin2 <- rank(RP.experiment.upin2$RP)

    RP.experiment.downin2 <- rankprod(data1, data2, logged,
                                     rev.sorting=TRUE, sum)
    rank.experiment.downin2 <- rank(RP.experiment.downin2$RP)

    RP.perm.upin2 <- bootRP(data1, data2, RP.experiment.upin2$RP,
                            logarithmic.data=logged,
                            rev.sorting=FALSE, nperms=num.perm, sum)$bootRP

    RP.perm.downin2 <- bootRP(data1, data2, RP.experiment.downin2$RP,
                              logarithmic.data=logged,
                              rev.sorting=TRUE, nperms=num.perm, sum)$bootRP
    ## C bootstrapping counts how many times the experimental RP is
    ## smaller than the bootstrapped RP, so exp.count is just that.
    ## We do the counting in C so as not to return a big array to R:
    ## we return a vector the size of the number of genes, rather than
    ## nperm * ngene.
    exp.count <- RP.perm.upin2/num.perm

    pval.upin2 <- RP.perm.upin2/(num.perm * ngene)

    pfp.upin2 <- exp.count/rank.experiment.upin2

    exp.count <- RP.perm.downin2/num.perm

    pval.downin2 <- RP.perm.downin2/(num.perm * ngene)

    pfp.downin2 <- exp.count/rank.experiment.downin2

    pfp <- data.frame(pfp.upin2, pfp.downin2)

    colnames(pfp) <- c("class1 < class2", "class1 > class2")
    pval <- data.frame(pval.upin2, pval.downin2)
    colnames(pval) <- c("class1 < class2", "class1 > class2")

    RPs <- data.frame(RP.experiment.upin2$RP, RP.experiment.downin2$RP)
    colnames(RPs) <- c("class1 < class2", "class1 > class2")

    RPrank <- data.frame(rank.experiment.upin2, rank.experiment.downin2)
    colnames(RPrank) <- c("class1 < class2", "class1 > class2")

    fold.change <- t(t(fold.change))
    colnames(fold.change) <- "log/unlog(class1/class2)"

    rownames(pfp) <- gene.names
    rownames(pval) <- gene.names
    rownames(RPs) <- gene.names
    rownames(RPrank) <- gene.names
    rownames(fold.change) <- gene.names

    ret <- list(pfp=pfp,
                pval=pval,
                RPs=RPs,
                RPrank=RPrank,
                AveFC=fold.change)

    if ( plot ) {
      if( !require("RankProd", quietly=TRUE) ) {
        warning("Cannot plot the rank product result - failed to load package \"RankProd\". 
                Please check that the package is installed and try again or run function with \'plot=FALSE\'.")
        return(NA)
        plotRP(ret, cutoff=NULL)
       }
    }
    
    reset.rng()
    return(ret)
  }

`pRSadvance` <- function(data, cl, origin, num.perm=100,
                         logged=TRUE, na.rm=FALSE, gene.names=NULL,
                         plot=FALSE, rand=NULL)
  {
    pRPadvance(data, cl, origin, num.perm, logged, na.rm, gene.names,
               plot, rand, sum=TRUE)
  }

`pRPadvance` <- function(data, cl, origin, num.perm=100,
                         logged=TRUE, na.rm=FALSE, gene.names=NULL,
                         plot=FALSE, rand=NULL, sum=FALSE)
  {
    init.rng(seed=rand)
    num.origins <- length(unique(origin))
    ngene <- nrow(data)

    if ( any(is.na(data)) ) {
      if ( !na.rm ) {
        stop("Cannot deal with NA values in data.  Try with na.rm=TRUE")
      }
      for ( i in unique(origin) ) {
        tmp <- matrix(data[, origin == i], nrow=ngene)
        data[, origin == i ] <- na.replace(tmp)
      }
    }
    if ( any(is.nan(data)) ) {
      stop("Data contain NaN values, please remove them")
    }

    data.pre <- xyCall.multi(data, cl, origin)

    data1.all <- data.pre$data1
    data2.all <- data.pre$data2


    if ( is.null(gene.names) )
      gene.names <- rownames(data)

    ## Ensure all list entries are in matrix form.
    for ( i in 1:num.origins ) {
      data1.all[[i]] <- matrix(data1.all[[i]], nrow=ngene)
      if ( !is.null(data2.all[[i]]) ) {
        data2.all[[i]] <- matrix(data2.all[[i]], nrow=ngene)
      }
    }

    if ( is.null(data2.all) ) {
      fc <- sapply(data1.all, function(x) apply(x, 1, mean))
      ave.fold.change <- apply(fc, 1, mean)
    } else {
      d1 <- sapply(data1.all, function(x) apply(x, 1, mean))
      d2 <- sapply(data2.all, function(x) apply(x, 1, mean))
      if ( logged ) {
        fc <- d1 - d2
      } else {
        fc <- d1/d2
      }
      ave.fold.change <- apply(fc, 1, mean)
    }
    RP.experiment.upin2 <- rankprod.multi(data1.all, data2.all,
                                          num.origins, logged,
                                          rev.sorting=FALSE, sum)
    rank.experiment.upin2 <- rank(RP.experiment.upin2$RP)
    RP.experiment.downin2 <- rankprod.multi(data1.all, data2.all,
                                          num.origins, logged,
                                          rev.sorting=TRUE, sum)
    rank.experiment.downin2 <- rank(RP.experiment.downin2$RP)

    RP.perm.upin2 <- bootRP.multi(data1.all, data2.all,
                                  RP.experiment.upin2$RP,
                                  num.origins,
                                  logarithmic.data=logged,
                                  rev.sorting=FALSE,
                                  nperms=num.perm, sum)
    RP.perm.downin2 <- bootRP.multi(data1.all, data2.all,
                                    RP.experiment.downin2$RP,
                                    num.origins,
                                    logarithmic.data=logged,
                                    rev.sorting=TRUE,
                                    nperms=num.perm, sum)$bootRP

    exp.count <- RP.perm.upin2$bootRP/num.perm

    pval.upin2 <- RP.perm.upin2$bootRP/(num.perm * ngene)

    pfp.upin2 <- exp.count/rank.experiment.upin2

    exp.count <- RP.perm.downin2/num.perm

    pval.downin2 <- RP.perm.downin2/(num.perm * ngene)

    pfp.downin2 <- exp.count/rank.experiment.downin2

    pfp <- data.frame(pfp.upin2, pfp.downin2)

    colnames(pfp) <- c("class1 < class2", "class1 > class2")
    pval <- data.frame(pval.upin2, pval.downin2)
    colnames(pval) <- c("class1 < class2", "class1 > class2")

    RPs <- data.frame(RP.experiment.upin2$RP, RP.experiment.downin2$RP)
    colnames(RPs) <- c("class1 < class2", "class1 > class2")

    RPrank <- data.frame(rank.experiment.upin2, rank.experiment.downin2)
    colnames(RPrank) <- c("class1 < class2", "class1 > class2")

    ave.fold.change <- t(t(ave.fold.change))
    colnames(ave.fold.change) <- "log/unlog(class1/class2)"

    colnames(fc) <- paste("(class1/class2)-data", 1:num.origins, sep="")

    rownames(pfp) <- gene.names
    rownames(pval) <- gene.names
    rownames(RPs) <- gene.names
    rownames(RPrank) <- gene.names
    rownames(ave.fold.change) <- gene.names
    rownames(fc) <- gene.names

    ret <- list(pfp=pfp,
                pval=pval,
                RPs=RPs,
                RPrank=RPrank,
                AveFC=ave.fold.change,
                all.FC=fc)

    if ( plot ) {
      if( !require("RankProd", quietly=TRUE) ) {
        warning("Cannot plot the rank product result - failed to load package \"RankProd\". 
                Please check that the package is installed and try again or run function with \'plot=FALSE\'.")
        return(NA)
        plotRP(ret, cutoff=NULL)
      }
    }

    reset.rng()
    return(ret)
  }

`rankprod.multi` <- function(data1, data2, num.origins, logged, rev.sorting,
                             sum=FALSE)
  {

    d1 <- NULL
    d2 <- NULL
    ngenes <- nrow(data1[[1]])
    for ( i in 1:num.origins ) {
      d1 <- cbind(d1, data1[[i]])
      d2 <- cbind(d2, data2[[i]])
    }
    nclass1 <- sapply(data1, ncol)
    if ( is.null(data2) ) {
      nclass2 <- NULL
    } else {
      nclass2 <- sapply(data2, ncol)
    }

    ret <- .C("rank_product_multi_",
              as.double(d1),
              as.integer(nclass1),
              as.double(d2),
              as.integer(nclass2),
              as.integer(ngenes),
              as.integer(num.origins),
              as.integer(logged),
              as.integer(rev.sorting),
              as.integer(sum),
              RP = double(ngenes))

    list(RP=ret$RP)
  }

`rankprod` <- function(data1, data2, logged, rev.sorting, sum=FALSE)
  {
    nclass1 <- ncol(data1)
    nclass2 <- ncol(data2)
    ngenes <- nrow(data1)
    ret <- .C("rank_product_",
              as.double(data1),
              as.integer(nclass1),
              as.double(data2),
              as.integer(nclass2),
              as.integer(ngenes),
              as.integer(logged),
              as.integer(rev.sorting),
              as.integer(sum),
              RP = double(ngenes))

    list(RP=ret$RP)
  }

`bootRP.multi` <- function(data1, data2, experimental.rp,
                           num.origins,
                           logarithmic.data=FALSE,
                           rev.sorting=FALSE, nperms=100, sum=FALSE)
  {
    d1 <- NULL
    d2 <- NULL
    ngenes <- nrow(data1[[1]])
    for ( i in 1:num.origins ) {
      d1 <- cbind(d1, data1[[i]])
      data1[[i]] <- matrix(data1[[i]], nrow=ngenes)
      if ( !is.null(data2) ) {
        d2 <- cbind(d2, data2[[i]])
        data2[[i]] <- matrix(data2[[i]], nrow=ngenes)
      }
    }
    nclass1 <- sapply(data1, ncol)
    if ( is.null(data2) ) {
      nclass2 <- NULL
    } else {
      nclass2 <- sapply(data2, ncol)
    }

    ret <- .C("boot_rank_product_multi_",
              as.double(d1),
              as.integer(nclass1),
              as.double(d2),
              as.integer(nclass2),
              as.double(experimental.rp),
              as.integer(ngenes),
              as.integer(num.origins),
              as.integer(logarithmic.data),
              as.integer(rev.sorting),
              as.integer(sum),
              as.integer(nperms),
              bootRP = double(ngenes),
              exitcode = double(1))

    if ( ret$exitcode == -1 ) {
      warning(paste("MPI is not initialized. bootRP aborted.\n"))
      return(NULL)
    }

    return(list(bootRP=ret$bootRP))

  }
`bootRP` <- function(data1, data2, experimental.rp, logarithmic.data=FALSE,
                     rev.sorting=FALSE, nperms=100, sum=FALSE)
  {
    nclass1 <- ncol(data1)
    nclass2 <- ncol(data2)
    ngenes <- nrow(data1)
    ret <- .C("boot_rank_product_",
              as.double(data1),
              as.integer(nclass1),
              as.double(data2),
              as.integer(nclass2),
              as.double(experimental.rp),
              as.integer(ngenes),
              as.integer(logarithmic.data),
              as.integer(rev.sorting),
              as.integer(sum),
              as.integer(nperms),
              bootRP = double(ngenes),
              exitcode = double(1))

    if ( ret$exitcode == -1 ) {
      warning(paste("MPI is not initialized. bootRP aborted.\n"))
      return(NULL)
    }

    return(list(bootRP=ret$bootRP))
  }

`na.replace` <- function(data)
  {
    na.genes <- which(is.na(data), arr.ind=T)
    data[na.genes] <- apply(matrix(data[na.genes[,1],], ncol=ncol(data)),
                            1, mean, na.rm=T)
    return(data)
  }

`xyCall` <- function(classes)
  {
    levels <- unique(classes)
    nclass <- length(levels)
    if ( nclass != 2 & nclass != 1 ) {
      stop("Invalid class labels specified, data must have 1 or 2 classes")
    }

    if ( nclass == 1 ) {
      x <- which(classes == levels)
      y <- NULL
    } else {
      x <- which(classes == min(levels))
      y <- which(classes == max(levels))
    }

    return(list(x=x, y=y))
  }

`xyCall.multi` <- function(data, classes, origins)
  {
    levels <- unique(classes)
    nclass <- length(levels)

    origin.levels <- unique(origins)
    norigin <- length(origin.levels)

    if ( nclass != 1 & nclass != 2 ) {
      stop("Invalid class labels specified, data must have 1 or 2 classes")
    }

    if ( nclass == 1 ) {
      data1 <- vector(mode="list", length=norigin)
      data2 <- NULL
      for ( i in 1:norigin ) {
        data1[[i]] <- data[, origins == origin.levels[i]]
      }
    } else {
      data1 <- vector(mode="list", length=norigin)
      data2 <- vector(mode="list", length=norigin)
      for ( i in 1:norigin ) {
        idx1 <- which(origins == origin.levels[i] & classes == min(levels))
        idx2 <- which(origins == origin.levels[i] & classes == max(levels))

        if ( length(idx1) == 0 | length(idx2) == 0 ) {
          stop(paste("Two class analysis specified but data from origin",
                     i,
                     "only has data from one class"))
        }
        data1[[i]] <- data[,idx1]
        data2[[i]] <- data[,idx2]
      }
    }

    return(list(data1=data1, data2=data2))
  }
