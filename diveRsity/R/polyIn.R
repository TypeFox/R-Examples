################################################################################
# polyIn                                                                       #
################################################################################
#' A function for calculating informativeness for the inference of ancestry
#' for loci of any ploidy.
#' Kevin Keenan 2014
#' @export 
polyIn <- function(infile = NULL, pairwise = FALSE, para = FALSE){
  if(is.null(infile)){
    stop("Please provide and input file!")
  }
  polyReader <- function(infile, para){
    # read data
    fastScan <- function(fname) {
      s <- file.info(fname)$size
      buf <- readChar(fname, s, useBytes = TRUE)
      # replace Mac encoded line endings
      if(length(grep("\r", buf)) != 0L){
        buf <- gsub("\r", "\n", buf)
        buf <- gsub("\n\n", "\n", buf)
      }
      return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
    }
    dat <- fastScan(infile)
    if(length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 1){
      locs <- strsplit(dat[2], split = "\\s+")[[1]]
      if(length(locs) == 1){
        locs <- strsplit(dat[2], split = ",")[[1]]
      }
      locs <- as.character(sapply(locs, function(x){
        x <- strsplit(x, split = "")[[1]]
        if(is.element(",", x)){
          x <- x[-(which(x == ","))]
        }
        return(paste(x, collapse = ""))
      }))
      dat <- c(dat[1], locs, dat[-(1:2)])
    }
    # npops
    popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
    npops <- length(popLoc)
    no_col <- popLoc[1] - 1
    # get genotypes
    strt <- popLoc + 1
    ends <- c(popLoc[-1] - 1, length(dat))
    genoRet <- function(strt, ends, x){
      out <- strsplit(x[strt:ends], split = "\\s+")
      return(do.call("rbind", out))
    }
    genos <- mapply(genoRet, strt = strt, ends = ends, 
                    MoreArgs = list(x = dat), SIMPLIFY = FALSE)
    indNames <- lapply(genos, function(x){
      return(x[,1])
    })
    indNames <- do.call("c", indNames)
    genos <- lapply(genos, function(x){
      x[x == "-9"] <- NA
      return(x[,-1])
    })
    # ploidy estimation
    ploidy <- round(mean(nchar(na.omit(genos[[1]][,1]))))
    # convert genotypes into allele arrays
    if(para){
      
      cl <- parallel::makeCluster(detectCores())
      parallel::clusterExport(cl, "ploidy", envir = environment())
      allArr <- parallel::parLapply(cl, genos, function(x){
        alls <- strsplit(x, split = "")
        # fix NAs
        alls <- lapply(alls, function(y){
          if(any(is.na(y))){
            y <- rep(NA, ploidy)
            return(y)
          } else {
            return(y)
          }
        })
        # alls array
        alls <- matrix(alls, ncol = ncol(x), nrow = nrow(x))
        alls <- lapply(1:ncol(alls), function(i){
          do.call("rbind", alls[,i])
        })
        af <- lapply(alls, function(al){
          ct <- table(al)
          return(as.vector(ct/sum(ct)))
        })
      })
      # create allele frequency matrices
      afMat <- parallel::parLapply(cl, 1:length(locs), function(i){
        dat <- lapply(allArr, "[[", i)
        return(do.call("cbind", dat))
      })
      parallel::stopCluster(cl)
    } else {
      allArr <- lapply(genos, function(x){
        alls <- strsplit(x, split = "")
        # fix NAs
        alls <- lapply(alls, function(y){
          if(any(is.na(y))){
            y <- rep(NA, ploidy)
            return(y)
          } else {
            return(y)
          }
        })
        # alls array
        alls <- matrix(alls, ncol = ncol(x), nrow = nrow(x))
        alls <- lapply(1:ncol(alls), function(i){
          do.call("rbind", alls[,i])
        })
        af <- lapply(alls, function(al){
          ct <- table(al)
          return(as.vector(ct/sum(ct)))
        })
      })
      # create allele frequency matrices
      afMat <- lapply(1:length(locs), function(i){
        dat <- lapply(allArr, "[[", i)
        return(do.call("cbind", dat))
      })
    }
    
    names(afMat) <- locs
    return(afMat)
  }
  
  inFunc <- function(af, pw = FALSE){
    if(pw){
      combs <- combn(ncol(af), 2)
      afComb <- lapply(1:ncol(combs), function(i){
        return(af[,combs[,i]])
      })
      Out <- lapply(afComb, function(x){
        sum(apply(x, 1, function(y){
          trm1 <- -mean(y) * log(mean(y))
          trm2 <- sum((y/length(y))*log(y))
          return(sum(sum(trm1 + trm2)))
        }))
      })
      inOut <- matrix(NA, ncol = ncol(af), nrow = ncol(af))
      for(i in 1:length(Out)){
        inOut[combs[2,i], combs[1,i]] <- Out[[i]]
      }
    } else {
      inOut <- apply(af, 1, function(x){
        trm1 <- -mean(x) * log(mean(x))
        trm2 <- sum((x/length(x))*log(x))
        return(sum(sum(trm1 + trm2)))
      })
      inOut <- sum(inOut)
    }
    return(inOut)
  }
  afs <- polyReader(infile, para)
  locNames <- names(afs)
  if(para){
    
    cl <- parallel::makeCluster(detectCores())
    parallel::clusterExport(cl, c("inFunc", "pairwise"), envir = environment())
    ins <- parallel::parLapply(cl, afs, inFunc, pw = pairwise)
    parallel::stopCluster(cl)
    
  } else {
    ins <- lapply(afs, inFunc, pw = pairwise)
  }
  if(!pairwise){
    ins <- unlist(ins)
    names(ins) <- locNames
  } else {
    names(ins) <- locNames
  }
  return(ins)
}
################################################################################
# polyIn                                                                       #
################################################################################