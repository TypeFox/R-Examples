itegeppXXR <-
  function(geno, des = 0, lim = 0.05)
  {
    if (!is.matrix(geno)) 
      geno <- as.matrix(geno)
    if ((!is.numeric(lim)) || (lim < 0.0) || (lim > 1.0)) 
      stop("Bad value for lim.")
    # gives out the design matrix for haplotypes(des=0) or
    # haplotype pairs(des=1, default)
    if(is.na(match(des, c(0, 1))))
      stop("des should be 0 or 1.")
    # max 16 SNPs
    ns <- dim(geno)[2]
    if(ns > 16) 
      stop("Number of SNPs should smaller than 17.")
    # replace NA as 0
    geno <- replace(geno, is.na(geno), 0)
    genotyp <- apply(geno, 1, paste, collapse = "")
    len <- dim(geno)[1]
    wid <- ns
    lest <- paste(rep(" ", 10000), collapse = "")
    dr <- rep(lest, len)
    fr <- rep(lest, (2^(ns)) + 1)
    hr <- rep(lest, len)
    # call C function
    zzz <- .C("itegeppXXR", 
              as.integer(des),
              as.double(lim), 
              as.character(genotyp), 
              as.double(1), 
              as.integer(len), 
              likeres = as.double(0.1, 0), 
              freqres = as.character(fr), 
              hapres = as.character(hr), 
              desres = as.character(dr), 
              PACKAGE = "HapEstXXR")
    hr <- zzz$hapres[zzz$hapres != lest]
    zzz$freqres <- as.matrix(zzz$freqres[zzz$freqres != lest])
    hp <- as.matrix(apply(zzz$freqres, 1, 
                          function(x) unlist(strsplit(x, split=" "))))
    if(dim(hp)[1]==0)
    {
      print(paste("Error in itegeppXXR: ",
                  "all inferred haplotypes with probability ", 
                  "below threshold(lim = ", lim, ")",  sep = ""))
      # abort of itegeppXXR
      return(list(hap = NA, freq = NA, hapres = NA, 
                  likeres = NA, desres = NA))
    }
    em.hap <- hp[1, , drop = FALSE]
    nhap <- length(em.hap)
    em.freq <- hp[2, ]
    # function to split for more than one blank
    fsplit <- function(x) {
      x <- unlist(strsplit(x, split= "[[:space:]]"))
      x <- as.numeric(x[nchar(x) > 0])
      return(x)
    }
    zzz$desres <- t(apply((as.matrix(zzz$desres)), 
                          1, 
                          fsplit))
    # Spalten|berschriften
    if(des == 0)
    { # haplotypes
      if(ncol(zzz$desres) == nhap)
      {
        colnames(zzz$desres) <- em.hap
      } else
      {
        #print("ACHTUNG seltene Haplotypen...")
        colnames(zzz$desres) <- c(em.hap, "R")
      }
    } else
    { # des=1   haplotype pairs
      nn <- rep(0, (nhap * (nhap + 1)) / 2)
      k <- 1
      for(i in 1:nhap)
      {
        for(j in i:nhap)
        {
          nn[k] <- paste(em.hap[i], em.hap[j], sep = ".")
          k <- k + 1
        }
      }
      if(ncol(zzz$desres) == length(nn))
      {
        colnames(zzz$desres) <- nn
      } else
      {
        #print("ACHTUNG seltene Haplotypen...")
        colnames(zzz$desres) <- c(nn, "R")
      }
    }
    zzz$desres <- zzz$desres[, colSums(zzz$desres) > 0, drop = FALSE]
    rownames(zzz$desres) <- NULL
    return(list(hap = as.character(em.hap), 
                freq = as.numeric(em.freq), 
                hapres = hr, 
                likeres = zzz$likeres, 
                desres = as.matrix(zzz$desres)))
  }
