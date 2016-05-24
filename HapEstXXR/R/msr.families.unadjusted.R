msr.families.unadjusted  <- 
  function(famid, patid, fid, mid, trait, snps, 
           pair.begin = TRUE, lim = 0.05, maxSNP = 3, nt = 10) {
    
    nloc <- dim(snps)[2]
    nind <- dim(snps)[1]
    snps <- as.matrix(snps)
    trait <- as.numeric(trait)
    if (!all(trait[!is.na(trait)] == 0 | trait[!is.na(trait)] == 1))
      stop("trait should be 0 for unaffected or 1 for affcted children.")
    ### error
    if (maxSNP > nloc) 
      stop("Error in stepwise.fam: maxSNP>nloc(dim(snps)[2]) ")
    if (length(famid) != nind) 
      stop("Error in stepwise.fam: unexpected length of famid")
    if (length(patid) != nind) 
      stop("Error in stepwise.fam: unexpected length of patid")
    if (length(fid) != nind) 
      stop("Error in stepwise.fam: unexpected length offid")
    if (length(mid) != nind) 
      stop("Error in stepwise.fam: unexpected length of mid")
    if (length(trait) != nind) 
      stop("Error in stepwise.fam: unexpected length of trait")
    if ((lim < 0) || (lim > 1)) 
      stop("Error in stepwise.fam: 0 <= lim <= 1.")
    # remove unaffected children
    not.rm.child <- ifelse((fid != 0) & (mid != 0) & (trait == 0), 
                           FALSE, TRUE)
    if (sum(!not.rm.child)>0) {
      cat(sum(!not.rm.child),
          " children removed because they were unaffected.\n",
          sep="")
    }
    famid <- famid[not.rm.child]
    patid <- patid[not.rm.child]
    fid   <- fid[not.rm.child]
    mid   <- mid[not.rm.child]
    snps  <- snps[not.rm.child, , drop = FALSE]
    trait <- trait[not.rm.child]
    if (length(famid) < 1) {
      stop("no families for TDT observed.")
    }
    # exclusion of nuclear families without two parents.
    excl.fam <- NULL
    #i <- unique(famid)[1]
    for(i in unique(famid)) {
      selfam <- famid == i
      selfid <- unique(fid[selfam])
      selfid <- selfid[(!is.na(selfid)) & (selfid != 0) ]
      selmid <- unique(mid[selfam])
      selmid <- selmid[(!is.na(selmid)) & (selmid != 0) ]
      if (length(selfid) != 1) {
        excl.fam <- c(excl.fam, i)
      } else {
        if (! any(patid[selfam] %in% selfid)) {
          excl.fam <- c(excl.fam, i)
        }
      }
      if (length(selmid) != 1) {
        excl.fam <- c(excl.fam, i)
      } else {
        if (! any(patid[selfam] %in% selmid)) {
          excl.fam <- c(excl.fam, i)
        }
      }
    }
    excl.fam <- unique(excl.fam)
    if (!is.null(excl.fam)) {
      print(paste("Exclusion of nuclear families without two parents: ", 
                  paste(excl.fam, collapse=" "), 
                  sep=""))
      selcond <-  !(famid %in% excl.fam)
      famid <- famid [selcond]
      patid <- patid [selcond]
      fid   <- fid   [selcond]
      mid   <- mid   [selcond]
      snps  <- snps  [selcond, , drop = FALSE]
      trait <- trait [selcond]
    }
    if (length(famid) < 1) {
      stop("no families for TDT observed.")
    }
    test <- function(x, a) { 
      all(x %in% a)
    }    
    nloc <- dim(snps)[2]
    lenx <- length(famid);
    lest <- paste(rep(" ", 1000), collapse = "")
    pr <- rep(lest, 1)
    res.list <- list()
    pval <- rep(NA, choose(nloc, 2))
    pos <- matrix(rep(0, choose(nloc, 2)*2), ncol=2)
    k <- 1
    if (pair.begin == FALSE) {
      # begin with single SNPs
      cat(paste("Iteration with 1 SNP  at same time.   System.time = ", Sys.time(), "\n", sep=""))
      rest <- single.snp.test.families(snps, 
                                       trait, 
                                       adj.var = NULL, 
                                       famid, 
                                       patid, 
                                       fid, 
                                       mid, 
                                       prt = FALSE)
      ind  <- (order(rest[, "p.value"]))[1:min(nt, length(rest[, "p.value"]))]
      pval <- as.numeric(rest[, "p.value"])[ind]
      pos  <- (rest[, "SNP"]) [ind]
      i <- 1
      res.list[[i]] <- data.frame(pos, 
                                  rep("families", length(pval)),
                                  pval, 
                                  stringsAsFactors = FALSE)
      rownames(res.list[[i]]) <- NULL
      colnames(res.list[[i]]) <- c(paste("SNP", 1:i, sep = ""),
                                   "type",
                                   "p.value")
      k <- k + 1
    } else {
      ### pair begin ###
      cat(paste("Iteration with 2 SNPs at same time.   System.time = ",
                Sys.time(), "\n", sep=""))
      cat(paste("Number of SNP pairs = ", choose(nloc, 2), "\n", sep=""))
      for(i in 1:(nloc - 1)) {
        for(j in (i + 1):nloc) {
          if ((k %% 5000) == 0) {
            cat(paste("Step =  ", k, 
                      "   System.time = ", (Sys.time()), "\n", sep=""))
          }
          xgeno <- snps[, c(i, j)];
          geno <- apply(xgeno, 1, paste, collapse = "");          
          zzz <- .C("haptdpn", 
                    as.character(famid),
                    as.character(patid), 
                    as.character(geno), 
                    as.integer((as.integer(trait) + 1)),
                    as.integer(lenx), 
                    as.double(lim), 
                    pvres = pr)[[7]];
          
          x <- strsplit(zzz, " ");
          x[[1]] <- x[[1]][x[[1]] != ""]
          x5 <- as.numeric(x[[1]][4])
          x1 <- as.numeric(x[[1]][1])
          if (x5 == 1){
            pv <- 1.0 - pchisq(x5, 1)
          } else {
            pv <- 1.0 - pchisq((x5 - 1) * x1 / x5, x5 - 1)
          }
          # save result
          pval[k] <- pv
          pos[k, ] <- c(i, j);
          k <- k + 1
        } # end of for j
      } # end of for i
      ind  <- (order(pval))[1:min(nt, length(pval))  ]
      pval <- pval[ind]
      pos <- pos [ind, , drop = FALSE]
      i <- 2
      res.list[[i]] <- data.frame(pos,
                                  rep("families", length(pval)), 
                                  pval, 
                                  stringsAsFactors = FALSE)
      rownames(res.list[[i]]) <- NULL
      colnames(res.list[[i]]) <- c(paste("SNP", 1:i, sep = ""),
                                   "type",
                                   "p.value")
    } # end of pair.begin
    i <- i + 1
    while((i <= nloc) && (i <= maxSNP)) {
      cat(paste("Iteration with ", i, 
                " SNPs at same time.   System.time = ", 
                Sys.time(), "\n", sep = ""))      
      # create data set with all possible SNP combinations
      BestPos <- as.matrix(res.list[[i - 1]] [, 1:(i - 1), drop = FALSE])
      storage.mode(BestPos) <- "integer"
      newdim <- as.integer(c(dim(BestPos)[1] * (nloc-1), dim(BestPos)[2] + 1))
      out <- .C("create_pattern_matrix", 
                pattern   = as.integer(BestPos),
                ndim      = dim(BestPos), 
                snps      = as.integer(1:nloc), 
                snplen    = nloc, 
                newpat    = as.integer(rep(0, newdim[1] * newdim[2])), 
                newpatdim = newdim, 
                len       = as.integer(0))
      BestPos  <- (matrix(out$newpat, 
                          nrow = newdim[1],
                          ncol = newdim[2],
                          byrow = FALSE))[1:out$len, , drop = FALSE]
      cat(paste("Iteration with ", i, " SNPs at same time. ---- ", 
                dim(BestPos)[1], " detected SNP combinations -----\n", 
                sep = ""))
      pval <- rep(NA, dim(BestPos)[1])
      for(k in 1:(dim(BestPos)[1])) {
        if ((k %% 5000) == 0) {
          cat(paste("Step =  ", k, 
                    "   System.time = ", (Sys.time()), "\n", sep="")) 
        }
        Pos <- sort(BestPos[k, ])
        xgeno <- snps[, Pos]
        geno <- apply(xgeno, 1, paste, collapse = "")
        lr <- rep(lest, 1)
        fr <- rep(lest, 2^length(nloc) + 1 + lenx)
        hr <- rep(lest, lenx)
        tr <- rep(lest, 2^length(nloc) + 1 + lenx)
        pr <- rep(lest, 1)
        zzz <- .C("haptdpnZR", 
                  famid = as.character(famid),
                  pid = as.character(patid),
                  geno = as.character(geno), 
                  trait = as.integer(trait), 
                  lenx = as.integer(lenx),
                  lim = as.double(lim),
                  likeres = lr, 
                  freqres = fr, 
                  hapres = hr,
                  tdtres = tr, 
                  pvres = pr, 
                  package = "HapEstXXR")
        x <- strsplit(zzz$pvres, " ")
        x[[1]] <- x[[1]][x[[1]] != ""]
        x5 <- as.numeric(x[[1]][4])
        x1 <- as.numeric(x[[1]][1])
        if (x5 == 1){
          pv <- 1.0 - pchisq(x5, 1)
        } else {
          pv <- 1.0 - pchisq((x5 - 1) * x1 / x5, x5 - 1)
        }
        pval[k] <- pv
      } # end of for k
      ind  <- (order(pval))[1:min(nt, length(pval))]
      pval <- pval[ind]
      BestPos <- BestPos [ind, , drop = FALSE]
      res.list[[i]] <- data.frame(BestPos,
                                  rep("families", length(pval)), 
                                  pval,
                                  stringsAsFactors = FALSE)
      rownames(res.list[[i]]) <- NULL
      colnames(res.list[[i]]) <- c(paste("SNP", 1:i, sep = ""), 
                                   "type",
                                   "p.value")
      i <- i + 1
    }  # while
    return(result = res.list)
  }
