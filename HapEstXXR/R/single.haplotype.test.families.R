single.haplotype.test.families  <- 
  function(snps, trait, famid, patid, 
           fid, mid, lim = 0.05, sort = FALSE) {
    
    snps <- as.matrix(snps)
    trait <- as.numeric(trait)
    if (!all(trait[!is.na(trait)] ==  0 | trait[!is.na(trait)] ==  1))
      stop("trait should be 0 for unaffected or 1 for affcted children.")
    # remove unaffected children
    not.rm.child <- ifelse((fid != 0) & (mid!= 0) & (trait == 0), 
                           FALSE, TRUE)
    if (sum(!not.rm.child) > 0) {
      cat(sum(!not.rm.child), 
          " children removed because they were unaffected.\n", 
          sep = "")
    }    
    famid <- famid[not.rm.child]
    patid <- patid[not.rm.child]
    fid   <- fid[not.rm.child]
    mid   <- mid[not.rm.child]
    snps  <- snps[not.rm.child, , drop = FALSE]
    trait <- trait[not.rm.child] + 1
    # exclusion of nuclear families without two parents.  
    excl.fam <- NULL
    #i <- unique(famid)[1]
    for(i in unique(famid)) {
      selfam <- famid == i
      selfid <- unique(fid[selfam])
      selfid <- selfid[(!is.na(selfid)) & (selfid!= 0)]
      selmid <- unique(mid[selfam])
      selmid <- selmid[(!is.na(selmid)) & (selmid!= 0)]
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
        if (!any(patid[selfam] %in% selmid)) {
          excl.fam <- c(excl.fam, i)
        }
      }
    }
    excl.fam <- unique(excl.fam)
    if (!is.null(excl.fam)) {
      print(paste("Exclusion of nuclear families without two parents: ",
                  paste(excl.fam, collapse = " "), 
                  sep = ""))
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
    if (sort == TRUE) {
      ordfam <- order.families(famid, patid, fid, mid)
      snps  <- snps[ordfam, , drop = FALSE]
      trait <- trait[ordfam]
      famid <- famid[ordfam]
      patid <- patid[ordfam]
      fid   <- fid[ordfam]
      mid   <- mid[ordfam]
    }
    lenx <- length(famid)
    #from now on genotypes
    lpi <- dim(snps)[2]
    geno <- apply(snps, 1, paste, collapse = "")
    nloc  <- dim(snps)[2]
    lest <- paste(rep(" ", 10000), collapse = "")
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
    ### RESULT: haplotypes
    zzz$freqres <- zzz$freqres[zzz$freqres != lest]
    splits <- function(x, ind) {
      xs <- unlist(strsplit(c(x[ind[1]], x[ind[2]]), split = " "))
      xs <- xs[ xs !=  "" ]
      return(xs)
    }
    haps <- as.data.frame(lapply(strsplit(zzz$freqres, ">>"), 
                                 splits, ind = c(2, 4)), 
                          stringsAsFactors  = FALSE)
    hap  <- as.character(c(haps[1, ], "rare"))
    freq <- as.numeric(haps[2, ])
    freq <- c(freq, 1 - sum(freq))
    ### RESULT: p value
    #  pv <- as.numeric(unlist(strsplit(zzz$pvres, ">>")))
    #  names(pv) <- c("Statistic", "pvalue")
    x <- strsplit(zzz$pvres, " ")
    x[[1]] <- x[[1]][x[[1]] != ""]
    x5 <- as.numeric(x[[1]][4])
    x1 <- as.numeric(x[[1]][1])
    if (x5 == 1){
      pv <- 1.0 - pchisq(x5, 1)
    } else{
      pv <- 1.0 - pchisq((x5 - 1) * x1 / x5, x5 - 1)
    }
    pv <- data.frame(as.numeric(x1), as.integer(x5), as.numeric(pv))
    colnames(pv) <- c("Statistic", "nhap", "pvalue")
    ### RESULT: TDT - haplotype specific test
    tdt <- zzz$tdtres [zzz$tdtres!= lest]
    tdtres <- matrix(, ncol = 5, nrow = length(tdt))
    si <- min(which(unlist(strsplit(tdt[1], split = "")) == ":")) + 1
    for(i in 1:(length(tdt))) {
      #hapi.test <- unlist(strsplit(tdt[i], split = " "))[c(6, 8)]
      hapi.test <- as.numeric(c(substr(tdt[i], si, si + 7),
                                substr(tdt[i], si + 14, si + 21)))
      trans  <- hapi.test[2]
      ntrans <- hapi.test[1]
      chi2  <- (trans-ntrans)^2 / (trans + ntrans)
      pvali <- 1 - pchisq(chi2, 1)
      tdtres[i, ] <- c(trans, ntrans, trans / ntrans, chi2, pvali)
    }
    rownames(tdtres) <- c(hap)  
    colnames(tdtres) <- c("trans", "non-trans", 
                          "trans:non-trans", 
                          "chi square", "p value")
    ### order by haplotype frequencies
    ord <- order(freq, decreasing = TRUE)
    hap <- hap[ord]
    freq <- freq[ord]
    tdtres <- tdtres[ord, ]
    hapres <- data.frame(Hap = as.character(hap),
                         Freq = as.numeric(freq), 
                         stringsAsFactors = FALSE)
    res.list <- list(haplotypes = hapres, 
                     global.test = pv, 
                     haplotype.i = tdtres)
    return(res.list)
  }
