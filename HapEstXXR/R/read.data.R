read.data <-
  function(filename, linkage = TRUE, map = NA)
  {
    if(linkage == TRUE) {
      f <- unlist(strsplit(readLines(filename, n = 1),
                           split = "[ \t\n\f\r\v]+"))
      f <- f[f != ""]
      lf <- length(f)
      nmarker <-(lf - 6)
      if((nmarker %% 2) != 0)
        stop("Line 1 did not have expected structure.")
      tab <- read.table(filename, 
                        colClasses = c(rep("character", 5),
                                       "numeric",
                                       rep("character", nmarker)))
      geno <- allele2toR(as.matrix(tab[, 7:lf, drop = FALSE]))
      famid <- tab[, 1]
      patid <- tab[, 2]
      fid <- tab[, 3]
      mid <- tab[, 4]
      sex <- tab[, 5]
      pheno <- tab[, 6]
    } else  #(linkage == FALSE)  ##=> Klaus format
    {
      f <- unlist(strsplit(readLines(filename, n = 1),
                           split = "[ \t\n\f\r\v]+"))
      f <- f[f != ""]
      lf <- length(f)
      if(!(lf == 3 || lf == 4))
        stop("Line 1 did not have 3 or 4 elements.")
      if(lf == 3) {
        tab <- read.table(filename, 
                          colClasses = c("character", 
                                         "character", 
                                         "numeric"))
        if(ncol(tab) != 3)
          stop("incorrect file.")
        patid <- tab[, 1]
        geno <- matrix(as.numeric(unlist(strsplit(tab[, 2], 
                                                  split = ""))), 
                       nrow = length(patid), 
                       byrow = TRUE)
        pheno <- tab[, 3]
        famid <- patid
        nind <- length(pheno)
        fid <- rep(0, nind)
        mid <- fid
        sex <- rep(0, nind)
      }
      else {  # lf=4
        tab <- read.table(filename, colClasses = "character" )
        if(ncol(tab) != 4)
          stop("incorrect file.")
        famid <- tab[, 1]
        patid <- tab[, 2]
        N <- length(famid)
        famid.unique <- unique(famid)
        nf <- length(famid.unique)
        fid <- rep(0, N)
        mid <- rep(0, N)
        k <- 1
        for(i in famid.unique) {
          fami <- tab[tab[, 1] == i, , drop = FALSE]
          nfi <- dim(fami)[1]
          if(nfi > 2) {
            for(j in(k + 2):(k + nfi - 1)) {
              fid[j] <- patid[k]
              mid[j] <- patid[k + 1]
            }
          }
          k <- k + nfi
        }
        geno <- matrix(as.numeric(unlist(strsplit(tab[, 3], 
                                                  split = ""))), 
                       nrow = length(patid),
                       byrow = TRUE)
        pheno <- tab[, 4]
        sex <- rep(0, length(pheno))
      }
    }
    rm(tab)
    chr  <- pos <- snp <- NULL
    if(!is.na(map)) {
      inmap <- read.table(file = map, 
                          colClasses = c("numeric", 
                                         "character", "numeric"), 
                          stringsAsFactors = FALSE)
      chr <- inmap[, 1]
      snp <- inmap[, 2]
      pos <- inmap[, 4]
      rm(inmap)
    }
    rownames(geno) <- patid
    colnames(geno) <- snp
    return(list(famid = famid, 
                patid = patid, 
                fid = fid,
                mid = mid, 
                sex = sex, 
                genotypes = geno,
                trait = pheno,
                chr = chr, 
                snp = snp,
                pos = pos))
  }

