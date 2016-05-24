################################################################################
# bigPreDiv - a function for the calculation of diff stats from big files
################################################################################
bigPreDiv <- function(prePopList, bs = FALSE, nloci, npops, 
                      popSizes, fstat){
  ps <- popSizes
  if (bs) {
    popList <- lapply(prePopList, function(x) {
      boot <- sample(1:length(x[, 1, 1]), replace = TRUE)
      return(x[boot, (2:(nloci + 1)), ])
    })
  } else {
    popList <- lapply(prePopList, function(x) {
      return(x[, (2:(nloci + 1)), ])
    })
  }
  indtyp <- lapply(popList, function(x) {
    apply(x, 2, function(y) {
      length(na.omit(y[, 1]))
    })
  })
  alls <- lapply(seq_along(popList), function(i) {
    apply(popList[[i]], 2, function(x) {
      return(unique(c(x[, 1], x[, 2])))
    })
  })
  all_alleles <- lapply(1:nloci, function(i) {
    alleles <- lapply(alls, function(x) {
      return(x[[i]])
    })
    return(sort(unique(unlist(alleles))))
  })
  #   obsAlls <- lapply(popList, function(x) {
  #     apply(x, 2, function(y) {
  #       als <- unique(c(na.omit(y[, 1]), na.omit(y[, 2])))
  #       counts <- sapply(als, function(z) {
  #         res <- length(which(y == z))
  #         return(res)
  #       })
  #     })
  #   })
  obsAlls <- lapply(popList, function(x) {
    lapply(1:ncol(x), function(i){
      alls <- c(x[,i,1], x[,i,2])
      return(table(alls))
    })
  })
  allele_freq <- lapply(1:nloci, function(i) {
    loc <- matrix(nrow = length(all_alleles[[i]]), ncol = npops)
    rownames(loc) <- all_alleles[[i]]
    for (j in 1:npops) {
      o <- obsAlls[[j]][[i]]
      n <- indtyp[[j]][i]
      loc[names(o), j] <- o/(2 * n)
    }
    loc[is.na(loc)] <- 0
    return(loc)
  })
  preLoc <- lapply(indtyp, function(x) {
    return(1/x)
  })
  loci_harm_N <- sapply(1:nloci, function(i) {
    loc <- sapply(1:npops, function(j) {
      return(preLoc[[j]][i])
    })
    return(npops/sum(loc))
  })
  loci_harm_N <- round(loci_harm_N, 2)
  indtypLoc <- lapply(1:nloci, function(i) {
    res <- sapply(1:npops, function(j) {
      return(indtyp[[j]][i])
    })
  })
  rm(indtyp)
  if (fstat) {
    badData <- sapply(indtypLoc, function(y) {
      is.element(0, y)
    })
    if (sum(badData) > 0) {
      nl <- nloci - (sum(badData))
    } else {
      nl <- nloci
    }
    gdData <- which(!badData)
    badData <- which(badData)
    all_genot <- matrix(data = NA, nrow = sum(ps),
                        ncol =  length(gdData))
    for (i in 1:npops) {
      if (i == 1) {
        res <- apply(popList[[i]], 2, function(y) {
          return(paste0(y[, 1], y[, 2]))
        })
        all_genot[1:ps[i],] <- res[, gdData]
        rm(res)
      } else {
        res <- apply(popList[[i]], 2, function(y) {
          return(paste0(y[, 1], y[, 2]))
        })
        all_genot[(sum(ps[1:(i - 1)]) + 1):sum(ps[1:i]), ] <- res[, gdData]
        rm(res)
      }
    }
    all_genot[all_genot == "NANA"] <- NA
    genoCount <- lapply(1:ncol(all_genot), function(i){
      table(all_genot[,i])
    })
    nameFormat <- function(x) {
      nms <- names(x)
      lgth <- nchar(nms[1])
      newNms <- sapply(nms, function(y) {
        paste(substr(y, 1, lgth/2), "/", substr(y, (lgth/2) + 
                                                  1, lgth), sep = "")
      })
      names(x) <- newNms
      return(x)
    }
    genoCount <- lapply(genoCount, nameFormat)
    h_sum <- list()
    for (i in 1:length(gdData)) {
      h_sum[[i]] <- vector()
      cnSplit <- strsplit(names(genoCount[[i]]), "/")
      for (j in 1:length(all_alleles[[gdData[i]]])) {
        het_id1 <- lapply(cnSplit, is.element, all_alleles[[gdData[i]]][j])
        het_id2 <- lapply(het_id1, sum)
        het_id1 <- which(het_id2 == 1)
        h_sum[[i]][j] <- sum(genoCount[[i]][het_id1])
      }
    }
    indtyp_tot <- lapply(indtypLoc, sum)
    kk_hsum <- lapply(1:ncol(all_genot), function(i) {
      list(h_sum[[i]], indtyp_tot[[gdData[i]]])
    })
    kk_hbar <- lapply(kk_hsum, function(x) {
      return(x[[1]]/x[[2]])
    })
    pdat <- lapply(1:length(all_genot[1, ]), function(i) {
      list(allele_freq[[gdData[i]]], indtypLoc[[gdData[i]]])
    })
    kk_p <- lapply(pdat, function(x) {
      if (is.null(x[[1]]) == FALSE) {
        apply(x[[1]], 1, function(y) {
          y * (2 * x[[2]])
        })
      }
    })
    res <- matrix(0, (nloci + 1), 2)
    colnames(res) <- c("Fst_WC", "Fit_WC")
    A <- vector()
    a <- vector()
    b <- vector()
    c <- vector()
    for (i in 1:length(gdData)) {
      kknbar <- indtyp_tot[[gdData[i]]]/npops
      kknC <- (indtyp_tot[[gdData[i]]] - sum(indtypLoc[[gdData[i]]]^2)/indtyp_tot[[gdData[i]]])/(npops - 1)
      kkptild <- kk_p[[i]]/(2 * indtypLoc[[gdData[i]]])
      kkptild[kkptild == "NaN"] <- NA
      kkpbar <- colSums(kk_p[[i]])/(2 * indtyp_tot[[gdData[i]]])
      kks2 <- colSums(indtypLoc[[gdData[i]]] * (kkptild - rep(kkpbar, each = npops))^2)/((npops - 1) * kknbar)
      kkA <- kkpbar * (1 - kkpbar) - (npops - 1) * kks2/npops
      kka <- kknbar * (kks2 - (kkA - (kk_hbar[[i]]/4))/(kknbar - 1))/kknC
      kkb <- (kknbar/(kknbar - 1))*(kkA-((2*kknbar-1)/(4*kknbar))*kk_hbar[[i]])
      #kkb <- kknbar * (kkA - (2 * (kknbar - 1)) * kk_hbar[[i]]/(4 * kknbar))/(kknbar - 1)
      kkc <- kk_hbar[[i]]/2
      A[i] <- sum(kkA, na.rm = TRUE)
      a[i] <- sum(kka, na.rm = TRUE)
      b[i] <- sum(kkb, na.rm = TRUE)
      c[i] <- sum(kkc, na.rm = TRUE)
      res[gdData[i], "Fst_WC"] <- round(sum(kka)/sum(kka + kkb + kkc), 4)
      res[gdData[i], "Fit_WC"] <- round(1 - sum(kkc)/sum(kka + kkb + kkc), 4)
    }
    res[res == "NaN"] <- NA
    res[res == 0] <- NA
    sumA <- sum(A, na.rm = TRUE)
    suma <- sum(a, na.rm = TRUE)
    sumb <- sum(b, na.rm = TRUE)
    sumc <- sum(c, na.rm = TRUE)
    res[(nloci + 1), "Fst_WC"] <- round(suma/(suma + sumb + sumc), 4)
    res[(nloci + 1), "Fit_WC"] <- round(1 - sumc/(suma + sumb + sumc), 4)
    z <- gc(reset = TRUE)
    rm(z)
    fst <- res
    rm(res)
  }
  ho <- lapply(popList, function(x) {
    apply(x, 2, function(y) {
      1 - (sum(na.omit(y[, 1] == y[, 2]))/length(na.omit(y[, 1])))
    })
  })
  he <- t(sapply(allele_freq, function(x) {
    apply(x, 2, function(y) {
      return(1 - sum(y^2))
    })
  }))
  mf <- lapply(allele_freq, function(x) {
    rowSums(x)/ncol(x)
  })
  ht <- sapply(mf, function(x) {
    1 - sum(x^2)
  })
  hs <- rowSums(he)/npops
  hs_est <- hs * ((2 * loci_harm_N)/((2 * loci_harm_N) - 1))
  ht_est <- ht + (hs_est/(2 * loci_harm_N * npops))
  ht_est[is.nan(ht_est)] <- NA
  hst <- round((ht - hs)/(1 - hs), 4)
  dst <- round(ht - hs, 4)
  gst <- round(dst/ht, 4)
  gst[is.nan(gst)] <- NA
  djost <- round((dst/(1 - hs)) * (npops/(npops - 1)), 4)
  djost[djost == 0] <- NA
  hst_est <- round((ht_est - hs_est)/(1 - hs_est), 4)
  dst_est <- round(ht_est - hs_est, 4)
  gst_est <- round(dst_est/ht_est, 4)
  gst_est[is.nan(gst_est)] <- NA
  gst_max <- ((npops - 1) * (1 - hs))/(npops - 1 + hs)
  gst_est_max <- (((npops - 1) * (1 - hs_est))/(npops - 1 + 
                                                  hs_est))
  gst_hedrick <- round(gst/gst_max, 4)
  gst_est_hedrick <- round(gst_est/gst_est_max, 4)
  gst_est_hedrick[gst_est_hedrick > 1] <- 1
  djost_est <- round((npops/(npops - 1)) * ((ht_est - hs_est)/(1 - hs_est)), 4)
  djost_est[djost_est == 0] <- NA
  ht_mean <- round(mean(ht, na.rm = TRUE), 4)
  hs_mean <- round(mean(hs), 4)
  gst_all <- round((ht_mean - hs_mean)/ht_mean, 4)
  gst_all_max <- round(((npops - 1) * (1 - hs_mean))/(npops - 
                                                        1 + hs_mean), 4)
  gst_all_hedrick <- round(gst_all/gst_all_max, 4)
  djost_all <- round(((ht_mean - hs_mean)/(1 - hs_mean)) * 
                       (npops/(npops - 1)), 4)
  hs_est_mean <- mean(hs_est, na.rm = TRUE)
  ht_est_mean <- mean(ht_est, na.rm = TRUE)
  gst_est_all <- round((ht_est_mean - hs_est_mean)/ht_est_mean, 
                       4)
  gst_est_all_max <- round((((npops - 1) * (1 - hs_est_mean))/(npops - 
                                                                 1 + hs_est_mean)), 4)
  gst_est_all_hedrick <- round(gst_est_all/gst_est_all_max, 
                               4)
  gst_est_all_hedrick[gst_est_all_hedrick > 1] <- 1
  if (nloci == 1) {
    djost_est_all <- round(djost_est, 4)
  } else {
    djost_est_all <- round(1/(1/mean(djost_est, na.rm = TRUE) + 
                                (var(djost_est, na.rm = TRUE) * (1/mean(djost_est, 
                                                                        na.rm = TRUE))^3)), 4)
  }
  djost_est[djost_est == 0] <- NaN
  djost[djost == 0] <- NaN
  if (fstat) {
    list(hst = hst, dst = dst, gst = gst, gst_hedrick = gst_hedrick, 
         djost = djost, locus_harmonic_N = loci_harm_N, hst_est = hst_est, 
         dst_est = dst_est, gst_est = gst_est, 
         gst_est_hedrick = gst_est_hedrick, 
         djost_est = djost_est, gst_all = gst_all, 
         gst_all_hedrick = gst_all_hedrick, 
         djost_all = djost_all, gst_est_all = gst_est_all, 
         gst_est_all_hedrick = gst_est_all_hedrick, 
         djost_est_all = djost_est_all, 
         fstats = fst)
  } else {
    list(hst = hst, dst = dst, gst = gst, gst_hedrick = gst_hedrick, 
         djost = djost, locus_harmonic_N = loci_harm_N, hst_est = hst_est, 
         dst_est = dst_est, gst_est = gst_est, 
         gst_est_hedrick = gst_est_hedrick, 
         djost_est = djost_est, gst_all = gst_all, 
         gst_all_hedrick = gst_all_hedrick, 
         djost_all = djost_all, gst_est_all = gst_est_all, 
         gst_est_all_hedrick = gst_est_all_hedrick, 
         djost_est_all = djost_est_all)
  }
}
################################################################################
# END - bigPreDiv
################################################################################