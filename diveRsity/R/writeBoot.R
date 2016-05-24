#' writeBoot: A function to write unprocessed bootstrapped matrices of 
#' differentiation and diversity partition statistics to file.
#' 
#' The function is based on fastDivPart from diveRsity.
#' 
#' Kevin Keenan 2014
#' 
#' 
writeBoot <- function(infile = NULL, outfile = NULL, gp = 3, bootstraps = 0, 
                      parallel = FALSE){
  ############################ Argument definitions ############################
  # define arguments for testing
  D <- infile
  on <- outfile
  fst <- TRUE
  bstrps <- bootstraps
  bsls <- FALSE
  bspw <- TRUE
  plt <- FALSE
  para <- parallel
  pWise <- TRUE
  
  ##############################################################################
  #Use pre.div to calculate the standard global and locus stats
  accDat <- pre.divLowMemory(list(infile = D,
                                  gp = gp,
                                  bootstrap = FALSE,
                                  locs = TRUE,
                                  fst = fst,
                                  min = FALSE))
  # create a directory for output
  if(!is.null(on)){
    suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                           "-[writeBoot]","/",sep="")))
  }
  of = paste(getwd(), "/", on, "-[writeBoot]", "/", sep = "")
  
  
  if(para){
    para_pack_inst<-is.element(c("parallel","doParallel","foreach","iterators"),
                               installed.packages()[,1])
  }
  
  para_pack <- all(para_pack_inst)
  ############################################################################
  ################################## Pairwise ################################
  ############################################################################
  # population pair combinations
  
  # define new functions
  ############################################################################
  ############################################################################
  # pwCalc
  ############################################################################
  # New optimised function for the calculation of pairwise statistics
  # Returns a 3D array where each 'slot' represents the pairwise matrix
  # for Gst_est, G_st_est_hed and D_jost_est respectively
  
  # Kevin Keenan
  # 2013
  
  pwCalc <- function(infile, fst,  bs = FALSE){
    
    
    #   # uncomment for testing
    #   infile <- "pw_test.txt"
    #   source("readGenepopX.R")
    #   # read pwBasicCalc function
    #   source("pwBasicCalc.R")
    # define baseline info
    dat <- readGenepopX(list(infile = infile,
                             bootstrap = bs))
    if(fst){
      # calculate all fst
      fstat <- pwFstWC(dat)
      # extract locus theta and variance components
      locTheta <- lapply(fstat, "[[", 1)
      # sum res
      aLoc <- Reduce(`+`, lapply(fstat, "[[", 2))
      bLoc <- Reduce(`+`, lapply(fstat, "[[", 3))
      cLoc <- Reduce(`+`, lapply(fstat, "[[", 4))
      # calculate pw Fst across loci
      pwTheta <- aLoc/(aLoc+bLoc+cLoc)
      # clean up
      rm(aLoc, bLoc, cLoc, fstat)
      z <- gc()
      rm(z)
    }
    # extract allele frequencies
    af <- dat$allele_freq
    # extract harmonic mean sample sizes
    
    # make space in RAM
    dat$allele_freq <- NULL
    z <- gc()
    rm(z)
    # extract npops and nloci
    npops <- dat$npops
    nloci <- dat$nloci
    # define pairwise relationships
    pw <- combn(dat$npops, 2)
    # generate pairwise locus harmonic mean sample sizes
    indtyp <- dat$indtyp
    pwHarm <- lapply(indtyp, pwHarmonic, pw = pw)
    
    
    # calculate pairwise ht and hs
    hths <- mapply(pwBasicCalc, af, pwHarm,
                   MoreArgs = list(pw = pw, npops = dat$npops),
                   SIMPLIFY = FALSE)
    # seperate ht and hs
    #   htLoc <- lapply(hths, "[[", 1)
    #   hsLoc <- lapply(hths, "[[", 2)
    # seperate ht_est and hs_est
    hsEstLoc <- lapply(hths, "[[", 1)
    htEstLoc <- lapply(hths, "[[", 2)
    
    # clean up
    rm(hths)
    z <- gc()
    rm(z)
    
    # Calculate locus stats
    # Standard locus stats
    # locus Gst
    #   gstLoc <- mapply(FUN = gstCalc, ht = htLoc, hs = hsLoc, 
    #                    SIMPLIFY = FALSE)
    #   # locus G'st
    #   gstHedLoc <- mapply(FUN = gstHedCalc, ht = htLoc, hs = hsLoc,
    #                       SIMPLIFY = FALSE)
    #   # locus D_jost
    #   dLoc <- mapply(FUN = djostCalc, ht = htLoc, hs = hsLoc,
    #                  SIMPLIFY = FALSE)
    
    # Estimated locus stats
    # locus Gst_est
    gstLocEst <- mapply(FUN = gstCalc, ht = htEstLoc, 
                        hs = hsEstLoc, 
                        SIMPLIFY = FALSE)
    # locus G'st_est
    gstHedLocEst <- mapply(FUN = gstHedCalc, ht = htEstLoc, 
                           hs = hsEstLoc,
                           SIMPLIFY = FALSE)
    # locus D_jost_est
    dLocEst <- mapply(FUN = djostCalc, ht = htEstLoc, 
                      hs = hsEstLoc,
                      SIMPLIFY = FALSE)
    
    #   # calculate mean ht and hs
    #   htMean <- Reduce(`+`, htLoc)/nloci
    #   hsMean <- Reduce(`+`, hsLoc)/nloci
    # calculate mean ht_est and hs_est
    htEstMean <- Reduce(`+`, htEstLoc)/nloci
    hsEstMean <- Reduce(`+`, hsEstLoc)/nloci
    
    # calculate standard stats (uncomment for loc stats)
    
    #   # overall dst
    #   dstAll <- htMean - hsMean
    #   # overall gst (Nei 1973)
    #   gstAll <- (dstAll)/htMean
    #   # overall max gst (Hedricks 2005)
    #   gstAllMax <- ((2 - 1)*(1 - hsMean)) / ((2 - 1) + hsMean)
    #   # overall Hedricks' Gst
    #   gstAllHedrick <- gstAll/gstAllMax
    #   # Overall D_jost (Jost 2008)
    #   djostAll <- (dstAll/(1-hsMean))*(2/(2-1))
    
    # Calculate estimated stats
    
    # Overall estimated dst
    dstEstAll <- htEstMean - hsEstMean
    # Overall estimated Gst (Nei & Chesser, 1983)
    gstEstAll <- dstEstAll/htEstMean
    # Overall estimated max Gst (Hedricks 2005)
    gstEstAllMax <- ((2-1)*(1-hsEstMean))/(2-1+hsEstMean)
    # Overall estimated Hedricks' Gst
    gstEstAllHed <- gstEstAll/gstEstAllMax
    # Overall estimated D_Jost (Chao et al., 2008)
    if(nloci == 1){
      djostEstAll <- (2/(2-1))*((dstEstAll)/(1 - hsEstMean))
    } else {
      dLocEstMn <- Reduce(`+`, dLocEst)/nloci
      # calculate variance (convert dLocEst to an array)
      dLocEst1 <- array(unlist(dLocEst), 
                        dim = c(nrow(dLocEst[[1]]), 
                                ncol(dLocEst[[1]]), 
                                length(dLocEst)))
      dLocEstVar <- apply(dLocEst1, c(1,2), var)
      djostEstAll <- 1/((1/dLocEstMn)+((dLocEstVar*((1/dLocEstMn)^3))))
      # tidy up
      rm(dLocEstMn, dLocEstVar)
      z <- gc()
      rm(z)
    }
    
    # define a function to arrange locus stats into arrays
    #   arrDef <- function(x){
    #     return(array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x))))
    #   }
    if(fst){
      resArr <- array(c(gstEstAll, gstEstAllHed, djostEstAll, pwTheta),
                      dim = c(nrow(gstEstAll),
                              ncol(gstEstAll),
                              4))
      lstats <- array(NA, dim = c(dat$npops, dat$npops, dat$nloci, 4))
      lstats[,,,1] <- unlist(gstLocEst)
      lstats[,,,2] <- unlist(gstHedLocEst)
      lstats[,,,3] <- unlist(dLocEst)
      lstats[,,,4] <- unlist(locTheta)
      #lstats <- array(list(gstLocEst, gstHedLocEst, dLocEst, locTheta)
      #lstats <- mapply(FUN = `list`, gstLocEst, gstHedLocEst, dLocEst, locTheta,
      #                 SIMPLIFY=FALSE)
    } else {
      resArr <- array(c(gstEstAll, gstEstAllHed, djostEstAll),
                      dim = c(nrow(gstEstAll),
                              ncol(gstEstAll), 3))
      lstats <- array(NA, dim = c(dat$npops, dat$npops, dat$nloci, 3))
      lstats[,,,1] <- unlist(gstLocEst)
      lstats[,,,2] <- unlist(gstHedLocEst)
      lstats[,,,3] <- unlist(dLocEst)
      #lstats <- list(gstLocEst, gstHedLocEst, dLocEst)
      #lstats <- mapply(FUN = `list`, gstLocEst, gstHedLocEst, dLocEst,
      #                 SIMPLIFY = FALSE)
    }
    # arrange loci into arrays
    #locOut <- lapply(lstats, arrDef)
    
    
    list(resArr = resArr,
         locOut = lstats)
  }
  ############################################################################
  # END - pwDivCalc
  ############################################################################
  # Calculate Weir & Cockerham's F-statistics (optimised)
  ############################################################################
  # pwFstWC: a function co calculate weir and cockerhams fis, fit, and fst
  ############################################################################
  pwFstWC<-function(rdat){
    #   rdat <- diveRsity::readGenepop("KK_test1v2.gen")
    pw <- combn(rdat$npops, 2)
    #   # account for loci with missing info for pops
    #   pwBadData <- function(indtyp, pw){
    #     out <- sapply(1:ncol(pw), function(i){
    #       is.element(0, indtyp[pw[,i]])
    #     })
    #   }
    #   badDat <- sapply(rdat$indtyp, pwBadData, pw = pw)
    #   if(any(badDat)){
    #     bd <- TRUE
    #   }
    #   # determine the number of loci per pw comparison
    #   nlocPw <- apply(badDat, 1, function(x){
    #     if(sum(x) > 0){
    #       nl <- rdat$nloci - sum(x)
    #     } else {
    #       nl <- rdat$nloci
    #     }
    #   })
    #   # define all good data
    #   gdDat <- lapply(1:nrow(badDat), function(i){
    #     which(!badDat[i,])
    #   })
    #   badDat <- lapply(1:nrow(badDat), function(i){
    #     which(badDat[i,])
    #   })
    # get all genotypes for each pw comparison
    allGenot <- apply(pw, 2, function(x){
      list(rdat$pop_list[[x[1]]], 
           rdat$pop_list[[x[2]]])
    })
    #   # filter bad data
    #   if(any(nlocPw != rdat$nloci)){
    #     idx <- which(nlocPw != rdat$nloci)
    #     for(i in idx){
    #       allGenot[[i]][[1]] <- allGenot[[i]][[1]][, gdDat[[i]]]
    #       allGenot[[i]][[2]] <- allGenot[[i]][[2]][, gdDat[[i]]]
    #     }
    #   }
    # unlist pw genotype data
    allGenot <- lapply(allGenot, function(x){
      return(do.call("rbind", x))
    })
    # identify unique genotypes
    #   genot <- lapply(allGenot, function(x){
    #     return(apply(x, 2, function(y){
    #       unique(na.omit(y))
    #     }))
    #   })
    # count number of genotypes per pw per loc
    
    genoCount <- lapply(allGenot, function(x){
      if(NCOL(x) == 1){
        return(list(table(x)))
      } else {
        lapply(1:ncol(x), function(i) table(x[,i]))
      }
    })
    
    
    #   genoCount <- lapply(allGenot, function(x){
    #     lapply(split(x,seq(NCOL(x))),table) # accounts for single loci
    #     #apply(x, 2, table)
    #   })
    
    
    # function to count heterozygotes
    htCount <- function(x){
      nms <- names(x)
      ncharGeno <- nchar(nms[1])
      alls <- cbind(substr(nms, 1, (ncharGeno/2)),
                    substr(nms, ((ncharGeno/2) + 1), ncharGeno))
      unqAlls <- unique(as.vector(alls))
      hetCounts <- sapply(unqAlls, function(a){
        idx <- which(rowSums(alls == a) == 1)
        return(sum(x[idx]))
      })
      return(hetCounts)
    }
    # hSum is the total observed hets per allele
    hSum <- lapply(genoCount, function(x){
      out <- lapply(x, htCount)
    })
    
    #   if(bd){
    #     # insert na for missing loci
    #     hSum <- lapply(seq_along(badDat), function(i){
    #       naPos <- badDat[[i]]
    #       idx <- c(seq_along(hSum[[i]]), (naPos - 0.5))
    #       return(c(hSum[[i]], rep(NA, length(naPos)))[order(idx)])
    #     }) 
    #   }
    # convert to locus orientated hSum
    hSum <- lapply(seq_along(hSum[[1]]), function(i){
      lapply(hSum, "[[", i)
    })
    
    # total ind typed per loc per pw
    indTypTot <- lapply(rdat$indtyp, function(x){
      return(apply(pw, 2, function(y){
        sum(x[y])
      }))
    })
    # nBar is the mean number of inds per pop
    nBar <- lapply(indTypTot, `/`, 2)
    
    # hbar per pw per loc
    hBar <- lapply(seq_along(hSum), function(i){
      divd <- indTypTot[[i]]
      return(mapply(`/`, hSum[[i]], divd, SIMPLIFY = FALSE))
    })
    
    # p per loc per pw
    pCalc <- function(x, y, pw){
      out <- lapply(seq_along(pw[1,]), function(i){
        return(cbind((x[,pw[1,i]]*(2*y[pw[1,i]])),
                     (x[,pw[2,i]]*(2*y[pw[2,i]]))))
      })
      return(out)
    }
    p <- mapply(FUN = pCalc, x = rdat$allele_freq, 
                y = rdat$indtyp, 
                MoreArgs = list(pw = pw), 
                SIMPLIFY = FALSE)
    
    #   # convert p elements into array structure
    #   pArr <- lapply(p, function(x){
    #     d3 <- length(x)
    #     d2 <- 2
    #     d1 <- nrow(x[[1]])
    #     return(array(unlist(x), dim = c(d1, d2, d3)))
    #   })
    
    fstatCal <- function(indT, indtyp, hBar, nBar, p, pw, npops){
      #         indT=indTypTot[[1]]
      #         indtyp=rdat$indtyp[[1]]
      #         hBar <- hBar[[1]]
      #         nBar <- nBar[[1]]
      #         p <- p[[1]]
      #         pw <- pw
      #         npops <- rdat$npops
      indLocPwSqSum <- sapply(seq_along(pw[1,]), function(i){
        return(sum(indtyp[pw[,i]]^2))
      })
      indtypPw <- lapply(1:ncol(pw), function(idx){
        return(indtyp[pw[,idx]])
      })
      nC <- indT - (indLocPwSqSum/indT)
      ptildCalc <- function(x,y){ 
        return(cbind((x[,1]/(2*y[1])),
                     (x[,2]/(2*y[2]))))
      }
      pTild <- mapply(FUN = ptildCalc, x = p, y = indtypPw,
                      SIMPLIFY = FALSE)
      pBar <- lapply(seq_along(p), function(i){
        return(rowSums((p[[i]])/(2*indT[i])))
      })
      s2 <- lapply(seq_along(pBar), function(i){
        pp <- (pTild[[i]]-pBar[[i]])^2
        pp <- cbind((pp[,1]*indtypPw[[i]][1]),
                    (pp[,2]*indtypPw[[i]][2]))
        pp <- rowSums(pp)
        return((pp/(1*nBar[i])))
      })
      A <- lapply(seq_along(pBar), function(i){
        return(pBar[[i]]*(1-pBar[[i]])-(1)*s2[[i]]/2)
      })
      # fix hBar for unequal lengths
      idx <- lapply(seq_along(A), function(i){
        out <- match(names(A[[i]]), names(hBar[[i]]))
        return(which(!is.na(out)))
      })
      A <- lapply(seq_along(A), function(i){
        return(A[[i]][idx[[i]]])
      })
      s2 <- lapply(seq_along(s2), function(i){
        return(s2[[i]][idx[[i]]])
      })
      a <- lapply(seq_along(s2), function(i){
        return(nBar[[i]]*(s2[[i]]-(A[[i]]-(hBar[[i]]/4))/(nBar[[i]]-1))/nC[[i]])
      })
      #     a <- lapply(seq_along(s2), function(i){
      #       return(a[[i]][idx[[i]]])
      #     })
      b <- lapply(seq_along(A), function(i){
        return((nBar[[i]]/(nBar[[i]]-1))*(A[[i]]-((2*nBar[[i]]-1)/(4*nBar[[i]]))*hBar[[i]]))
        #return((nBar[[i]]/(nBar[[i]]-1))*(A[[i]]-(2*(nBar[[i]]-1))*hBar[[i]]/(4*nBar[[i]])))
      })
      #     b <- lapply(seq_along(A), function(i){
      #       return(b[[i]][idx[[i]]])
      #     })
      cdat <- lapply(seq_along(A), function(i){
        return(hBar[[i]]/2)
      })
      #     cdat <- lapply(seq_along(A), function(i){
      #       return(cdat[[i]][idx[[i]]])
      #     })
      A <- sapply(A, sum)
      a <- sapply(a, sum)
      b <- sapply(b, sum)
      cdat <- sapply(cdat, sum)
      theta <- a/(a+b+cdat)
      pwMat <- matrix(ncol = npops, nrow = npops)
      aMat <- matrix(ncol = npops, nrow = npops)
      bMat <- matrix(ncol = npops, nrow = npops)
      cMat <- matrix(ncol = npops, nrow = npops)
      for(i in 1:ncol(pw)){
        pwMat[pw[2,i], pw[1,i]] <- theta[i]
        aMat[pw[2,i], pw[1,i]] <- a[i]
        bMat[pw[2,i], pw[1,i]] <- b[i]
        cMat[pw[2,i], pw[1,i]] <- cdat[i]
      }
      pwMat[is.nan(pwMat)] <- NA
      aMat[is.nan(aMat)] <- NA
      cMat[is.nan(bMat)] <- NA
      bMat[is.nan(bMat)] <- NA
      
      list(pwMat, aMat, bMat, cMat)
    }
    
    # run fstatCal for each locus
    pwLoc <- mapply(FUN = fstatCal, indT = indTypTot,
                    indtyp = rdat$indtyp, hBar = hBar,
                    nBar = nBar, p = p, 
                    MoreArgs = list(pw = pw, npops = rdat$npops),
                    SIMPLIFY = FALSE)
    return(pwLoc)
  }
  ############################################################################
  # END - pwDivCalc
  ############################################################################
  # pwBasicCalc: a small function for calculating pairwise ht and hs 
  ############################################################################
  pwBasicCalc <- function(af, sHarm, pw, npops){
    ht <- matrix(ncol = npops, nrow = npops)
    hs <- matrix(ncol = npops, nrow = npops)
    htEst <- matrix(ncol = npops, nrow = npops)
    hsEst <- matrix(ncol = npops, nrow = npops)
    for(i in 1:ncol(pw)){
      id1 <- pw[1,i]
      id2 <- pw[2,i]
      # locus ht
      ht[id2, id1] <- 1 - sum(((af[,id1] + af[,id2])/2)^2)
      # locus hs
      hs[id2, id1] <- 1 - sum((af[,id1]^2 + af[,id2]^2)/2)
      # locus hs_est
      hsEst[id2, id1] <- hs[id2, id1]*((2*sHarm[id2,id1])/(2*sHarm[id2,id1]-1))
      # locus ht_est
      htEst[id2, id1] <- ht[id2, id1] + (hsEst[id2, id1]/(4*sHarm[id2, id1]))
    }
    #   ht[is.nan(ht)] <- 0
    #   hs[is.nan(hs)] <- 0
    htEst[is.nan(htEst)] <- 0
    hsEst[is.nan(hsEst)] <- 0
    list(hsEst = hsEst,
         htEst = htEst)
  }
  ############################################################################
  # END - pwBasicCalc
  ############################################################################
  
  # define locus stat calculators
  gstCalc <- function(ht, hs){
    return((ht - hs)/ht)
  }
  
  gstHedCalc <- function(ht, hs){
    gstMax <- ((2-1)*(1-hs))/(2-1+hs)
    return(((ht-hs)/ht)/gstMax)
  }
  
  djostCalc <- function(ht, hs){
    return((2/1)*((ht-hs)/(1-hs)))
  }
  
  # calculate pairwise locus harmonic mean
  pwHarmonic <- function(lss, pw){
    np <- length(lss)
    lhrm <- matrix(ncol = np, nrow = np)
    pwSS <- cbind(lss[pw[1,]], lss[pw[2,]])
    lhrmEle <- (0.5 * ((pwSS[,1]^-1) + (pwSS[,2]^-1)))^-1
    for(i in 1:ncol(pw)){
      idx1 <- pw[1,i]
      idx2 <- pw[2,i]
      lhrm[idx2, idx1] <- lhrmEle[i]
    }
    return(lhrm)
  }
  ############################################################################
  # pwDivCalc: a small function for calculating pairwise ht and hs 
  ############################################################################
  pwDivCalc <- function(x, pw, npops){
    ht <- matrix(ncol = npops, nrow = npops)
    hs <- matrix(ncol = npops, nrow = npops)
    for(i in 1:ncol(pw)){
      gamma <- sum(sqrt(abs(x[,pw[1,i]] * x[,pw[2,i]])))^-1 
      f <- gamma * sqrt(x[,pw[1,i]] * x[,pw[2,i]])
      ht[pw[1,i],pw[2,i]] <- 1 - sum(((f + x[,pw[1,i]])/2)^2)
      ht[pw[2,i],pw[1,i]] <- 1 - sum(((f + x[,pw[2,i]])/2)^2)
      hs[pw[1,i],pw[2,i]] <- 1 - sum((f^2 + x[,pw[1,i]]^2)/2)
      hs[pw[2,i],pw[1,i]] <- 1 - sum((f^2 + x[,pw[2,i]]^2)/2)
    }
    ht[is.nan(ht)] <- 0
    hs[is.nan(hs)] <- 0
    list(ht = ht, 
         hs = hs)
  }
  ############################################################################
  # END - pwDivCalc
  ############################################################################
  
  ############################################################################
  ############################################################################
  # working well 24/10/13
  if(pWise || bspw){
    # get pw names
    pw <- combn(accDat$npops, 2)
    popNms <- accDat$pop_names
    # for pw bootstrap table
    pw_nms <- paste(popNms[pw[1,]], popNms[pw[2,]], sep = " vs. ")
    
    pwStats <- pwCalc(D, fst, bs = FALSE)
    # extract stats
    gstPW <- pwStats$resArr[,,1]
    gstHPW <- pwStats$resArr[,,2]
    dPW <- pwStats$resArr[,,3]
    if(fst){
      thetaPW <- pwStats$resArr[,,4]
    }
    # clean up
    locstats <- pwStats$locOut
    rm(pwStats)
    z <- gc()
    rm(z)
    #       spc1 <- rep("", ncol(gstPW))
    #       if(fst){
    #         statNms <- c("Gst_est", "G'st_est", "Djost_est", "Fst_WC")
    #         outobj <- rbind(c(statNms[1], spc1), 
    #                         c("", popNms),
    #                         cbind(popNms, round(gstPW, 4)),
    #                         c(statNms[2], spc1),
    #                         c("", popNms),
    #                         cbind(popNms, round(gstHPW, 4)), 
    #                         c(statNms[3], spc1),
    #                         c("", popNms),
    #                         cbind(popNms, round(dPW, 4)), 
    #                         c(statNms[4], spc1),
    #                         c("", popNms),
    #                         cbind(popNms, round(thetaPW, 4)))
    #         outobj[is.na(outobj)] <- ""
    #         pwMatListOut <- list(gstPW, gstHPW, dPW, thetaPW)
    #         # add names to pwMatListOut
    #         names(pwMatListOut) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
    #         # tidy up
    #         rm(gstPW, gstHPW, dPW, thetaPW)
    #         z <- gc()
    #         rm(z)
    #       } else {
    #         statNms <- c("Gst_est", "G'st_est", "Djost_est")
    #         outobj <- rbind(c(statNms[1], spc1), 
    #                         c("", popNms),
    #                         cbind(popNms, round(gstPW, 4)),
    #                         c(statNms[2], spc1),
    #                         c("", popNms),
    #                         cbind(popNms, round(gstHPW, 4)), 
    #                         c(statNms[3], spc1),
    #                         c("", popNms),
    #                         cbind(popNms, round(dPW, 4)))
    #         outobj[is.na(outobj)] <- ""
    #         pwMatListOut <- list(gstPW, gstHPW, dPW)
    #         # add names to pwMatListOut
    #         names(pwMatListOut) <- c("gstEst", "gstEstHed", "djostEst")
    #         # tidy up
    #         rm(gstPW, gstHPW, dPW)
    #         z <- gc()
    #         rm(z)
    #       }
    # 
    #       for(i in 1:length(pwMatListOut)){
    #         dimnames(pwMatListOut[[i]]) <- list(popNms, popNms)
    #       }
    
    # convert locstats in to list format
    #       locstats <- lapply(apply(locstats, 4, list), function(x){
    #         lapply(apply(x[[1]], 3, list), function(y){
    #           out <- y[[1]]
    #           dimnames(out) <- list(popNms, popNms)
    #           return(out)
    #         })
    #       })
    # add names etc
    #       if(fst){
    #         # prepare locstats for output
    #         names(locstats) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
    #       } else {
    #         names(locstats) <- c("gstEst", "gstEstHed", "djostEst")
    #       }
    #       for(i in 1:length(locstats)){
    #         names(locstats[[i]]) <- accDat$locus_names
    #       }
    
  }
  
  #Bootstrap
  if(bspw == TRUE){
    if (para && para_pack) {
      
      cl <- parallel::makeCluster(detectCores())
      parallel::clusterExport(cl, c("pwCalc", "fst", "D", "readGenepopX",
                                    "fileReader", "pwFstWC", "pwHarmonic",
                                    "pwBasicCalc", "djostCalc", "gstCalc",
                                    "gstHedCalc"), 
                              envir = environment())
      pwBsStat <- parallel::parLapply(cl, 1:bstrps, function(...){
        return(pwCalc(infile = D, fst, bs = TRUE))
      })
      parallel::stopCluster(cl)
    } else {
      pwBsStat <- lapply(1:bstrps, function(...){
        return(pwCalc(D, fst, bs = TRUE))
      })
    }
    
    
    # seperate each stat
    
    gstEst <- lapply(pwBsStat, function(x){
      x$resArr[,,1]
    })
    
    gstEstHed <- lapply(pwBsStat, function(x){
      x$resArr[,,2]
    })
    
    dEst <- lapply(pwBsStat, function(x){
      x$resArr[,,3]
    })
    
    if(fst){
      theta <- lapply(pwBsStat, function(x){
        x$resArr[,,4]
      })
    }
    #pwBsLoc <- lapply(pwBsStat, "[[", 2)
    # tidy up
    rm(pwBsStat)
    z <- gc()
    rm(z)
    
    # convert bs lists to arrays for calculations
    if(fst){
      stats <- list(gstEst = array(unlist(gstEst),
                                   dim = c(nrow(gstEst[[1]]),
                                           nrow(gstEst[[1]]),
                                           bstrps)),
                    gstEstHed = array(unlist(gstEstHed),
                                      dim = c(nrow(gstEstHed[[1]]),
                                              nrow(gstEstHed[[1]]),
                                              bstrps)),
                    dEst = array(unlist(dEst),
                                 dim = c(nrow(dEst[[1]]),
                                         nrow(dEst[[1]]),
                                         bstrps)),
                    theta = array(unlist(theta),
                                  dim = c(nrow(theta[[1]]),
                                          nrow(theta[[1]]),
                                          bstrps)))
      
    } else {
      stats <- list(gstEst = array(unlist(gstEst),
                                   dim = c(nrow(gstEst[[1]]),
                                           nrow(gstEst[[1]]),
                                           bstrps)),
                    gstEstHed = array(unlist(gstEstHed),
                                      dim = c(nrow(gstEstHed[[1]]),
                                              nrow(gstEstHed[[1]]),
                                              bstrps)),
                    dEst = array(unlist(dEst),
                                 dim = c(nrow(dEst[[1]]),
                                         nrow(dEst[[1]]),
                                         bstrps)))
    }
    # tidy up
    if(fst){
      rm(dEst, gstEst, gstEstHed, theta)
      z <- gc()
      rm(z) 
    } else {
      # tidy up
      z <- gc()
      rm(z) 
    }
    #       # convert locus stats into arrays for CI calculations
    #       npops <- accDat$npops
    #       nloci <- accDat$nloci
    #       if(fst){
    #         locStats <- list(
    #           # nei's Gst
    #           gstLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,1])
    #           })), dim = c(npops, npops, nloci, bstrps)),
    #           # Hedrick's Gst
    #           gstHedLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,2])
    #           })), dim = c(npops, npops, nloci, bstrps)),
    #           # Jost's D
    #           dJostLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,3])
    #           })), dim = c(npops, npops, nloci, bstrps)),
    #           # Weir & Cockerham's Fst
    #           thetaLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,4])
    #           })), dim = c(npops, npops, nloci, bstrps))
    #         )
    #       } else {
    #         locStats <- list(
    #           # nei's Gst
    #           gstLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,1])
    #           })), dim = c(npops, npops, nloci, bstrps)),
    #           # Hedrick's Gst
    #           gstHedLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,2])
    #           })), dim = c(npops, npops, nloci, bstrps)),
    #           # Jost's D
    #           dJostLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,3])
    #           })), dim = c(npops, npops, nloci, bstrps))
    #         )
    #       }
    #       # convert locus bs stats into seperate loci
    #       locStats <- lapply(locStats, function(x){
    #         lapply(apply(x, 3, list), function(y){
    #           return(y[[1]])
    #         })
    #       })
    #       
    #       # calculate bias corrected CI
    #       
    #       biasCor <- function(param, bs_param){
    #         mnBS <- apply(bs_param, c(1,2), mean, na.rm = TRUE)
    #         mnBS[is.nan(mnBS)] <- NA
    #         mnBS <- mnBS - param
    #         bs_param <- sweep(bs_param, c(1:2), mnBS, "-")
    #         return(bs_param)
    #       }
    #       biasFix <- lapply(1:length(locstats), function(i){
    #         mapply(FUN = biasCor, param = locstats[[i]], bs_param = locStats[[i]],
    #                SIMPLIFY = FALSE)
    #       })
    #       # works well
    #       if (para && para_pack) {
    #         library(parallel)
    #         cl <- makeCluster(detectCores())
    #         bcLocLCI <- parLapply(cl, biasFix, function(x){
    #           loc <- lapply(x, function(y){
    #             apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
    #           })
    #           return(loc)
    #         })
    #         bcLocUCI <- parLapply(cl, biasFix, function(x){
    #           loc <- lapply(x, function(y){
    #             apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
    #           })
    #           return(loc)
    #         })
    #         stopCluster(cl)
    #       } else {
    #         bcLocLCI <- lapply(biasFix, function(x){
    #           loc <- lapply(x, function(y){
    #             apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
    #           })
    #           return(loc)
    #         })
    #         bcLocUCI <- lapply(biasFix, function(x){
    #           loc <- lapply(x, function(y){
    #             apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
    #           })
    #           return(loc)
    #         })
    #       }
    #       # organise data into output format
    #       # define function
    #       dfSort <- function(act, Low, High, pw_nms){
    #         df <- data.frame(actual = as.vector(act[lower.tri(act)]),
    #                          lower = as.vector(Low[lower.tri(Low)]),
    #                          upper = as.vector(High[lower.tri(High)]))
    #         rownames(df) <- pw_nms
    #         df[is.nan(as.matrix(df))] <- NA
    #         return(df)
    #       }
    #       pwLocOutput <- lapply(1:length(locstats), function(i){
    #         mapply(FUN = dfSort, act = locstats[[i]], Low = bcLocLCI[[i]],
    #                High = bcLocUCI[[i]], MoreArgs = list(pw_nms = pw_nms),
    #                SIMPLIFY = FALSE)
    #       })
    #       if(fst){
    #         names(pwLocOutput) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
    #       } else {
    #         names(pwLocOutput) <- c("gstEst", "gstEstHed", "djostEst")
    #       }
    #       for(i in 1:length(pwLocOutput)){
    #         names(pwLocOutput[[i]]) <- accDat$locus_names
    #       }
    
    
    pwMatListOut <- list(gstPW, gstHPW, dPW, thetaPW)
    # add names to pwMatListOut
    names(pwMatListOut) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
    
    
    
    
    # organise data
    # calculate bias for cis
    biasCalc <- function(param, bs_param, pw){
      #bias <- param
      for(i in 1:ncol(pw)){
        dat <- bs_param[pw[2,i], pw[1,i], ]
        t0 <- param[pw[2,i], pw[1,i]]
        mnBS <- mean(dat , na.rm = TRUE) - t0
        bs_param[pw[2,i], pw[1,i], ] <- bs_param[pw[2,i], pw[1,i], ] - mnBS
      }
      return(bs_param)
    }
    
    # try adjusting bootstrapped estimate using bias
    
    bcStats <- mapply(biasCalc, param = pwMatListOut, bs_param = stats, 
                      MoreArgs = list(pw = pw), SIMPLIFY = FALSE)
    
    ## Write results
    # we need three files for each of the four statistics calculated
    fileNms <- list()
    for(i in 1:4){
      fileNms[[i]] <- vector()
      fileNms[[i]][1] <- paste(of, names(stats)[i], 
                               "-actual.txt", sep = "")
      fileNms[[i]][2] <- paste(of, names(stats)[i], 
                               "-uncorrected.txt", sep = "")
      fileNms[[i]][3] <- paste(of, names(stats)[i], 
                               "-corrected.txt", sep = "")
    }
    # open file connection and write results
    for(i in 1:4){
      actual <- file(fileNms[[i]][1], "w")
      uncor <- file(fileNms[[i]][2], "w")
      corr <- file(fileNms[[i]][3], "w")
      # write standard statistics
      pwMatListOut[[i]][is.na(pwMatListOut[[i]])] <- ""
      for(j in 1:nrow(pwMatListOut[[i]])){
        cat(pwMatListOut[[i]][j,], "\n", file = actual, sep = "\t")
      }
      close(actual)
      # write standard bootstraps
      # generate a bootstrap list
      out1 <- lapply(apply(stats[[i]], 3, list), "[[", 1)
      out1 <- do.call("rbind", out1)
      out1[is.na(out1)] <- ""
      out1 <- out1[,-ncol(out1)]
      for(j in 1:nrow(out1)){
        cat(out1[j,], "\n", file = uncor, sep = "\t")
      }
      close(uncor)
      rm(out1)
      out2 <- lapply(apply(bcStats[[i]], 3, list), "[[", 1)
      out2 <- do.call("rbind", out2)
      out2[is.na(out2)] <- ""
      out2 <- out2[,-ncol(out2)]
      for(j in 1:nrow(out2)){
        cat(out2[j,], "\n", file = corr, sep = "\t")
      }
      close(corr)
    }
  }
}
################################################################################
# writeBoot end                                                                #
################################################################################