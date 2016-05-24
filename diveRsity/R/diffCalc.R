################################################################################
# NEW STYLE FUNCTIONS                                                          #
################################################################################
## Uncomment for non c++ version
#' ### diffCalcRcpp: an Rcpp version of diveRsity::diffCalc
#'
#' __Kevin Keenan__ (2014)
#'
#'
#' ```{r echo = FALSE, message = TRUE, cache = FALSE}
#'  library(rCharts, warn = FALSE)
#'  knitr::opts_knit$set(self.contained=TRUE)
#'  knitr::opts_chunk$set(tidy = FALSE, message = TRUE, echo=TRUE)
#' ```
#'
#' This code will be the basis for a faster more efficinet pairwise bootstrap
#' routine. The major difference will be that all of the data will be resampled
#' n time and the pairwise statistics will be calculated following this step.
#' The original method pairs the data then bootstraps each pair seperatly. This
#' is much less efficient.
#'
#' __Kevin Keenan__ (2014)
#'
diffCalc <- function(infile = NULL, outfile = NULL, fst = FALSE,
                     pairwise = FALSE, bs_locus = FALSE,
                     bs_pairwise = FALSE, boots = NULL,
                     ci_type = "individuals", alpha = 0.05,
                     para = FALSE){
  #' # Calculate diversity statistics functions
  #'
  #' __Kevin Keenan__ (2014)
  #'
  #' ### D_Jost
  dCalc <- function(ht, hs, n = NULL){
    if(!is.null(n)){
      return(((ht-hs)/(1-hs)) * (n/(n-1)))
    } else {
      return(2* ((ht-hs)/(1-hs)))
    }
  }
  #' ### Gst (Nei)
  gstCalc <- function(ht, hs){
    return((ht - hs)/ht)
  }
  #' ### G'st (Hedrick)
  GstCalc <- function(ht, hs, n = NULL){
    if(!is.null(n)){
      htmax <- ((n - 1) + hs)/n
    } else {
      htmax <- (1+hs)/2
    }
    return(((ht-hs)/ht)/((htmax-hs)/htmax))
  }
  
  #' ### G''st (Meirmans & Hedrick, 2011)
  GgstCalc <- function(ht, hs, n = NULL){
    if(!is.null(n)){
      return((n*(ht-hs))/(((n*ht) - hs)*(1-hs)))
    } else {
      return((2*(ht-hs))/(((2*ht) - hs)*(1-hs)))
    }
  }
  
  #' ### Locus and overall WC stats
  thetaCalc <- function(a, b, cdat){
    return(a/(a+b+cdat))
  }
  #' ### bias correction
  #'
  #' Function for correcting the bias associated with bootstrapped
  #' statistics
  #'
  #' __Kevin Keenan__ (2014)
  #'
  # Calculate CIs (bias corrected)
  diffCalcbCor <- function(act, bs){
    if(is.matrix(bs)){
      mn <- rowMeans(bs, na.rm = TRUE) - act
      mn[is.nan(mn)] <- NA
      bc <- mapply(`-`, split(bs, row(bs)), mn)
      return(bc)
    } else {
      mn <- mean(bs, na.rm = TRUE) - act
      mn[is.nan(mn)] <- NA
      return(bs - mn)
    }
  }
  #' # Calculate the harmonic sample size for pair of populations
  #'
  #' __Kevin Keenan__ (2014)
  #' As of 04/10/14 implemented as a C++ function
  #diffCalcHarm <- function(idt, pw){
  #  ps <- apply(pw, 2, function(x){
  #    return((0.5*((idt[x[1]]^-1) + idt[x[2]]^-1))^-1)
  #  })
  #  return(ps)
  #}
  #   infile <- "Test_data.txt"#Test_data
  #   outfile <- NULL
  #   fst = TRUE
  #   pairwise = TRUE
  #   bs_locus = TRUE
  #   bs_pairwise = TRUE
  #   boots = 100
  #   alpha = 0.05
  #   ci_type = "individuals"
  #   para = TRUE
  #   rgp <- diveRsity:::rgp
  #   statCalc <- diveRsity:::statCalc
  #   glbWCcpp <- diveRsity:::glbWCcpp
  #   varFunc <- diveRsity:::varFunc
  #   myTab <- diveRsity:::myTab
  #   pwHCalc <- diveRsity:::pwHCalc
  #   pwWCcpp <- diveRsity:::pwWCcpp
  #   diffCalcHarm <- diveRsity:::diffCalcHarm
  #   tabMerge <- diveRsity:::tabMerge
  #   pwTabMerge <- diveRsity:::pwTabMerge
  
  
  bs <- boots
  # read infile
  ip <- rgp(infile)
  # define some parameters
  
  # pop sizes
  ps = sapply(ip$genos, function(x){
    dim(x)[1]
  })
  
  # npops
  np = ncol(ip$af[[1]])
  
  # number of loci
  if(ci_type == "loci"){
    nl <- dim(ip$genos[[1]])[2]
  }
  # pairwise matrix index
  pw <- combn(np, 2)
  
  # set up parallel cluster
  if(para){
    
    ncor <- parallel::detectCores()
  }
  
  # create resample indexes
  if(ci_type == "individuals"){
    if(!is.null(bs)){
      idx <- lapply(1:bs, function(i){
        lapply(1:np, function(j){
          sample(ps[j], size = ps[j], replace = TRUE)
        })
      })
    }
  } else {
    idx <- lapply(1:bs, function(i){
      sample(nl, nl, replace = TRUE)
    })
  }
  
  #' ### Calculate point estimates (locus and global)
  #' Calculate gst, Gst, Fst and D for loci (across samples) and
  #' global (across loci and samples).
  #'
  
  #####################################
  # Point calculations
  #####################################
  preStats <- statCalc(rsDat = ip$genos, al = ip$af, fst = fst, bs = FALSE)
  # identify loci with unscored samples and fix for glbWCcpp
  #lapply
  # tabMerge [ silenced in place of C++ version, since 04/10/14 ]
  #tabMerge <- function(...){
  #  ip <- unlist(list(...))
  #  idx <- names(ip) != "NA"
  #  ip <- ip[idx]
  #  out <- sapply(split(ip, names(ip)), sum)
  #  if(length(out) == 0L){
  #    ot <- NA
  #    names(ot) <- "NA"
  #    return(ot)
  #  } else {
  #    return(out)
  #  }
  #}
  if(fst){
    # locus variance components (working)
    hsum <- lapply(preStats$hsum, tabMerge)
    # fix missing data loci
    hsum <- lapply(hsum, function(x){
      x[!(names(x) == "NA")]
    })
    varC <- mapply("glbWCcpp", hsum = hsum, af = preStats$alOut,
                   indtyp = preStats$indtyp, SIMPLIFY = FALSE)
    rm(hsum)
    # Locus F-statistics
    locFst <- sapply(varC, function(x){
      return(x$a/(x$a+x$b+x$c))
    })
    locFit <- sapply(varC, function(x){
      return(1 - (x$c/(x$a + x$b + x$c)))
    })
    locFis <- sapply(varC, function(x){
      return(1 - (x$c/(x$b + x$c)))
    })
    
    # global variance components
    glba <- mean(sapply(varC, "[[", "a"), na.rm = TRUE)
    glbb <- mean(sapply(varC, "[[", "b"), na.rm = TRUE)
    glbc <- mean(sapply(varC, "[[", "c"), na.rm = TRUE)
    
    # F-statistics
    glbTheta <- glba/(glba + glbb + glbc)
    glbFit <- 1 - (glbc/(glba + glbb + glbc))
    glbFis <- 1 - (glbc/(glbb + glbc))
  }
  
  # calculate harmonic mean sample size per locus
  hrm <- function(x){
    return(length(x)/(sum(1/x)))
  }
  harmLoc <- sapply(preStats$indtyp, hrm)
  # locus global stats (replace with Rcpp function)
  # calculate heterozygosities
  hthsLoc <- mapply(varFunc, af = preStats$alOut, sHarm = harmLoc,
                    SIMPLIFY = FALSE)
  # Jost's D (loci + global) #
  dLoc <- dCalc(ht = sapply(hthsLoc, "[[", "htEst"),
                hs = sapply(hthsLoc, "[[", "hsEst"),
                n = np)
  mnD <- mean(dLoc, na.rm = TRUE)
  varD <- var(dLoc, na.rm = TRUE)
  dGlb <- 1/((1/mnD)+varD*(1/mnD)^3)
  
  # gst (loci + global) #
  gLoc <- gstCalc(sapply(hthsLoc, "[[", "htEst"),
                  sapply(hthsLoc, "[[", "hsEst"))
  gGlb <- gstCalc(mean(sapply(hthsLoc, "[[", "htEst"), na.rm = TRUE),
                  mean(sapply(hthsLoc, "[[", "hsEst"), na.rm = TRUE))
  
  # Gst (loci + global) #
  GLoc <- GstCalc(sapply(hthsLoc, "[[", "htEst"),
                  sapply(hthsLoc, "[[", "hsEst"),
                  n = np)
  GGlb <- GstCalc(mean(sapply(hthsLoc, "[[", "htEst"), na.rm = TRUE),
                  mean(sapply(hthsLoc, "[[", "hsEst"), na.rm = TRUE),
                  n = np)
  
  # G''st (loci + global) #
  GGLoc <- GgstCalc(sapply(hthsLoc, "[[", "htEst"),
                    sapply(hthsLoc, "[[", "hsEst"),
                    n = np)
  GGGlb <- GgstCalc(mean(sapply(hthsLoc, "[[", "htEst"), na.rm = TRUE),
                    mean(sapply(hthsLoc, "[[", "hsEst"), na.rm = TRUE),
                    n = np)
  
  #####################################
  # END
  #####################################
  
  
  #' ### Create Est structure
  #'
  
  #####################################
  # Arange points estimates for output
  #####################################
  # non mem prob
  if(!fst){
    est <- data.frame(loci = c(ip$locs, "Global"),
                      gst = round(c(gLoc, gGlb), 4),
                      Gst = round(c(GLoc, GGlb), 4),
                      GGst = round(c(GGLoc, GGGlb), 4),
                      D = round(c(dLoc, dGlb), 4))
    rownames(est) <- NULL
  } else {
    est <- data.frame(loci = c(ip$locs, "Global"),
                      gst = round(c(gLoc, gGlb), 4),
                      Gst = round(c(GLoc, GGlb), 4),
                      GGst = round(c(GGLoc, GGGlb), 4),
                      D = round(c(dLoc, dGlb), 4),
                      Fis = round(c(locFis, glbFis), 4),
                      Fst = round(c(locFst, glbTheta), 4),
                      Fit = round(c(locFit, glbFit), 4))
    rownames(est) <- NULL
  }
  
  #####################################
  # END
  #####################################
  
  
  #' #### Output format (preview)
  #' ```{r, cache=FALSE, echo=FALSE, results='asis'}
  #'  dt <- Datatables$new()
  #'  dt$addTable(est)
  #'  dt$print("table1", include_assets = TRUE, cdn = TRUE)
  #'
  #'```
  #'
  #' <br>
  #'
  #' ### Calculate locus and global bootstraps
  #' Calculate the bootstrapped \(95\%\) Confidence intervals for locus and
  #' global gst, Gst, Fst and D.
  #'
  #'
  #####################################
  # Bootstrap data
  ####################################
  if(bs_locus || bs_pairwise){
    # pre-statistics
    # to run the bootstraps
    # replace al with 0
    al <- lapply(ip$af, function(x){
      nms <- rownames(x)
      x[x != 0] <- 0
      rownames(x) <- nms
      return(x)
    })
    if(para){
      cl <- parallel::makeCluster(ncor)
      #clusterExport(cl, c("myTab", "al"), envir = environment())
      bsDat <- parallel::parLapply(cl, idx, statCalc, rsDat = ip$genos, al = al,
                         fst = fst, ci_type = ci_type)
      parallel::stopCluster(cl)
    } else {
      system.time({
        bsDat <- lapply(idx, statCalc, rsDat = ip$genos, al = al, fst = fst,
                        ci_type = ci_type)
      })
    }
    indtyp <- lapply(bsDat, "[[", "indtyp")
    # clean up
    rm(idx)
  }
  #####################################
  # END
  ####################################
  
  #####################################
  # Locus and global bootstrapping
  #####################################
  
  # added Rcpp
  if(bs_locus){
    
    #####################################
    # W&C locus stats
    #####################################
    # Low memory overhead (Takes 12.1 sec for 1000 bs or "pw_test.gen")
    if(fst){
      
      bsVarC <- lapply(bsDat, function(x){
        hsum <- lapply(x$hsum, tabMerge)
        # fix missing data loci
        hsum <- lapply(hsum, function(x){
          x[!(names(x) == "NA")]
        })
        stat <- mapply(glbWCcpp, hsum = hsum, af = x$alOut,
                       indtyp = x$indtyp, SIMPLIFY = FALSE)
        bsFst <- sapply(stat, function(y){
          return(y$a / (y$a + y$b + y$c))
        })
        bsFit <- sapply(stat, function(y){
          return(1 - (y$c / (y$a + y$b + y$c)))
        })
        bsFis <- sapply(stat, function(y){
          return(1 - (y$c / (y$b + y$c)))
        })
        glba <- mean(sapply(stat, "[[", "a"), na.rm = TRUE)
        glbb <- mean(sapply(stat, "[[", "b"), na.rm = TRUE)
        glbc <- mean(sapply(stat, "[[", "c"), na.rm = TRUE)
        glbst <- glba/(glba + glbb + glbc)
        glbit <- 1 - (glbc/(glba + glbb + glbc))
        glbis <- 1 - (glbc/(glbb + glbc))
        list(bsFstLoc = bsFst, bsFitLoc = bsFit, bsFisLoc = bsFis,
             bsFstAll = glbst, bsFitAll = glbit, bsFisAll = glbis)
      })
      
      
      # Extract bootstrap statistics
      bsFstL <- sapply(bsVarC, "[[", "bsFstLoc")
      bsFisL <- sapply(bsVarC, "[[", "bsFisLoc")
      bsFitL <- sapply(bsVarC, "[[", "bsFitLoc")
      bsFstA <- sapply(bsVarC, "[[", "bsFstAll")
      bsFisA <- sapply(bsVarC, "[[", "bsFisAll")
      bsFitA <- sapply(bsVarC, "[[", "bsFitAll")
      # Calculate bias corrected bs values
      bsFstL <- diffCalcbCor(locFst, bsFstL)
      bsFisL <- diffCalcbCor(locFis, bsFisL)
      bsFitL <- diffCalcbCor(locFit, bsFitL)
      bsFstA <- diffCalcbCor(glbTheta, bsFstA)
      bsFisA <- diffCalcbCor(glbFis, bsFisA)
      bsFitA <- diffCalcbCor(glbFit, bsFitA)
      # Calculate 95% CIs
      LCI <- data.frame(Fst = apply(bsFstL, 2, quantile, probs = alpha/2,
                                    na.rm = TRUE),
                        Fis = apply(bsFisL, 2, quantile, probs = alpha/2,
                                    na.rm = TRUE),
                        Fit = apply(bsFitL, 2, quantile, probs = alpha/2,
                                    na.rm = TRUE))
      UCI <- data.frame(Fst = apply(bsFstL, 2, quantile, probs = 1-(alpha/2),
                                    na.rm = TRUE),
                        Fis = apply(bsFisL, 2, quantile, probs = 1-(alpha/2),
                                    na.rm = TRUE),
                        Fit = apply(bsFitL, 2, quantile, probs = 1-(alpha/2),
                                    na.rm = TRUE))
      # clean up
      #rm(bsFstL, bsFisL, bsFitL, bsVarC)
      #z <- gc(reset = TRUE)
      #rm(z)
    }
    #####################################
    # END W&C locus stats
    #####################################
    
    
    #####################################
    # Het based locus stats
    #####################################
    # Calculate heterozygosity based stats
    # harmonic sample sizes
    
    allHarm <- lapply(bsDat, function(x){
      sapply(x$indtyp, function(y){
        return(((1/length(y))*(sum(y^-1)))^-1)
      })
    })
    
    # add allHarm to bsDat
    listAdd <- function(bsDat, sHarm){
      bsDat$nHarm <- sHarm
      return(bsDat)
    }
    
    bsDat <- mapply(listAdd, bsDat = bsDat, sHarm = allHarm, SIMPLIFY = FALSE)
    
    # Very low memory over head (takes 1.6 sec for 1000 boostraps using
    # "pw_test.gen")
    # Calculate stats
    # Parallel is not worth the memory overhead.
    hetVar <- lapply(bsDat, function(x){
      stat <- mapply("varFunc", af = x$alOut, sHarm = x$nHarm,
                     SIMPLIFY = FALSE)
      ht <- sapply(stat, "[[", "htEst")
      hs <- sapply(stat, "[[", "hsEst")
      n <- ncol(x$alOut[[1]])
      bsDLoc <- dCalc(ht = ht, hs = hs, n = n)
      mnD <- mean(bsDLoc, na.rm = TRUE)
      varD <- var(bsDLoc, na.rm = TRUE)
      bsDAll <- 1/((1/mnD) + varD * (1/mnD)^3)
      bsgLoc <- gstCalc(ht = ht, hs = hs)
      bsgAll <- gstCalc(ht = mean(ht, na.rm = TRUE),
                        hs = mean(hs, na.rm = TRUE))
      bsGLoc <- GstCalc(ht = ht, hs = hs, n = n)
      bsGAll <- GstCalc(ht = mean(ht, na.rm = TRUE),
                        hs = mean(hs, na.rm = TRUE), n = n)
      bsGGLoc <- GgstCalc(ht = ht, hs = hs, n = n)
      bsGGAll <- GgstCalc(ht = mean(ht, na.rm = TRUE),
                          hs = mean(hs, na.rm = TRUE), n = n)
      list(bsDLoc = bsDLoc, bsgLoc = bsgLoc, bsGLoc = bsGLoc,
           bsDAll = bsDAll, bsgAll = bsgAll, bsGAll = bsGAll,
           bsGGLoc = bsGGLoc, bsGGAll = bsGGAll)
    })
    
    
    # Extract bootstrap statistics
    bsDL <- sapply(hetVar, "[[", "bsDLoc")
    bsgL <- sapply(hetVar, "[[", "bsgLoc")
    bsGL <- sapply(hetVar, "[[", "bsGLoc")
    bsGGL <- sapply(hetVar, "[[", "bsGGLoc")
    bsDA <- sapply(hetVar, "[[", "bsDAll")
    bsgA <- sapply(hetVar, "[[", "bsgAll")
    bsGA <- sapply(hetVar, "[[", "bsGAll")
    bsGGA <- sapply(hetVar, "[[", "bsGGAll")
    # Calculate bias corrected bs values
    bsDL <- diffCalcbCor(dLoc, bsDL)
    bsgL <- diffCalcbCor(gLoc, bsgL)
    bsGL <- diffCalcbCor(GLoc, bsGL)
    bsGGL <- diffCalcbCor(GGLoc, bsGGL)
    bsDA <- diffCalcbCor(dGlb, bsDA)
    bsgA <- diffCalcbCor(gGlb, bsgA)
    bsGA <- diffCalcbCor(GGlb, bsGA)
    bsGGA <- diffCalcbCor(GGGlb, bsGGA)
    
    # Calculate 95% CIs
    if(fst){
      # lower
      LCI$D <- apply(bsDL, 2, quantile, probs = alpha/2, na.rm = TRUE)
      LCI$gst <- apply(bsgL, 2, quantile, probs = alpha/2, na.rm = TRUE)
      LCI$Gst <- apply(bsGL, 2, quantile, probs = alpha/2, na.rm = TRUE)
      LCI$GGst <- apply(bsGGL, 2, quantile, probs = alpha/2, na.rm = TRUE)
      # upper
      UCI$D <- apply(bsDL, 2, quantile, probs = 1-(alpha/2), na.rm = TRUE)
      UCI$gst <- apply(bsgL, 2, quantile, probs = 1-(alpha/2), na.rm = TRUE)
      UCI$Gst <- apply(bsGL, 2, quantile, probs = 1-(alpha/2), na.rm = TRUE)
      UCI$GGst <- apply(bsGGL, 2, quantile, probs = 1-(alpha/2), na.rm = TRUE)
    } else {
      # lower
      LCI <- data.frame(D = apply(bsDL, 2, quantile,
                                  probs = alpha/2, na.rm = TRUE))
      LCI$gst <- apply(bsgL, 2, quantile, probs = alpha/2, na.rm = TRUE)
      LCI$Gst <- apply(bsGL, 2, quantile, probs = alpha/2, na.rm = TRUE)
      LCI$GGst <- apply(bsGGL, 2, quantile, probs = alpha/2, na.rm = TRUE)
      # upper
      UCI <- data.frame(D = apply(bsDL, 2, quantile,
                                  probs = 1-(alpha/2), na.rm = TRUE))
      UCI$gst <- apply(bsgL, 2, quantile, probs = 1-(alpha/2), na.rm = TRUE)
      UCI$Gst <- apply(bsGL, 2, quantile, probs = 1-(alpha/2), na.rm = TRUE)
      UCI$GGst <- apply(bsGGL, 2, quantile, probs = 1-(alpha/2), na.rm = TRUE)
    }
    
    #####################################
    # END: Het based locus stats
    #####################################
    
    # global CIs
    if(fst){
      glbBS <- data.frame(gst = bsgA, Gst = bsGA, GGst = bsGGA,
                          d = bsDA, fst = bsFstA, fit = bsFitA,
                          fis = bsFisA)
      
      rm(bsFstA, bsFisA, bsFitA, bsDA, bsgA, bsGA, bsGGA)
    } else {
      glbBS <- data.frame(gst = bsgA, Gst = bsGA, GGst = bsGGA, d = bsDA)
      rm(bsDA, bsgA, bsGA, bsGGA)
    }
    glbLCI <- apply(glbBS, 2, quantile, probs = alpha/2, na.rm = TRUE)
    glbUCI <- apply(glbBS, 2, quantile, probs = 1-(alpha/2), na.rm = TRUE)
    if(fst){
      glbOut <- data.frame(stat = c("gst", "Gst", "GGst", "D", "Fst",
                                    "Fit", "Fis"),
                           actual = c(gGlb, GGlb, GGGlb, dGlb, glbTheta,
                                      glbFit, glbFis),
                           lower = glbLCI, upper = glbUCI, row.names = NULL)
    } else {
      glbOut <- data.frame(stat = c("gst", "Gst", "GGst", "D"),
                           actual = c(gGlb, GGlb, GGGlb, dGlb),
                           lower = glbLCI, upper = glbUCI, row.names = NULL)
    }
  }
  
  
  #####################################
  # END
  #####################################
  
  popnms <- sapply(ip$indnms, "[[", 1)
  pwpops <- paste(popnms[pw[1,]], " vs ", popnms[pw[2,]])
  
  ######################################
  # Pairwise bootstrapping
  #####################################
  #' ### Calculate PW bootstrap statistics
  
  if(pairwise || bs_pairwise){
    #' ### Define heterozygosity based stats (Gst, G'st, D_Jost)
    
    # Calculate non-bootstrap parameters
    
    # calculate harmonic sample sizes from preStats
    preNharm <- lapply(preStats$indtyp, diffCalcHarm, pw = pw-1)
    # add to preStats
    preStats$sHarm <- preNharm
    # calculate heterozygosities
    nbshths <- mapply(pwHCalc, af = preStats$alOut, sHarm = preStats$sHarm,
                      MoreArgs = list(pw = pw-1), SIMPLIFY = FALSE)
    # convert to locus focus
    nbshths <- lapply(c("hsEst", "htEst"), function(x){
      lapply(nbshths, "[[", x)
    })
    names(nbshths) <- c("hsEst", "htEst")
    # calculate mean ht and hs
    hsMn <- vapply(1:ncol(pw), function(i){
      mean(vapply(nbshths$hsEst, "[[", i, FUN.VALUE = numeric(1)), na.rm = TRUE)
    }, FUN.VALUE = numeric(1))
    
    htMn <- vapply(1:ncol(pw), function(i){
      mean(vapply(nbshths$htEst, "[[", i, FUN.VALUE = numeric(1)), na.rm = TRUE)
    }, FUN.VALUE = numeric(1))
    
    # calculate non-bootstrap pairwise stats
    
    # Jost's D
    pwDLoc <- mapply(dCalc, ht = nbshths$htEst, hs = nbshths$hsEst)
    pwDall <- apply(pwDLoc, 1, function(x){
      mnD <- mean(x, na.rm = TRUE)
      vrD <- var(x, na.rm = TRUE)
      1/((1/mnD)+((vrD*((1/mnD)^3))))
    })
    
    # gst
    pwgLoc <- mapply(gstCalc, ht = nbshths$htEst, hs = nbshths$hsEst)
    pwgAll <- gstCalc(ht = htMn, hs = hsMn)
    
    # Gst
    pwGLoc <- mapply(GstCalc, ht = nbshths$htEst, hs = nbshths$hsEst)
    pwGAll <- GstCalc(ht = htMn, hs = hsMn)
    
    # G''st
    pwGGLoc <- mapply(GgstCalc, ht = nbshths$htEst, hs = nbshths$hsEst)
    pwGGAll <- GgstCalc(ht = htMn, hs = hsMn)
    
    if(fst){
      # theta
      # calculate standard statistics (non-bootstrap)
      hsum <- lapply(preStats$hsum, function(x){
        pwTabMerge(x, pw-1)
      })
      # fix missing data loci
      # fix missing data loci
      hsum <- lapply(hsum, function(x){
        x <- lapply(x, function(y){
          y <- y[!(names(y) == "NA")]
          return(y)
        })
      })
      pwVar <- mapply("pwWCcpp", hsum1 = hsum, af1 = preStats$alOut,
                      indtyp1 = preStats$indtyp, MoreArgs = list(pw = pw-1),
                      SIMPLIFY = FALSE)
      # convert to locus focus
      pwVar <- lapply(c("a", "b", "c"), function(x){
        return(lapply(pwVar, "[[", x))
      })
      names(pwVar) <- c("a", "b", "cdat")
      
      # loci
      pwFstLoc <- mapply(thetaCalc, a = pwVar$a, b = pwVar$b, cdat = pwVar$cdat)
      
      # Global
      mnA <- vapply(1:ncol(pw), function(i){
        mean(vapply(pwVar$a, "[", i, FUN.VALUE = numeric(1)), na.rm = TRUE)
      }, FUN.VALUE = numeric(1))
      mnB <- vapply(1:ncol(pw), function(i){
        mean(vapply(pwVar$b, "[", i, FUN.VALUE = numeric(1)), na.rm = TRUE)
      }, FUN.VALUE = numeric(1))
      mnC <- vapply(1:ncol(pw), function(i){
        mean(vapply(pwVar$cdat, "[", i, FUN.VALUE = numeric(1)), na.rm = TRUE)
      }, FUN.VALUE = numeric(1))
      pwFstAll <- mnA/(mnA + mnB + mnC)
    }
  }
  
  if(bs_pairwise){
    # calculate harmonic sample sizes (list[bs]:list[loc]:vector[pw])
    sHarm <- lapply(indtyp, function(x){
      lapply(x, diffCalcHarm, pw = pw-1)
    })
    rm(indtyp)
    # combine sHarm with bsDat
    listAdd <- function(bsDat, sHarm){
      bsDat$sHarm <- sHarm
      bsDat$nHarm <- NULL
      return(bsDat)
    }
    
    bsDat <- mapply(listAdd, bsDat = bsDat, sHarm = sHarm, SIMPLIFY = FALSE)
    rm(sHarm)
    
    # calculate data
    hths <- lapply(bsDat, function(x){
      lapply(1:length(x$alOut), function(i){
        pwHCalc(x$alOut[[i]], sHarm = x$sHarm[[i]], pw = pw-1)
      })
    })
    
    # convert hths to stat focus
    hths <- lapply(hths, function(x){
      hsEst <- lapply(x, "[[", 1)
      htEst <- lapply(x, "[[", 2)
      #hs <- lapply(x, "[[", 3)
      #ht <- lapply(x, "[[", 4)
      return(list(hsEst = hsEst, htEst = htEst))#, hs = hs, ht = ht))
    })
    
    
    # calculate locus estimator (parallel is slower!)
    pwDLocbs <- sapply(hths, function(x){
      return(mapply(dCalc, ht = x$htEst, hs = x$hsEst, SIMPLIFY = TRUE))
    }, simplify = "array")
    # overall d bootstraps
    # Silence para since the benefits are minimal
    #if(para){
    #  cl <- makeCluster(ncor)
    #  pwDAllbs <- parApply(cl, pwDLocbs, c(1,3), function(x){
    #    mn <- mean(x, na.rm = TRUE)
    #    vr <- var(x, na.rm = TRUE)
    #    return(1/((1/mn)+vr*(1/mn)^3))
    #  })
    #  stopCluster(cl)
    #} else {
    pwDAllbs <- apply(pwDLocbs, c(1,3), function(x){
      mn <- mean(x, na.rm = TRUE)
      vr <- var(x, na.rm = TRUE)
      1/((1/mn)+((vr*((1/mn)^3))))
    })
    #}
    rm(pwDLocbs)
    
    # calculate locus estimtor
    #     pwgLocbs <- sapply(hths, function(x){
    #       return(mapply(gstCalc, ht = x$htEst, hs = x$hsEst, SIMPLIFY = TRUE))
    #     }, simplify = "array")
    # overall gst (returns pw values in rows and bootstraps reps in cols)
    pwgAllbs <- sapply(hths, function(x){
      ht <- rowMeans(matrix(unlist(x$htEst), ncol = length(x$htEst)))
      hs <- rowMeans(matrix(unlist(x$hsEst), ncol = length(x$hsEst)))
      return(gstCalc(ht, hs))
    }, simplify = "array")
    
    # calculate locus estimator
    #     pwGLocbs <- sapply(hths, function(x){
    #       return(mapply(GstCalc, ht = x$htEst, hs = x$hsEst, SIMPLIFY = TRUE))
    #     }, simplify = "array")
    # overall Gst (returns pw values in rows and bootstraps reps in cols)
    pwGAllbs <- sapply(hths, function(x){
      ht <- rowMeans(matrix(unlist(x$htEst), ncol = length(x$htEst)))
      hs <- rowMeans(matrix(unlist(x$hsEst), ncol = length(x$hsEst)))
      return(GstCalc(ht, hs))
    }, simplify = "array")
    
    # overall G''st (returns pw values in rows and bootstraps reps in cols)
    pwGGAllbs <- sapply(hths, function(x){
      ht <- rowMeans(matrix(unlist(x$htEst), ncol = length(x$htEst)))
      hs <- rowMeans(matrix(unlist(x$hsEst), ncol = length(x$hsEst)))
      return(GgstCalc(ht, hs))
    }, simplify = "array")
    
    #' ### Calculate Weir & Cockerham variance components (v. slow, try to
    #' write a C++ tabMerge function)
    if(fst){
      # Calculate bootstrap F-statistics
      wcVar <- lapply(bsDat, function(x){
        x$hsum <- lapply(x$hsum, pwTabMerge, pw = pw - 1)
        # Calculate pairwise Fst
        return(mapply("pwWCcpp", hsum1 = x$hsum, indtyp1 = x$indtyp,
                      af1 = x$alOut, MoreArgs = list(pw = pw-1),
                      SIMPLIFY = FALSE))
      })
      # Calculate bootstrap theta
      #       pwFstLocbs <- sapply(wcVar, function(x){
      #         a <- lapply(x, "[[", "a")
      #         b <- lapply(x, "[[", "b")
      #         cdat <- lapply(x, "[[", "cdat")
      #         return(mapply(thetaCalc, a = a, b = b, cdat = cdat,
      #                SIMPLIFY = TRUE))
      #       }, simplify = "array")
      pwFstAllbs <- sapply(wcVar, function(x){
        a <- rowMeans(sapply(x, "[[", "a"))
        b <- rowMeans(sapply(x, "[[", "b"))
        cdat <- rowMeans(sapply(x, "[[", "c"))
        return(thetaCalc(a, b, cdat))
      }, simplify = "array")
      # clean up
      rm(wcVar, bsDat)
      z <- gc()
      
      
      
      #################################
      # calculate WC bias corrected CIs
      #################################
      
      # loci
      #             pwFstLocbs <- mapply(diffCalcbCor,
      #                                  act = split(pwFstLoc, row(pwFstLoc)),
      #                                  bs = lapply(1:dim(pwFstLocbs)[1],
      #                                              function(i){
      #                                                return(pwFstLocbs[i,,])}),
      #                                  SIMPLIFY = "array")
      #             # transpose array
      #             pwFstLocbs <- aperm(pwFstLocbs)
      #
      #             pwFstLCI <- apply(pwFstLocbs, c(1,2), quantile, prob = 0.025,
      #                               na.rm = TRUE)
      #             pwFstUCI <- apply(pwFstLocbs, c(1,2), quantile, prob = 0.975,
      #                               na.rm = TRUE)
      
      # global
      pwFstAllbs <- t(mapply(diffCalcbCor, act = pwFstAll,
                             bs = split(pwFstAllbs, row(pwFstAllbs))))
      pwFstAllLCI <- apply(pwFstAllbs, 1, quantile, prob = 0.025,
                           na.rm = TRUE)
      pwFstAllUCI <- apply(pwFstAllbs, 1, quantile, prob = 0.975,
                           na.rm = TRUE)
      
      #################################
      # END
      #################################
    }
    
    #' ### Heterozygosity based CIs
    
    #################################
    # Calculate het based bs stats
    #################################
    
    # Jost's D
    # loci
    #     pwDLocbs <- mapply(diffCalcbCor, act = split(pwDLoc, row(pwDLoc)),
    #                        bs = lapply(1:dim(pwDLocbs)[1], function(i){
    #                          return(pwDLocbs[i,,])}), SIMPLIFY = "array")
    #     pwDLocbs <- aperm(pwDLocbs)
    #     # CIs
    #     pwDLocLCI <- apply(pwDLocbs, c(1,2), quantile, prob = 0.025,
    #                        na.rm = TRUE)
    #     pwDLocUCI <- apply(pwDLocbs, c(1,2), quantile, prob = 0.975,
    #                        na.rm = TRUE)
    
    # global
    pwDAllbs <- t(mapply(diffCalcbCor, act = pwDall,
                         bs = split(pwDAllbs, row(pwDAllbs))))
    # CIs
    pwDAllLCI <- apply(pwDAllbs, 1, quantile, prob = 0.025, na.rm = TRUE)
    pwDAllUCI <- apply(pwDAllbs, 1, quantile, prob = 0.975, na.rm = TRUE)
    
    # gst
    # loci
    #     pwgLocbs <- mapply(diffCalcbCor, act = split(pwgLoc, row(pwgLoc)),
    #                        bs = lapply(1:dim(pwgLocbs)[1], function(i){
    #                          return(pwgLocbs[i,,])}), SIMPLIFY = "array")
    #     pwgLocbs <- aperm(pwgLocbs)
    #     # CIs
    #     pwgLocLCI <- apply(pwgLocbs, c(1,2), quantile, prob = 0.025,
    #                        na.rm = TRUE)
    #     pwgLocUCI <- apply(pwgLocbs, c(1,2), quantile, prob = 0.975,
    #                        na.rm = TRUE)
    
    # global
    pwgAllbs <- t(mapply(diffCalcbCor, act = pwgAll,
                         bs = split(pwgAllbs, row(pwgAllbs))))
    # CIs
    pwgAllLCI <- apply(pwgAllbs, 1, quantile, prob = 0.025, na.rm = TRUE)
    pwgAllUCI <- apply(pwgAllbs, 1, quantile, prob = 0.975, na.rm = TRUE)
    
    # Gst
    # loci
    #     pwGLocbs <- mapply(diffCalcbCor, act = split(pwGLoc, row(pwGLoc)),
    #                        bs = lapply(1:dim(pwGLocbs)[1], function(i){
    #                          return(pwGLocbs[i,,])}), SIMPLIFY = "array")
    #     pwGLocbs <- aperm(pwGLocbs)
    #     # CIs
    #     pwGLocLCI <- apply(pwGLocbs, c(1,2), quantile, prob = 0.025,
    #                        na.rm = TRUE)
    #     pwGLocUCI <- apply(pwGLocbs, c(1,2), quantile, prob = 0.975,
    #                        na.rm = TRUE)
    
    # global
    pwGAllbs <- t(mapply(diffCalcbCor, act = pwGAll,
                         bs = split(pwGAllbs, row(pwGAllbs))))
    # CIs
    pwGAllLCI <- apply(pwGAllbs, 1, quantile, prob = 0.025, na.rm = TRUE)
    pwGAllUCI <- apply(pwGAllbs, 1, quantile, prob = 0.975, na.rm = TRUE)
    
    # G''st
    
    # global
    pwGGAllbs <- t(mapply(diffCalcbCor, act = pwGAll,
                          bs = split(pwGGAllbs, row(pwGGAllbs))))
    # CIs
    pwGGAllLCI <- apply(pwGGAllbs, 1, quantile, prob = 0.025, na.rm = TRUE)
    pwGGAllUCI <- apply(pwGGAllbs, 1, quantile, prob = 0.975, na.rm = TRUE)
    
    
    
    #################################
    # END
    #################################
    # stop cluster
    #if(para){
    #  stopCluster(cl)
    #}
  }
  
  
  
  
  #' ### organise outputs
  # est is already available from above
  
  #################################
  # Global bootstrap results
  #################################
  op <- list(std_stats = est)
  if(bs_locus){
    op$global_bs <- data.frame(stat = glbOut[,1], round(glbOut[,-1], 4))
  }
  #################################
  # END
  #################################
  
  
  #################################
  # Locus Confidence intervals
  #################################
  if(bs_locus){
    if(fst){
      statnms <- c("gst", "Gst", "GGst", "D", "Fst", "Fis", "Fit")
    } else {
      statnms <- c("gst", "Gst", "GGst", "D")
    }
    locCI <- lapply(statnms, function(x){
      return(data.frame(locus = ip$locs,
                        actual = est[-nrow(est) ,x],
                        lower = round(LCI[, x], 4),
                        upper = round(UCI[, x], 4),
                        row.names = NULL))
    })
    names(locCI) <- statnms
    op$bs_locus <- locCI
    
    rm(locCI)
    z <- gc(reset = TRUE)
  }
  #################################
  # END
  #################################
  
  
  #################################
  # Pairwise stats
  #################################
  if(pairwise){
    # locus
    if(fst){
      # gst
      pwgLoc <- data.frame(round(t(pwgLoc), 4))
      dimnames(pwgLoc) <- list(ip$locs, pwpops)
      op$pw_locus$gst <- pwgLoc
      rm(pwgLoc)
      # Gst
      pwGLoc <- data.frame(round(t(pwGLoc), 4))
      dimnames(pwGLoc) <- list(ip$locs, pwpops)
      op$pw_locus$Gst <- pwGLoc
      rm(pwGLoc)
      #G''st
      pwGGLoc <- data.frame(round(t(pwGGLoc), 4))
      dimnames(pwGGLoc) <- list(ip$locs, pwpops)
      op$pw_locus$GGst <- pwGGLoc
      rm(pwGGLoc)
      # D
      pwDLoc <- data.frame(round(t(pwDLoc), 4))
      dimnames(pwDLoc) <- list(ip$locs, pwpops)
      op$pw_locus$D <- pwDLoc
      rm(pwDLoc)
      # Fst
      pwFstLoc <- data.frame(round(t(pwFstLoc), 4))
      dimnames(pwFstLoc) <- list(ip$locs, pwpops)
      op$pw_locus$Fst <- pwFstLoc
      rm(pwFstLoc)
    } else {
      # gst
      pwgLoc <- data.frame(round(t(pwgLoc), 4))
      dimnames(pwgLoc) <- list(ip$locs, pwpops)
      op$pw_locus$gst <- pwgLoc
      rm(pwgLoc)
      # Gst
      pwGLoc <- data.frame(round(t(pwGLoc), 4))
      dimnames(pwGLoc) <- list(ip$locs, pwpops)
      op$pw_locus$Gst <- pwGLoc
      rm(pwGLoc)
      #G''st
      pwGGLoc <- data.frame(round(t(pwGGLoc), 4))
      dimnames(pwGGLoc) <- list(ip$locs, pwpops)
      op$pw_locus$GGst <- pwGGLoc
      rm(pwGGLoc)
      # D
      pwDLoc <- data.frame(round(t(pwDLoc), 4))
      dimnames(pwDLoc) <- list(ip$locs, pwpops)
      op$pw_locus$D <- pwDLoc
      rm(pwDLoc)
    }
    # global
    if(fst){
      # gst
      op$pairwise <- list(gst = matrix(NA, nrow = np, ncol = np))
      op$pairwise$gst[lower.tri(op$pairwise$gst)] <- round(pwgAll, 4)
      dimnames(op$pairwise$gst) <- list(popnms, popnms)
      #rm(pwgAll)
      # Gst
      op$pairwise$Gst <- op$pairwise$gst
      op$pairwise$Gst[lower.tri(op$pairwise$Gst)] <- round(pwGAll, 4)
      dimnames(op$pairwise$Gst) <- list(popnms, popnms)
      # G''st
      op$pairwise$GGst <- op$pairwise$gst
      op$pairwise$GGst[lower.tri(op$pairwise$GGst)] <- round(pwGGAll, 4)
      dimnames(op$pairwise$GGst) <- list(popnms, popnms)
      #rm(pwGAll)
      # D
      op$pairwise$D <- op$pairwise$gst
      op$pairwise$D[lower.tri(op$pairwise$D)] <- round(pwDall, 4)
      dimnames(op$pairwise$D) <- list(popnms, popnms)
      #rm(pwDall)
      # Fst
      op$pairwise$Fst <- op$pairwise$gst
      op$pairwise$Fst[lower.tri(op$pairwise$Fst)] <- round(pwFstAll, 4)
      dimnames(op$pairwise$Fst) <- list(popnms, popnms)
      #rm(pwFstAll)
    } else {
      # gst
      op$pairwise <- list(gst = matrix(NA, nrow = np, ncol = np))
      op$pairwise$gst[lower.tri(op$pairwise$gst)] <- round(pwgAll, 4)
      dimnames(op$pairwise$gst) <- list(popnms, popnms)
      #rm(pwgAll)
      # Gst
      op$pairwise$Gst <- op$pairwise$gst
      op$pairwise$Gst[lower.tri(op$pairwise$Gst)] <- round(pwGAll, 4)
      dimnames(op$pairwise$Gst) <- list(popnms, popnms)
      #rm(pwGAll)
      # G''st
      op$pairwise$GGst <- op$pairwise$gst
      op$pairwise$GGst[lower.tri(op$pairwise$GGst)] <- round(pwGGAll, 4)
      dimnames(op$pairwise$GGst) <- list(popnms, popnms)
      # D
      op$pairwise$D <- op$pairwise$gst
      op$pairwise$D[lower.tri(op$pairwise$D)] <- round(pwDall, 4)
      dimnames(op$pairwise$D) <- list(popnms, popnms)
      #rm(pwDall)
    }
  }
  #################################
  # END
  #################################
  
  
  #################################
  # Pairwise CI
  #################################
  if(bs_pairwise){
    if(fst){
      # gst
      op$bs_pairwise <- list(gst = data.frame(populations = pwpops,
                                              actual = round(pwgAll, 4),
                                              lower = round(pwgAllLCI, 4),
                                              upper = round(pwgAllUCI, 4),
                                              row.names = NULL))
      rm(pwgAll, pwgAllLCI, pwgAllUCI)
      # Gst
      op$bs_pairwise$Gst <- data.frame(populations = pwpops,
                                       actual = round(pwGAll, 4),
                                       lower = round(pwGAllLCI, 4),
                                       upper = round(pwGAllUCI, 4),
                                       row.names = NULL)
      rm(pwGAll, pwGAllLCI, pwGAllUCI)
      # G''st
      op$bs_pairwise$GGst <- data.frame(populations = pwpops,
                                        actual = round(pwGGAll, 4),
                                        lower = round(pwGGAllLCI, 4),
                                        upper = round(pwGGAllUCI, 4),
                                        row.names = NULL)
      rm(pwGGAll, pwGGAllLCI, pwGGAllUCI)
      # D
      op$bs_pairwise$D <- data.frame(populations = pwpops,
                                     actual = round(pwDall, 4),
                                     lower = round(pwDAllLCI, 4),
                                     upper = round(pwDAllUCI, 4),
                                     row.names = NULL)
      rm(pwDall, pwDAllLCI, pwDAllUCI)
      # Fst
      op$bs_pairwise$Fst <- data.frame(populations = pwpops,
                                       actual = round(pwFstAll, 4),
                                       lower = round(pwFstAllLCI, 4),
                                       upper = round(pwFstAllUCI, 4),
                                       row.names = NULL)
      rm(pwFstAll, pwFstAllLCI, pwFstAllUCI)
    } else {
      # gst
      op$bs_pairwise <- list(gst = data.frame(populations = pwpops,
                                              actual = round(pwgAll, 4),
                                              lower = round(pwgAllLCI, 4),
                                              upper = round(pwgAllUCI, 4),
                                              row.names = NULL))
      rm(pwgAll, pwgAllLCI, pwgAllUCI)
      # Gst
      op$bs_pairwise$Gst <- data.frame(populations = pwpops,
                                       actual = round(pwGAll, 4),
                                       lower = round(pwGAllLCI, 4),
                                       upper = round(pwGAllUCI, 4),
                                       row.names = NULL)
      rm(pwGAll, pwGAllLCI, pwGAllUCI)
      # G''st
      op$bs_pairwise$GGst <- data.frame(populations = pwpops,
                                        actual = round(pwGGAll, 4),
                                        lower = round(pwGGAllLCI, 4),
                                        upper = round(pwGGAllUCI, 4),
                                        row.names = NULL)
      rm(pwGGAll, pwGGAllLCI, pwGGAllUCI)
      # D
      op$bs_pairwise$D <- data.frame(populations = pwpops,
                                     actual = round(pwDall, 4),
                                     lower = round(pwDAllLCI, 4),
                                     upper = round(pwDAllUCI, 4),
                                     row.names = NULL)
      rm(pwDall, pwDAllLCI, pwDAllUCI)
      z <- gc(reset = TRUE)
    }
  }
  #################################
  # END
  #################################
  
  
  #################################
  # Write results to file
  #################################
  if(!is.null(outfile)){
    # set up an output folder
    opf <- paste(getwd(), "/", outfile, "-[diffCalc]/", sep = "")
    dir.create(opf, showWarnings = FALSE)
    
    outnms <- names(op)
    # define a write function
    out <- sapply(outnms, function(x){
      if(x == "std_stats" || x == "global_bs"){
        ot <- paste(colnames(op[x][[1]]), collapse = "\t")
        preot <- apply(op[x][[1]], 1, paste, collapse = "\t")
        ot <- c(ot, preot)
        #fl <- file(paste(x, ".txt", sep = ""), "w")
        #cat(ot, sep = "\n", file = fl)
        #close(fl)
        writeLines(paste(ot, collapse = "\n"),
                   paste(opf, x, ".txt", sep = ""))
        ot <- NULL
      } else if(x == "pairwise"){
        statnms <- names(op[x][[1]])
        ot <- lapply(statnms, function(y){
          dat <- op[x][[1]][y][[1]]
          dat[is.na(dat)] <- ""
          dimnames(dat) <- list(NULL, NULL)
          opt <- apply(dat, 1, paste0, collapse = "\t", na.rm = "")
          popnmsOut <- paste(popnms, "\t", sep = "")
          opt <- mapply(paste, popnmsOut, opt,
                        MoreArgs = list(collapse = "\t"))
          opt <- c(y, "", paste("pops", paste(popnms, collapse = "\t"),
                                sep = "\t"), opt, "")
          return(opt)
        })
        if(fst){
          ot <- c("Pairwise stats", "gst = Nei & Chesser, (1983)",
                  "Gst = Hedrick, (2005)", "GGst = Meirmans & Hedrick, (2011)",
                  "D = Jost, (2008)",
                  "Fst = Weir & Cockerham's theta, (1984)",
                  "", unlist(ot))
        } else {
          ot <- c("Pairwise stats",
                  "gst = Nei & Chesser, (1983)",
                  "Gst = Hedrick, (2005)",
                  "GGst = Meirmans & Hedrick, (2011)",
                  "D = Jost, (2008)",
                  "", unlist(ot))
        }
        writeLines(paste(ot, collapse = "\n"),
                   paste(opf, x, ".txt", sep = ""))
        ot <- NULL
      } else if(x == "bs_locus"){
        statnms <- names(op[x][[1]])
        ot <- lapply(statnms, function(y){
          ot1 <- c("", y, "", paste(colnames(op[x][[1]][y][[1]]),
                                    collapse = "\t"))
          ot2 <- apply(op[x][[1]][y][[1]], 1, paste, collapse = "\t")
          return(c(ot1, ot2))
        })
        if(fst){
          ot <- c("Locus 95% CIs", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "GGst = Meirmans & Hedrick, (2011)",
                  "D = Jost, 2008",
                  "Fst = Weir & Cockerham, 1984",
                  "Fis = Weir & Cockerham, 1984",
                  "Fit = Weir & Cockerham, 1984", "",
                  unlist(ot))
        } else {
          ot <- c("Locus 95% CIs", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "GGst = Meirmans & Hedrick, (2011)",
                  "D = Jost, 2008",
                  unlist(ot))
        }
        writeLines(paste(ot, collapse = "\n"),
                   paste(opf, x, ".txt", sep = ""))
        ot <- NULL
      } else if(x == "pw_locus"){
        statnms <- names(op[x][[1]])
        ot <- lapply(statnms, function(y){
          ot1 <- c("", y, "", paste("Loci", paste(colnames(op[x][[1]][y][[1]]),
                                                  collapse = "\t"), sep = "\t"))
          ot2 <- apply(op[x][[1]][y][[1]], 1, paste, collapse = "\t")
          ot2 <- mapply(`paste`, rownames(op[x][[1]][y][[1]]), ot2,
                        MoreArgs = list(sep = "\t"))
          return(c(ot1, ot2))
        })
        if(fst){
          ot <- c("Locus Pairwise estimates", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "GGst = Meirmans & Hedrick, (2011)",
                  "D = Jost, 2008",
                  "Fst = Weir & Cockerham, 1984",
                  unlist(ot))
        } else {
          ot <- c("Locus Pairwise estimates", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "GGst = Meirmans & Hedrick, (2011)",
                  "D = Jost, 2008",
                  unlist(ot))
        }
        writeLines(paste(ot, collapse = "\n"),
                   paste(opf, x, ".txt", sep = ""))
        ot <- NULL
      } else if(x == "bs_pairwise"){
        statnms <- names(op[x][[1]])
        ot <- lapply(statnms, function(y){
          ot1 <- c("", y, "", paste(colnames(op[x][[1]][y][[1]]),
                                    collapse = "\t"))
          ot2 <- apply(op[x][[1]][y][[1]], 1, paste, collapse = "\t")
          return(c(ot1, ot2))
        })
        if(fst){
          ot <- c("Pairwise 95% CIs", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "GGst = Meirmans & Hedrick, (2011)",
                  "D = Jost, 2008",
                  "Fst = Weir & Cockerham, 1984",
                  unlist(ot))
        } else {
          ot <- c("Pairwise 95% CIs", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "GGst = Meirmans & Hedrick, (2011)",
                  "D = Jost, 2008",
                  unlist(ot))
        }
        writeLines(paste(ot, collapse = "\n"),
                   paste(opf, x, ".txt", sep = ""))
        ot <- NULL
      }
    })
    rm(out)
    #z <- gc(reset = TRUE)
  }
  #################################
  # END
  #################################
  
  return(op)
}
################################################################################
# end diffCalc function                                                        #
################################################################################