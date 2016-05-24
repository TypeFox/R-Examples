c##############################################################################
#' Calculate basic stats
#' 
#' Kevin Keenan, 2015
#' 
#' This revamped version of divBasic now allow for more rigorous tests of HWP following recommendations by Waples, 2014 J.ofHeredity, and the methods of Engels, 2009, Genetics. HWP test are carried out using the HWxtest package by Engels.

##############################################################################
#' @export
#' HWE exact testing added 29/10/2014
basicStats <- function (infile = NULL, outfile = NULL, fis_ci = FALSE,
                        ar_ci = FALSE, fis_boots = NULL, ar_boots = NULL, 
                        mc_reps = 9999, rarefaction = TRUE, 
                        ar_alpha = 0.05, fis_alpha = 0.05) {
  if("HWxtest" %in% rownames(installed.packages()) == FALSE) {
    stop("Please install HWxtest: devtools::install_github('wrengels/HWxtest', 
         subdir='pkg')", call. = FALSE)
  }
  #   infile <- "SNP_test.gen"
  #   outfile <- "NULL"
  #   HWPexact <- TRUE
  #   fis_ci = TRUE
  #   fis_alpha = 0.05
  #   ar_ci = TRUE
  #   ar_alpha = 0.05
  #   mc_reps = 9999
  #   fis_boots = 999
  #   ar_boots = 999
  #   rarefaction = FALSE
  #   para = FALSE
  #   rgp <- diveRsity:::rgp
  #   source("rarefactor.R")
  #   source("arSample.R")
  #   Rcpp::sourceCpp("obsHet.cpp")
  #   Rcpp::sourceCpp("expHet.cpp")
  #   Rcpp::sourceCpp("bsHetCalc.cpp")
  #   Rcpp::sourceCpp("hweTab.cpp")
  
  dat <- rgp(infile)
  gp <- dat$gp
  nl <- length(dat$locs)
  np <- ncol(dat$af[[1]])
  ps_tot <- sapply(dat$genos, function(x) dim(x)[1])
  
  # calculate allelic richness
  if(rarefaction){
    if(ar_ci){
      warning("Rarefaction method used to calculate allelic richness no CIs returned")
    }
    out <- list(ar = rarefactor(dat))
    mn_ar <- colMeans(out$ar[,-1], na.rm = TRUE)
    sd_ar <- apply(out$ar[,-1], 2, sd, na.rm = TRUE)
  } else {
    res <- arSample(dat = dat, nrep = ar_boots, alpha = ar_alpha, ci = ar_ci)
    out <- list(ar = res$ar)
    mn_ar <- colMeans(out$ar[,-1], na.rm = TRUE)
    sd_ar <- apply(out$ar[,-1], 2, sd, na.rm = TRUE)
  }
  
  # loci_pop_sizes
  
  out$lps = data.frame(locus = dat$locs, do.call("rbind", dat$ps))
  colnames(out$lps) <- c("Locus", gsub(",", "", sapply(dat$indnms, "[", 1)))
  
  # observed heterozygosity
  out$obs_het <- sapply(dat$genos, function(x){
    apply(x, 2, function(y){
      obsHet(y)
    })
  })
  out$obs_het <- data.frame(locus = dat$locs, out$obs_het)
  colnames(out$obs_het) <- c("Locus", gsub(",", "", sapply(dat$indnms, "[", 1)))
  
  # expected heterozygosity
  out$exp_het <- lapply(dat$af, expHet)
  out$exp_het <- do.call(rbind, out$exp_het)
  out$exp_het <- data.frame(locus = dat$locs, out$exp_het)
  colnames(out$exp_het) <- c("Locus", gsub(",", "", sapply(dat$indnms, "[", 1)))
  
  # unbiased expected heterozygosity
  out$uexp_het <- (2*out$lps[,-1]/(2*out$lps[,-1] - 1)) * out$exp_het[,-1]
  out$uexp_het <- data.frame(locus = dat$locs, round(out$uexp_het, 3))
  colnames(out$uexp_het) <- c("Locus", 
                              gsub(",", "", sapply(dat$indnms, "[", 1)))
  
  # Fis
  out$fis <- 1 - (out$obs_het[,-1]/out$exp_het[,-1])
  out$fis <- data.frame(locus = dat$locs, out$fis)
  
  # Fis bootstrapping
  if(!is.null(fis_boots) && fis_ci) {
    # create resample indexes
    idx <- lapply(1:fis_boots, function(i){
      lapply(ps_tot, function(x) sample(x, x, TRUE))
    })
    # fis wrapper
    fiscalc <- function(genos, idx){
      ho <- apply(genos[idx,,], 2, obsHet)
      he <- apply(genos[idx,,], 2, bsHetCalc)
      return(1-(ho/he))
    }
    
    fis_bs <- simplify2array(lapply(idx, function(x){
      mapply(fiscalc, genos = dat$genos, idx = x, SIMPLIFY = TRUE)
    }))
    
    loc_cis <- apply(fis_bs, c(1,2), quantile, 
                     prob = c(fis_alpha/2, 1 - (fis_alpha/2)), 
                     na.rm = TRUE)
    glb_cis <- apply(fis_bs, c(2,3), mean, na.rm = TRUE)
    glb_cis <- apply(glb_cis, 1, quantile, 
                     prob = c(fis_alpha/2, 1 - (fis_alpha/2)), na.rm = TRUE)
  }
  
  # HWP testing
  # exact
  hwe_res <- lapply(dat$genos, function(y){
    apply(y, 2, function(x){
      tx <- hweTab(x)
      if(length(tx) == 0L){
        Pvalues = rep(NA, 4)
        names(Pvalues) <- c("LLR", "Prob", "U", "Chisq")
        list(Pvalues = Pvalues, observed = Pvalues, method = NA, statName = NA)
      } else {
        HWxtest::hwx.test(tx, detail = 0, method = "auto", B = mc_reps)
      }
    })
  })
  hwe_p <- lapply(hwe_res, function(x){
    lapply(x, "[[", "Pvalues")
  })
  # Likelihood ratio
  out$hwe_llr_p <- t(do.call(rbind, lapply(hwe_p, function(x){
    sapply(x, "[", "LLR")
  })))
  rownames(out$hwe_llr_p) <- NULL
  out$hwe_llr_p <- data.frame(Locus = dat$locs, out$hwe_llr_p)
  colnames(out$hwe_llr_p) <- c("Locus", 
                               gsub(",", "", sapply(dat$indnms, "[", 1)))
  
  # Directional HWE
  # wrapper to only return U for homozygosity excess tests
  homU <- function(loc){
    if(is.na(loc$observed["U"])){
      return(NA)
    } else if(loc$observed["U"] <= 0L){
      return(loc$Pvalues["U"])
    } else {
      return(NA)
    }
  }
  # extract pvalues for homoygosity excess
  out$hwe_hom <- t(do.call(rbind, lapply(hwe_res, function(x){
    sapply(x, homU)
  })))
  rownames(out$hwe_hom) <- NULL
  out$hwe_hom <- data.frame(Locus = dat$locs, out$hwe_hom)
  colnames(out$hwe_hom) <- c("Locus", 
                             gsub(",", "", sapply(dat$indnms, "[", 1)))
  
  # wrapper to only return U for heterozygosity excess
  hetU <- function(loc){
    if(is.na(loc$observed["U"])){
      return(NA)
    } else if(loc$observed["U"] > 0L){
      return(loc$Pvalues["U"])
    } else {
      return(NA)
    }
  }
  # extract pvalues for heterozygosity excess
  out$hwe_het <- t(do.call(rbind, lapply(hwe_res, function(x){
    sapply(x, hetU)
  })))
  rownames(out$hwe_het) <- NULL
  out$hwe_het <- data.frame(Locus = dat$locs, out$hwe_het)
  colnames(out$hwe_het) <- c("Locus", 
                             gsub(",", "", sapply(dat$indnms, "[", 1)))
  
  # Create the main table
  out$main_tab <- lapply(2:(np+1), function(i){
    data.frame(t(sapply(out, function(x){
      return(x[,i])
    })))
  })
  # add locus names (rounding is really slow for large data)
  out$main_tab <- lapply(out$main_tab, function(x){
    colnames(x) <- dat$locs
  #  x <- round(x, 3)
    return(x)
  })
  # add multi-locus estimates
  out$main_tab <- lapply(out$main_tab, function(x){
    mns <- apply(x[1:6, ], 1, mean, na.rm = TRUE)
    pvals <- apply(x[7:9, ], 1, function(y){
      df <- 2 * length(y)
      chi <- -2 * sum(log(y), na.rm = TRUE)
      p <- pchisq(chi, df = df, lower.tail = FALSE)
    })
    x$overall <- round(c(mns, pvals), 3)
    return(x)
  })
  
  # add cis to main tab
  if(fis_ci){
    out$main_tab <- lapply(1:np, function(x){
      x <- rbind(out$main_tab[[x]], round(c(loc_cis[1, ,x], glb_cis[1,x]), 3), 
                 round(c(loc_cis[2, ,x], glb_cis[2,x]), 3))
      rownames(x)[c(nrow(x)-1, nrow(x))] <- c("fis_lo", "fis_hi")
      return(x)
    })
  }
  if(ar_ci && !rarefaction){
    out$main_tab <- lapply(1:np, function(x){
      x <- rbind(out$main_tab[[x]], 
                 round(as.vector(c(res$lo[,x+1], res$mean_ci[[x]][1])), 3),
                 round(as.vector(c(res$hi[,x+1], res$mean_ci[[x]][2])), 3))
      rownames(x)[c(nrow(x)-1, nrow(x))] <- c("ar_lo", "ar_hi")
      return(x)
    })
  } else if(ar_ci && rarefaction) {
    warning("Rarefaction method used to calculate allelic richness, no CIs returned")
  }
  
  out$main_tab <- lapply(out$main_tab, function(x){
    rownames(x)[2] <- "size"
    rownames(x)[which(rownames(x) == "hwe_llr_p")] <- "hwe_glb"
    return(x)
  })
  
  names(out$main_tab) <- gsub(",", "", sapply(dat$indnms, "[", 1))
  
  if(!is.null(outfile)){
    main_tab <- lapply(1:np, function(i){
      x <- apply(out$main_tab[[i]], 1, paste, collapse = "\t")
      x <- mapply(paste, rownames(out$main_tab[[i]]), x, 
                  MoreArgs = list(sep = "\t"))
      x <- c("", gsub(",", "", dat$indnms[[i]][1]), 
             paste("stat", paste(colnames(out$main_tab[[i]]), collapse = "\t"),
                   sep = "\t"), x)
      return(x)
    })
    main_tab <- paste(unlist(main_tab), collapse = "\n")
    writeLines(main_tab, paste(outfile, ".txt", sep = ""))
  }
  return(out)
}