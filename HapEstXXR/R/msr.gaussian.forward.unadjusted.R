msr.gaussian.forward.unadjusted <-
  function(
    snps, trait, lim, maxSNP, 
    nt, sort.by, selection, p.threshold,
    pair.begin, pattern.begin.mat, 
    baseline.hap, min.count) {
    
    N <- dim(snps)[1]
    ns <- dim(snps)[2]
    nt <- as.integer(nt)  
    ## Begin stepwise regression
    res <- list(NA)
    if (any(!is.na(pattern.begin.mat))) {
      if (pair.begin == TRUE) {
        stop(paste("choose either", 
                   " begin with defined pattern or pair.wise.", 
                   sep = ""))
      }
      if ((dim(pattern.begin.mat)[2] >=  maxSNP) ||
            (dim(pattern.begin.mat)[2] >=  dim(snps)[2]))
        stop("dimension of begin.pattern.mat not adequate.")
    }
    if (pair.begin == FALSE) { 
      #########################################################################
      # begin single SNP test
      cat(paste("Iteration with 1 SNP  at same time.   System.time = ", 
                Sys.time(), "\n", sep = ""))
      nind <- df <- pval <- aic <- aicc <- rep(NA, dim(snps)[2])
      single.test <- single.snp.test(snps, 
                                     trait,
                                     prt = FALSE,
                                     type = "gaussian")
      nind  <- as.integer(single.test$N)
      df    <- as.integer(rep(1, length(single.test$N)))
      pval  <- as.numeric(single.test$p.value)
      aic   <- as.numeric(single.test$AIC)
      aicc  <- as.numeric(single.test$AICc)
      i <- 1
      if (all(pval > p.threshold[1])) {
        cat("None of the single test results fulfilled the p value threshold.")        
        return(res)
      } 
      res[[1]] <- data.frame(a = as.integer(single.test$SNP), 
                             b = rep("gaussian", length(nind)), 
                             c = as.integer(nind), 
                             d = as.integer(df), 
                             e = as.numeric(pval), 
                             f = as.numeric(aic),
                             g = as.numeric(aicc),
                             stringsAsFactors = FALSE,
                             row.names = NULL)
      colnames(res[[1]]) <- c(paste("snp", 1:i, sep = ""),
                              "type", "nSubj", "df", "p.value",
                              "AIC", "AICc")
      res[[i]] <- res[[i]][pval < p.threshold[1], ]
      res[[i]] <- (res[[i]])[
        order((res[[i]])[, sort.by])[1:min(nt, nrow(res[[i]]))], ,
        drop = FALSE]
      rownames(res[[1]]) <- NULL
    } else {
      #########################################################################
      # start with all two pair haplotypes!
      # begin pair wise
      cat(paste("Iteration with 2 SNPs at same time.   System.time = ",
                Sys.time(), "\n", sep = ""))
      # construct all pairs
      Z <- 1:ns
      X <- rep(Z, rep.int(length(Z), length(Z)))
      Y <- rep(Z, times = ceiling(length(X) / length(Z)))
      cont <- ifelse(Y > X, TRUE, FALSE)
      snp.pos <- cbind(X[cont], Y[cont])
      rm(X, Y, Z)
      nind <- df <- pval <- aic <- aicc <- rep(NA, dim(snp.pos)[1])
      # Analysis of all pairs
      cat(paste("Number of SNP pairs = ", dim(snp.pos)[1], "\n\n", sep = ""))
      for(j in 1:dim(snp.pos)[1]) {
        if ((j %% 5000) == 0) {
          cat(paste("Step =  ", j,
                    "   System.time = ",
                    (Sys.time()), "\n", sep = "")) 
        }
        geno.pair <- snps[, snp.pos [j, ], drop = FALSE]
        hap.test <- msr.gaussian.haplotype.test.unadjusted(
          geno.pair, 
          trait, 
          lim  = lim, 
          baseline.hap = baseline.hap, 
          min.count = min.count)
        # found no haplotypes with probability over lim
        # p is greater then the threshold
        if (all(is.na(hap.test$haplotypes)) |
              ifelse(is.na(hap.test$global.p.value), TRUE, 
                     hap.test$global.p.value > p.threshold[2])) {
          nind[j] <- NA
          df[j]   <- NA
          pval[j] <- NA
          aic[j]  <- NA
          aicc[j] <- NA
        } else {
          nind[j]  <- hap.test$nSubj
          df[j]    <- hap.test$df
          pval[j]  <- hap.test$global.p.value
          aic[j]   <- hap.test$AIC
          aicc[j]   <- hap.test$AICc
        }
      }
      i <- 2
      if(all(is.na(nind))) {
        cat(paste("None of the SNP pairwise test results",
                  "fulfilled the p value threshold."), sep = "") 
        return(res)
      } 
      ii <- !is.na(nind)
      res[[i]] <- data.frame(snp.pos[ii, , drop = FALSE], 
                             "gaussian", 
                             nind[ii],
                             df[ii], 
                             pval[ii],
                             aic[ii],
                             aicc[ii],
                             stringsAsFactors = FALSE,
                             row.names = NULL)
      colnames(res[[i]]) <- c(paste("snp", 1:i, sep = ""),
                              "type", "nSubj", "df", "p.value",
                              "AIC", "AICc")
      res[[i]] <- (res[[i]])[
        order((res[[i]])[, sort.by])[1:min(nt, nrow(res[[i]]))], ,
        drop = FALSE]
      rownames(res[[i]]) <- NULL
    }
    if(length(res[[i]]) < 1) return(res)
    #################################################
    # Next iteration: from 2, 3, 4 ...(or col number of 
    # pattern.begin.mat) to maxSNP
    i <- i + 1
    # i No. of SNPs for every haplotype pattern
    while(i <=  maxSNP) {
      cat(paste("Iteration with ", i, 
                " SNPs at same time.   System.time = ",
                (Sys.time()), "\n", sep = ""))
      # prepare selection
      sel.thres <- 0
      if (selection == 1) {
        sel.thres <- min(as.numeric(res[[i - 1]][ , "AICc"]))
      } else {
        if (selection == 2) {
          sel.thres <- min(as.numeric(res[[i - 1]][ , "AIC"]))
        } else {
          if (selection == 3) {
            sel.thres <- min(as.numeric(res[[i - 1]][ , "p.value"]))
          } else {
            if (selection == 4) {
              nsel <- dim(res[[i - 1]])[1]
              sel.thres <- mean(log10(sort(
                as.numeric((res[[i - 1]])[ , "p.value"])[1:(min(nsel, 10))])))
            }
          }
        }
      }
      # create matrix with all possible combinations
      snp.pos <- as.matrix(res[[i - 1]] [, 1:(i - 1), drop = FALSE])
      storage.mode(snp.pos) <- "integer"
      newdim <- as.integer(c(dim(snp.pos)[1] * (ns - 1), dim(snp.pos)[2] + 1))
      out <- .C("create_pattern_matrix", 
                pattern   = as.integer(snp.pos),
                ndim      = dim(snp.pos), 
                snps      = as.integer(1:ns), 
                snplen    = ns, 
                newpat    = as.integer(rep(0, newdim[1] * newdim[2])), 
                newpatdim = newdim,
                len       = as.integer(0))
      snp.pos <- (matrix(out$newpat, 
                         nrow = newdim[1],
                         ncol = newdim[2], 
                         byrow = F))[1:out$len, , drop = FALSE]
      cat(paste("Iteration with ", i, 
                " SNPs at same time. ---- ", 
                dim(snp.pos)[1], 
                " detected SNP combinations -----\n",
                sep = ""))
      nind <- df <- pval <- aic <- aicc <- snp.sel <- NULL 
      # evaluate all combinations in the step before
      k <- 0
      for(j in 1:(dim(snp.pos)[1])) {
        Pos <- as.integer(snp.pos[j, ])
        if ((j %% 5000) == 0) { 
          cat(paste("Step =  ", j, 
                    "   System.time = ", 
                    (Sys.time()), "\n", sep = "")) 
        }
        geno <- matrix(snps[, Pos], N, i) 
        hap.test <- msr.gaussian.haplotype.test.unadjusted(
          geno,
          trait,
          lim  = lim, 
          baseline.hap = baseline.hap, 
          min.count = min.count)
        if (all(is.na(hap.test$haplotypes))) {
          # found no haplotypes with probability over lim 
          txt.pos <- paste(Pos, collapse = " ")
          cat("Step ", i, 
              ": SNPs ", txt.pos, 
              " all inferred haplotypes with probability ", 
              "below threshold(lim = ", lim, ")\n", sep = "")
        } else {
          save.yes <- FALSE
          if (hap.test$global.p.value < p.threshold[i]) {     
            if (selection == 0 | selection == 4) {
              save.yes <- TRUE
            } else {
              if (selection == 1 & (as.numeric(hap.test$AICc) < sel.thres)) {
                save.yes <- TRUE
              } else {
                if (selection == 2 & (hap.test$AIC < sel.thres)) {
                  save.yes <- TRUE
                } else {
                  if (selection == 3 & (hap.test$global.p.value < sel.thres)) {
                    save.yes <- TRUE
                  }
                }
              }
            }
          }
        } 
        if (save.yes) {
          snp.sel <- rbind(snp.sel, Pos)
          nind <- c(nind, hap.test$nSubj) 
          df   <- c(df, hap.test$df)
          pval <- c(pval, hap.test$global.p.value)
          aic  <- c(aic, hap.test$AIC)
          aicc <- c(aicc, hap.test$AICc)
        }       
      } # for j
      # save results
      ntest.i <- length(nind)
      if (length(nind) < 1) {
        return(res)
      }
      if (selection == 0) {
        ii <- 1:length(nind)
      } else {
        if (selection == 1) {
          ii <- order(aicc)[1:(min(nt, ntest.i))]         
        } else {
          if (selection == 2) {
            ii <- order(aic)[1:(min(nt, ntest.i))]         
          } else {
            if (selection == 3) {
              ii <- order(pval)[1:(min(nt, ntest.i))]
            } else {
              if (selection == 4) {
                mean.log10.p10 <- mean(log10(
                  as.numeric(pval)[1:(min(ntest.i, 10))]))
                if (mean.log10.p10 < (sel.thres * 1.10)) {
                  ii <- 1:ntest.i
                } else {
                  return(res)
                }
              }
            }
          }
        }
      }
      res[[i]] <- data.frame(snp.sel[ii, , drop = FALSE], 
                             "gaussian", 
                             nind[ii], 
                             df[ii], 
                             pval[ii],
                             aic[ii],
                             aicc[ii],
                             stringsAsFactors = FALSE,
                             row.names = NULL)
      colnames(res[[i]]) <- c(paste("snp", 1:i, sep = ""),
                              "type", "nSubj", "df", "p.value",
                              "AIC", "AICc")
      res[[i]] <- (res[[i]])[order((res[[i]])[, sort.by])
                             [1:min(nt, ntest.i)], , drop = FALSE]
      rownames(res[[i]]) <- 1:dim(res[[i]])[1]
      i <- i + 1
    } # while i
    return(res)
  }