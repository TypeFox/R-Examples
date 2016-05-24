################################################################################
# New inCalc                                                                   #
################################################################################
#
#
#
#
#' New inCalc function for diveRsity package
#' 
#' Kevin Keenan 2014
#' 
# infile <- "Test.txt"
# boots = 10
# pairwise = FALSE
# outfile <- "out"
# para = TRUE
# xlsx <- FALSE
inCalc <- function(infile = NULL, outfile = NULL, pairwise = FALSE, 
                   xlsx = FALSE, boots = NULL, para = FALSE){
  #source("rgp.R")
  #source("inFunc.R")
  # calculate basic information
  alf <- rgp(infile)
  locs <- alf$locs
  pnms <- sapply(alf$indnms, function(x){return(x[1])})
  af <- lapply(alf$af, function(x){
    colnames(x) <- pnms
    return(x)
  })
  names(af) <- locs
  if(pairwise){
    inStatPW <- lapply(af, inFunc, pw = TRUE) 
  }
  inStatGLB <- lapply(af, inFunc, pw = FALSE)
  ####-- Global bootstrap --####
  if(!is.null(boots)){
    # population sizes
    ps <- sapply(alf$genos, function(x){
      return(dim(x)[1])
    })
    # generate resample indexes
    idx <- lapply(1:boots, function(x){
      lapply(ps, function(y){
        sample(y, y, replace = TRUE)
      })
    })
    # define a bootstrap function
    paraFunc <- function(genos, idx, af){
      # generate resamples
      rsFun <- function(x, y){
        return(x[y,,])
      }
      rsDat <- mapply(rsFun, x = genos, y = idx, SIMPLIFY = FALSE)
      
      # calculate allele frequecies
      alf <- lapply(rsDat, function(x){
        apply(x, 2, function(y){
          table(c(y[,1], y[,2]))/(length(na.omit(y[,1]))*2)
        })
      })
      # count loci
      nloci <- dim(rsDat[[1]])[2]
      # sort frequencies by loci
      locAl <- lapply(1:nloci, function(i){
        lapply(alf, "[[", i)
      })
      
      alSort <- function(x, y){
        idx <- lapply(x, function(z){
          match(names(z), rownames(y))
        })
        for(i in 1:length(idx)){
          y[idx[[i]], i] <- x[[i]]
        }
        return(y)
      }
      
      # generate allele frequency output
      alOut <- mapply(alSort, x = locAl, y = af, SIMPLIFY = FALSE)
      return(alOut)
    }
    
    if(para){
      
      cl <- parallel::makeCluster(detectCores())
      parallel::clusterExport(cl, c("inFunc", "paraFunc", "alf", "pairwise"), 
                              envir = environment())
      glbbs <- parallel::parLapply(cl, idx, function(x){
        af <- paraFunc(alf$genos, x, alf$af)
        return(lapply(af, inFunc, pw = FALSE))
      })
      if(!pairwise){
        parallel::stopCluster(cl)
      }
    } else {
      glbbs <- lapply(idx, function(x){
        af <- paraFunc(alf$genos, x, alf$af)
        return(lapply(af, inFunc, pw = FALSE))
      })
    }
    ####-- calculate global 95% CIs --####
    glbbs <- sapply(glbbs, function(x){
      unlist(x)
    })
    glbLowCI <- apply(glbbs, 1, quantile, probs = 0.025)
    glbUpCI <- apply(glbbs, 1, quantile, probs = 0.975) 
  }
  # organise data
  # global In
  glbDat <- data.frame(Locus = locs, Global_In = round(unlist(inStatGLB), 4))
  rownames(glbDat) <- NULL
  if(!is.null(boots)){
    glbDat$lower_ci <- round(glbLowCI, 4)
    glbDat$upper_ci <- round(glbUpCI, 4)
  }
  ####-- Pairwise data --####
  if(pairwise){
    # pairwise In
    pw <- combn(ncol(af[[1]]), 2)
    pwmat <- round(do.call("rbind", inStatPW), 4)
    opHeader <- paste(colnames(af[[1]])[pw[1,]], " vs ", 
                      colnames(af[[1]])[pw[2,]], sep = "")
    colnames(pwmat) <- opHeader
    
    pwDat <- data.frame(Locus = locs)
    pwDat <- cbind(pwDat, pwmat)
    rownames(pwDat) <- NULL
    row1 <- paste(colnames(pwDat), collapse = "\t")
  }
  ### Bootstrap code ###
  if(!is.null(boots) && pairwise){
    if(para){
      inbs <- parallel::parLapply(cl, idx, function(x){
        af <- paraFunc(alf$genos, x, alf$af)
        return(lapply(af, inFunc, pw = pairwise))
      })
      parallel::stopCluster(cl)
    } else {
      inbs <- lapply(idx, function(x){
        af <- paraFunc(alf$genos, x, alf$af)
        return(lapply(af, inFunc, pw = pairwise))
      })
    }
    ####-- Calculate pairwise CIs --####
    inbs <- sapply(inbs, function(x){
      return(do.call("rbind", x))
    }, simplify = "array")
    # calculate CIs
    lowCI <- round(as.data.frame(apply(inbs, c(1,2), quantile, probs = 0.025)),
                   4)
    upCI <- round(as.data.frame(apply(inbs, c(1,2), quantile, probs = 0.975)),
                  4)
    # add names
    lowCI <- cbind(glbDat$Locus, lowCI)
    upCI <- cbind(glbDat$Locus, upCI)
    rownames(lowCI) <- NULL
    colnames(lowCI) <- c("Locus", opHeader)
    rownames(upCI) <- NULL
    colnames(upCI) <- c("Locus", opHeader)
  }
  
  
  ####-- Write the data to file --####
  
  # Set up directory
  if(!is.null(outfile)){
    outDir <- paste(getwd(), "/", outfile, "-[inCalc]/", sep = "")
    if(!file.exists(outDir)){
      dir.create(path = outDir, showWarnings = FALSE)
    }
    if(xlsx){
      if ("xlsx" %in% rownames(installed.packages()) == FALSE) {
        stop("You must install the 'xlsx' package to write results in this format.")
      } else {
        outF <- paste(outDir, "outfile-[in.Calc].xlsx", sep = "")
        xlsx::write.xlsx(glbDat, file = outF, sheetName = "Global In", 
                         col.names = TRUE, row.names = FALSE, append = FALSE)
        if(pairwise){
          xlsx::write.xlsx(pwDat, file = outF, sheetName = "Pairwise In", 
                           col.names = TRUE, row.names = FALSE, append = TRUE)
        }
        if(!is.null(boots)){
          xlsx::write.xlsx(lowCI, file = outF, sheetName = "Lower CI (PW)", 
                           col.names = TRUE, row.names = FALSE, append = TRUE)
          xlsx::write.xlsx(upCI, file = outF, sheetName = "Upper CI (PW)", 
                           col.names = TRUE, row.names = FALSE, append = TRUE)
        }
      }
    } else {
      rn <- paste(colnames(glbDat), collapse = "\t")
      glbDatOut <- apply(glbDat, 1, paste, collapse = "\t")
      glbDatOut <- c(rn, glbDatOut)
      of1 <- paste(outDir, "Global-[in.Calc].txt", sep = "")
      fl1 <- file(of1, "w")
      for(i in 1:length(glbDatOut)){
        cat(glbDatOut[i], sep = "\n", file = fl1)
      }
      close(fl1)
      if(pairwise){
        pwDatOut <- apply(pwDat, 1, paste, collapse = "\t")
        pwDatOut <- c(row1, pwDatOut)
        of2 <- paste(outDir, "pairwise-[in.Calc].txt", sep = "")
        fl2 <- file(of2, "w")
        for(i in 1:length(pwDatOut)){
          cat(pwDatOut[i], sep = "\n", file = fl2)
        }
        close(fl2)
      }
      if(!is.null(boots) && pairwise){
        of3 <- paste(outDir, "Lower_CI-[in.Calc].txt", sep = "")
        fl3 <- file(of3, "w")
        of4 <- paste(outDir, "Upper_CI-[in.Calc].txt", sep = "")
        fl4 <- file(of4, "w")
        lowCIout <- c(row1, apply(lowCI, 1, paste, collapse = "\t"))
        upCIout <- c(row1, apply(upCI, 1, paste, collapse = "\t"))
        for(i in 1:length(lowCIout)){
          cat(lowCIout[i], sep = "\n", file = fl3)
          cat(upCIout[i], sep = "\n", file = fl4)
        }
        close(fl3)
        close(fl4)
      }
    }
  }
  ####-- Function outputs --####
  if(!pairwise && is.null(boots)){
    list(global = glbDat)
  } else if(pairwise && is.null(boots)){
    list(global = glbDat,
         pairwise = pwDat)
  } else if (pairwise && !is.null(boots)){
    list(global = glbDat,
         pairwise = pwDat,
         lower_CI = lowCI,
         upper_CI = upCI)    
  } else {
    list(global = glbDat)
  }
}
#' actual inCalc sub-funciton
#' 
#' Kevin Keenan 2014
inFunc <- function(af, pw = FALSE){
  if(pw){
    combs <- combn(ncol(af), 2)
    afComb <- lapply(1:ncol(combs), function(i){
      return(af[,combs[,i]])
    })
    np <- length(unique(c(combs[1,], combs[2,])))
    out <- apply(combs, 2, function(x){
      y <- af[,x]
      y <- y[rowSums(y) != 0,]
      inOut <- apply(y, 1, function(z){
        p_j <- sum(z)/length(z)
        trm1 <- (-p_j*(log(p_j)))
        lgx <- log(z)
        lgx[is.infinite(lgx)] <- 0
        trm2 <- sum(z/(length(z))*lgx) 
        return(trm1 + trm2)
      })
      return(sum(inOut))
    })
    return(out)
  } else {
    inOut <- apply(af, 1, function(x){
      p_j <- sum(x)/length(x)
      trm1 <- (-p_j*(log(p_j)))
      lgx <- log(x)
      lgx[is.infinite(lgx)] <- 0
      trm2 <- sum(x/(length(x))*lgx)  
      return(trm1 + trm2)
    })
    inOut <- sum(inOut)
  }
  return(inOut)
}