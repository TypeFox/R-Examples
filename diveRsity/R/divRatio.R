################################################################################
# divRatio: calculates diversity standardised to yardstick popukation
################################################################################
#' @export
divRatio <- function(infile = NULL, outfile = NULL, gp = 3, pop_stats =  NULL, 
                     refPos = NULL, boots = 1000,  para = FALSE) {
  #data(Test_data, package = "diveRsity")
  #infile <- Test_data
  #outfile = NULL
  #refPos = 5
  #gp = 3
  popStats = pop_stats
  #boots = 99#boots
  #fileReader <- diveRsity::fileReader
  #AR <- diveRsity:::AR
  #Hex <- diveRsity:::Hex
  #arHex <- diveRsity:::arHex
  #para = TRUE
  # create a directory for output
  if(!is.null(outfile)){
    suppressWarnings(dir.create(path=paste(getwd(),"/", outfile,
                                           "-[diveRsity]","/",sep="")))
    of = paste(getwd(), "/", outfile, "-[diveRsity]", "/", sep = "")
    write_res <- is.element("xlsx", installed.packages()[, 1])
  }
  # read the allelic richness and heterozygosity functions
  data1 <- fileReader(infile)
  data1[data1==0]<-NA;data1[data1=="999999"]<-NA;data1[data1=="000000"]<-NA
  #raw_data<-data1
  npops<-length(which(toupper(data1[,1]) == "POP"))
  pop_pos<- c(which(toupper(data1[,1]) == "POP"), (nrow(data1)+1))
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  # Calculate the minimum sample size
  pop_sizes <- sapply(1:npops, function(i){
    pop_pos[(i+1)] - pop_pos[i]-1
  })
  #minSize <- min(pop_sizes) 
  pop_names<-as.character(data1[(pop_pos[1:npops]+1),1])
  nloci<- (pop_pos[1]-2)
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  pop_list <- lapply(1:npops, function(i){
    return(as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                           2:(nloci+1)]))
  })  
  if (gp==3) {
    plMake<-function(x){
      out <- matrix(sprintf("%06g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "0000NA"] <- "    NA"
      }
      return(out)
    }
  } else if (gp==2) {
    plMake<-function(x){
      out <- matrix(sprintf("%04g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "00NA"] <- "  NA"
      }
      return(out)
    }
  }
  suppressWarnings(pop_list<-lapply(pop_list, plMake))
  # deal with missing data
  if (gp == 3){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "    NA"]<-NA
    }
  } else if (gp == 2){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "  NA"] <-NA
    }
  }
  ##############################################################################
  # if only the refpop raw data is given
  if(npops == 1 && !is.null(pop_stats)){
    refPop <- pop_list[[1]]
    refPos <- 1
    # read subject population stats
    trypopDF <- try(read.table(popStats, header = TRUE), silent = TRUE)
    if(is(trypopDF, "try-error")){ 
      message <- paste("[ERROR]",
                       "",
                       "There is a problem with 'pop_stats' file format",
                       "",
                       "See the package user manual for more details",
                       "of file format requirements",
                       "",
                       sep = "\n")
      cat(message)
      stop()
    } else {
      popDF <- read.table(popStats, header = TRUE)
    }
    # calculate refpop standard stats
    refalr <- rowMeans(replicate(boots, AR(refPop)))
    refhe <- Hex(refPop)
    refStd <- data.frame(alr = refalr, hexp = refhe)
    refStd <- data.frame(pops = paste(pop_names, "-(ref)", sep = ""),
                         n = nrow(refPop),
                         alr = mean(refStd$alr, na.rm = TRUE),
                         alrse = sd(refStd$alr, na.rm = TRUE) / 
                           sqrt(length(na.omit(refStd$alr))),
                         he = mean(refStd$hexp, na.rm = TRUE),
                         hese = sd(refStd$hexp, na.rm = TRUE) / 
                           sqrt(length(na.omit(refStd$hexp)))
    )
    if(is.element("validloci", names(popDF))){
      refStd$validloci <- paste(loci_names, sep = "\t", collapse = "\t")
    }
    # extract subject population sizes
    popDF <- rbind(refStd, popDF)
    popSizes <- as.numeric(popDF$n)
    # extract valid locus information
    if(is.element("validloci", names(popDF))){
      vlocs <- as.character(popDF$validloci)
      vlocs <- sapply(vlocs, function(x){
        return(strsplit(x, split = "\\s+"))
      })
      validLoc <- lapply(vlocs, function(x){
        return(sapply(x, function(y){
          as.numeric(which(loci_names == y))
        }))
      })
    } else {
      validLoc <- lapply(1:length(popSizes), function(...){
        locs <- 1:nloci
        names(locs) <- loci_names 
        return(locs)
      })
    }    
    
    # calculate refPopStats based on each subject pop sample size
    ###########################################################################
    if(para){
      
        cores <- parallel::detectCores()
        cl <- parallel::makeCluster(cores)
        parallel::clusterExport(cl, c("arHex", "boots", "refPop", "popSizes", 
                                      "gp", "nloci", "validLoc"), 
                                envir = environment())
        refPopStats <- parallel::parLapply(cl, seq_along(popSizes), function(i){
          inner <- list(ref = refPop[,validLoc[[i]]], size = popSizes[i],
                        gp = gp)
          outer <- replicate(boots, arHex(inner), simplify = FALSE)
          alrPre <- sapply(outer, function(x){
            return(x$alls)
          })
          alr <- rowMeans(alrPre)
          hexpPre <- sapply(outer, function(x){
            return(x$hexp)
          })
          hexp <- rowMeans(hexpPre)
          return(data.frame(alr = alr, hexp = hexp))
        })
        parallel::stopCluster(cl) 
      
    } else {
      refPopStats <- lapply(1:length(popSizes), function(i){
        inner <- list(ref = refPop, size = popSizes[i],
                      gp = gp)
        outer <- replicate(boots, arHex(inner), simplify = FALSE)
        alrPre <- sapply(outer, function(x){
          return(x$alls)
        })
        alr <- rowMeans(alrPre)
        hexpPre <- sapply(outer, function(x){
          return(x$hexp)
        })
        hexp <- rowMeans(hexpPre)
        return(data.frame(alr = alr, hexp = hexp))
      })
    }
    ###########################################################################
    # calculate the means and s.e. for all refPopStats
    refs <- lapply(refPopStats, function(x){
      meanAlr <- mean(x$alr, na.rm = TRUE)
      seAlr <- sd(x$alr, na.rm = TRUE)/sqrt(length(na.omit(x$alr)))
      meanHexp <- mean(x$hexp, na.rm = TRUE)
      seHexp <- sd(x$hexp, na.rm = TRUE)/sqrt(length(na.omit(x$hexp)))
      list(alr = data.frame(mean = meanAlr, se = seAlr),
           hexp = data.frame(mean = meanHexp, se = seHexp))
    })
    
    # extract all subject pops info
    subs <- lapply(1:nrow(popDF), function(i){
      return(list(alr = data.frame(mean = popDF$alr[i], se = popDF$alrse[i]),
                  hexp = data.frame(mean = popDF$he[i], se = popDF$hese[i])))
    })
    ##############################################################################
    # calculate the ratio stats
    seRatCalc <- function(subs, refs){
      rat <- subs$mean/refs$mean
      seRat <- sqrt((rat^2) * (((subs$se/subs$mean)^2) + ((refs$se/refs$mean)^2)))
      return(seRat)
    }
    tsapply <- function(...) t(sapply(...))  
    divRatio <- tsapply(1:length(subs), function(i){
      # create a subset variable for refs
      alrRat <- subs[[i]]$alr$mean / refs[[i]]$alr$mean
      hexpRat <- subs[[i]]$hexp$mean / refs[[i]]$hexp$mean
      alrSErat <- seRatCalc(subs[[i]]$alr, refs[[i]]$alr) 
      hexpSErat <- seRatCalc(subs[[i]]$hexp, refs[[i]]$hexp)
      res <- c(pop = as.character(popDF$pops[i]),
               n = popSizes[i],
               alr = round(subs[[i]]$alr$mean, 4),
               alrSE = round(subs[[i]]$alr$se, 4),
               He = round(subs[[i]]$hexp$mean, 4),
               HeSE = round(subs[[i]]$hexp$se, 4),
               alrRatio = round(alrRat, 4), 
               alrSEratio = round(alrSErat, 4),
               heRatio = round(hexpRat, 4),
               heSEratio = round(hexpSErat, 4))
      return(res)
    })  
  } else {
    ############################################################################
    # Subset pop_list into subject populations and reference population
    # reference population
    refPop <- pop_list[[refPos]]
    
    # Run AR and Hexpected function for each population other than the refpop
    # call them subject populations
    if(para){
      
        cores <- parallel::detectCores()
        cl <- parallel::makeCluster(cores)
        parallel::clusterExport(cl, c("Hex", "AR", "boots", "gp", "nloci",
                                      "pop_list"), 
                                envir = environment())
        subPopStats <- parallel::parLapply(cl, pop_list, function(x){
          # Calculate allelic richness
          # bootstrap first
          alrbs <- replicate(boots, AR(x))
          # calculate the mean of the bootstraps per locus
          alr <- rowMeans(alrbs)
          # Calculate expected Het
          hex <- Hex(x)
          # create return obj
          return(data.frame(alr = alr, hexp = hex))
        }) 
        
    } else {
      subPopStats <- lapply(pop_list, function(x){
        # Calculate allelic richness
        # bootstrap first
        alrbs <- replicate(boots, AR(x))
        # calculate the mean of the bootstraps per locus
        alr <- rowMeans(alrbs)
        # Calculate expected Het
        hex <- Hex(x)
        # create return obj
        return(data.frame(alr = alr, hexp = hex))
      })
    }
    
    # Check if any loci in each population is missing data
    validLocs <- lapply(subPopStats, function(x){
      which(!is.na(x[,1]))
    })
    
    # calculate the standardized alr and hex for the ref pop
    if(para){
      parallel::clusterExport(cl, c("arHex", "refPop", "boots", "gp", 
                                    "pop_sizes", "validLocs"), 
                              envir = environment())
      refPopStats <- parallel::parLapply(cl, seq_along(pop_list), function(i){
        inner <- list(ref = refPop[,validLocs[[i]]], size = pop_sizes[i],
                      gp = gp)
        outer <- replicate(boots, arHex(inner), simplify = FALSE)
        alrPre <- sapply(outer, function(x){
          return(x$alls)
        })
        alr <- rowMeans(alrPre)
        hexpPre <- sapply(outer, function(x){
          return(x$hexp)
        })
        hexp <- rowMeans(hexpPre)
        return(data.frame(alr = alr, hexp = hexp))
      })
      parallel::stopCluster(cl)
    } else {
      refPopStats <- lapply(seq_along(pop_list), function(i){
        inner <- list(ref = refPop, size = pop_sizes[i],
                      gp = gp)
        outer <- replicate(boots, arHex(inner), simplify = FALSE)
        alrPre <- sapply(outer, function(x){
          return(x$alls)
        })
        alr <- rowMeans(alrPre)
        hexpPre <- sapply(outer, function(x){
          return(x$hexp)
        })
        hexp <- rowMeans(hexpPre)
        return(data.frame(alr = alr, hexp = hexp))
      })
    }
    # calculate the means and s.e. for all refPopStats
    refs <- lapply(refPopStats, function(x){
      meanAlr <- mean(x$alr, na.rm = TRUE)
      seAlr <- sd(x$alr, na.rm = TRUE)/sqrt(length(na.omit(x$alr)))
      meanHexp <- mean(x$hexp, na.rm = TRUE)
      seHexp <- sd(x$hexp, na.rm = TRUE)/sqrt(length(na.omit(x$hexp)))
      list(alr = data.frame(mean = meanAlr, se = seAlr),
           hexp = data.frame(mean = meanHexp, se = seHexp))
    })
    # calculate the means and s.e. for all subPopStats
    #if()
    subs <- lapply(subPopStats, function(x){
      meanAlr <- mean(x$alr, na.rm = TRUE)
      seAlr <- sd(x$alr, na.rm = TRUE)/sqrt(length(na.omit(x$alr)))
      meanHexp <- mean(x$hexp, na.rm = TRUE)
      seHexp <- sd(x$hexp, na.rm = TRUE)/sqrt(length(na.omit(x$hexp)))
      list(alr = data.frame(mean = meanAlr, se = seAlr),
           hexp = data.frame(mean = meanHexp, se = seHexp))
    })
    # Account for refpop calculation variation
    subs[refPos] <- refs[refPos]
    ##############################################################################
    # calculate the ratio stats
    seRatCalc <- function(subs, refs){
      rat <- subs$mean/refs$mean
      seRat <- sqrt((rat^2) * (((subs$se/subs$mean)^2) + ((refs$se/refs$mean)^2)))
      return(seRat)
    }
    tsapply <- function(...) t(sapply(...))  
    divRatio <- tsapply(seq_along(pop_list), function(i){
      # create a subset variable for refs
      alrRat <- subs[[i]]$alr$mean / refs[[i]]$alr$mean
      hexpRat <- subs[[i]]$hexp$mean / refs[[i]]$hexp$mean
      alrSErat <- seRatCalc(subs[[i]]$alr, refs[[i]]$alr) 
      hexpSErat <- seRatCalc(subs[[i]]$hexp, refs[[i]]$hexp)
      res <- c(pop = pop_names[i],
               n = pop_sizes[i],
               alr = round(subs[[i]]$alr$mean, 4),
               alrSE = round(subs[[i]]$alr$se, 4),
               He = round(subs[[i]]$hexp$mean, 4),
               HeSE = round(subs[[i]]$hexp$se, 4),
               alrRatio = round(alrRat, 4), 
               alrSEratio = round(alrSErat, 4),
               heRatio = round(hexpRat, 4),
               heSEratio = round(hexpSErat, 4))
      return(res)
    })
  }
  # add reference data to divRatio
  popnms <- divRatio[,1]
  stst <- divRatio[,-1]
  class(stst) <- "numeric"
  refPop <- stst[refPos,]
  divRatio <- stst[-refPos,]
  popnms[refPos] <- paste(popnms[refPos], "-(ref)", sep = "")
  divRatio <- as.data.frame(rbind(refPop, divRatio))
  row.names(divRatio) <- NULL
  rp <- popnms[refPos]
  popnms <- c(rp, popnms[-refPos])
  divRatio <- cbind(popnms, divRatio)
  colnames(divRatio)[1] <- "pops"
  #refPop <- divRatio[refPos,]
  #divRatio <- divRatio[-refPos,]
  #refPop[1] <- paste(refPop[1], "-(ref)", sep = "")
  #divRatio <- as.data.frame(rbind(refPop, divRatio))
  if(!is.null(outfile)){
    if(write_res){
      # standard stats
      xlsx::write.xlsx(divRatio, file = paste(of, "[divRatio].xlsx",sep = ""),
                       sheetName = "Diversity_ratios", col.names = TRUE,
                       row.names = FALSE, append = FALSE)
    } else {
      write.table(divRatio, "divRatio-out.txt", col.names = TRUE, 
                  row.names = FALSE, append = FALSE, sep = "\t", 
                  quote = FALSE)
    }
  }
  divRatio[,-1] <- apply(divRatio[,-1], 2, function(x){
    return(as.numeric(as.character(x)))
  })
  divRatio[,1] <- as.character(divRatio[,1])
  return(divRatio)
}
################################################################################
# end divRatio
################################################################################
################################################################################
# AR: calculates the number of allele per locus from pop_list (divRatio)
################################################################################
# This function accepts a list containing a single population sample from
# a standard pop_list object, usually the ref sample and an interger 
# representing the resample size used to bootstrap allelic richness
# Define Allelic richness function for a single population to
# be bootstrapped for a given sample size
AR <- function(x){
  if(length(x) == 2L){
    pl <- x$ref
    mSize <- x$size
    bser<-function(x){
      return(matrix(x[sample(nrow(x), mSize, replace = TRUE), ], ncol = ncol(x)))
    }
    pop_list <- bser(pl) # resample
  } else {
    pl <- x
    mSize <- nrow(pl)
    bser<-function(x){
      return(matrix(x[sample(nrow(x), mSize, replace = TRUE), ], ncol = ncol(x)))
    }
    pop_list <- bser(pl) # resample
  }
  nloci <- ncol(pop_list)
  gp = nchar(pop_list[1,1])/2
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss <- function(x){  # where x is object pop_list
      pl <- list()
      pl[[1]] <- matrix(substr(x, 1, 2), ncol = nloci)
      pl[[2]] <- matrix(substr(x, 3, 4), ncol = nloci)
      return(pl)
    }
  }
  pop_alleles <- pl_ss(pop_list)
  alln <- function(x){
    res <- sapply(1:ncol(x[[1]]), function(i){
      list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
    })
  }
  allele_names <- alln(pop_alleles)
  Alls <- sapply(allele_names, function(x){
    length(x)
  })
  return(Alls)
}
################################################################################
# End AR
################################################################################
################################################################################
# Hex: calcuates expected heterozygosity from pop_list (divRatio)
################################################################################
# This function accepts a list containing a single population sample from
# a standard pop_list object, usually the ref sample and an interger 
# representing the resample size used to bootstrap expected heterozygosity
# Define a hexp function
Hex <- function(x){
  if(length(x) == 2L){
    pl <- x$ref
    mSize <- x$size
    bser<-function(x){
      return(matrix(x[sample(nrow(x), mSize, replace = TRUE), ],ncol=ncol(x)))
    }
    pop_list <- bser(pl) # resample
  } else {
    pop_list <- x
  }
  nloci = ncol(pop_list)
  gp = nchar(pop_list[1,1])/2
  # split string genotypes
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3), ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6), ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,2), ncol=nloci)
      pl[[2]]<-matrix(substr(x,3,4), ncol=nloci)
      return(pl)
    }
  }
  pop_alleles <- pl_ss(pop_list)
  alln <- function(x){
    res <- sapply(1:ncol(x[[1]]), function(i){
      list(sort(unique(c(x[[1]][, i], x[[2]][, i])), decreasing = FALSE))
    })
  }
  allele_names <- alln(pop_alleles)
  # Calculate expected He
  #if(npops == 1){
  #  loci_combi <- allele_names[,1]
  #} else {
  #  loci_combi <- apply(allele_names, 1, FUN = 'unlist')
  #}
  aaList <- function(x){
    return(sort(unique(x, decreasing = FALSE)))
  }
  # finter out unique alleles
  all_alleles <- lapply(allele_names, aaList)
  # Create allele frequency holders
  allele_freq <- lapply(1:ncol(pop_list), function(i){
    Nrow <- length(all_alleles[[i]])
    Ncol <- length(pop_list)
    mat <- matrix(rep(0,(Ncol * Nrow)), ncol = Ncol)
    rownames(mat) <- all_alleles[[i]]
    return(mat)
  })
  # rbind pop_alleles
  pa1 <- rbind(pop_alleles[[1]], pop_alleles[[2]])
  
  # Count alleles
  actabPre <- function(x){
    lapply(1:ncol(x), function(i){
      table(x[,i])
    })
  }
  actab <- actabPre(pa1) 
  # Count the number of individuals typed per locus per pop
  indtyppop1 <- function(x){
    apply(x, 2, function(y){
      length(na.omit(y))/2
    })
  }
  indtyppop <- indtyppop1(pa1)
  #calculate allele frequencies
  afCalcpop <- sapply(1:length(actab), function(i){
    actab[[i]]/(indtyppop[i] * 2)
  })
  # calculate heterozygosities
  Hexp <- sapply(afCalcpop, function(x){
    1 - (sum(x^2))
  })
  return(Hexp)
}
################################################################################
# End Hex
################################################################################
#
#
#
#
#
#
#
################################################################################
# arHex: calculates bootstrapped allelic richness and He (divRatio)
################################################################################
# This function will calculate bootstrapped allelic richness and expected
# heterozygosity for use in the Skrbinsek diversity standardization method
arHex <- function(x){
  gp = x$gp
  nloci = ncol(x$ref)
  if(length(x) == 4L){
    pl <- x$ref
    mSize <- x$size
    bser<-function(x){
      return(matrix(x[sample(nrow(x), mSize, replace = TRUE), ], 
                    ncol = ncol(x)))
    }
    pop_list <- bser(pl) # resample
  } else {
    pl <- x$ref
    mSize <- nrow(pl)
    bser<-function(x){
      return(matrix(x[sample(nrow(x), mSize, replace = TRUE), ], 
                    ncol = ncol(x)))
    }
    pop_list <- bser(pl) # resample
  }
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss <- function(x){  # where x is object pop_list
      pl <- list()
      pl[[1]] <- matrix(substr(x, 1, 2), ncol = nloci)
      pl[[2]] <- matrix(substr(x, 3, 4), ncol = nloci)
      return(pl)
    }
  }
  pop_alleles <- pl_ss(pop_list)
  alln <- function(x){
    res <- sapply(1:ncol(x[[1]]), function(i){
      list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
    })
  }
  allele_names <- alln(pop_alleles)
  Alls <- sapply(allele_names, function(x){
    length(x)
  })
  
  #return(Alls)
  #############################################################################
  # Heterozygosity
  aaList <- function(x){
    return(sort(unique(x, decreasing = FALSE)))
  }
  # finter out unique alleles
  all_alleles <- lapply(allele_names, aaList)
  # Create allele frequency holders
  allele_freq <- lapply(1:ncol(pop_list), function(i){
    Nrow <- length(all_alleles[[i]])
    Ncol <- length(pop_list)
    mat <- matrix(rep(0,(Ncol * Nrow)), ncol = Ncol)
    rownames(mat) <- all_alleles[[i]]
    return(mat)
  })
  # rbind pop_alleles
  pa1 <- rbind(pop_alleles[[1]], pop_alleles[[2]])
  
  # Count alleles
  actabPre <- function(x){
    lapply(1:ncol(x), function(i){
      table(x[,i])
    })
  }
  actab <- actabPre(pa1) 
  # Count the number of individuals typed per locus per pop
  indtyppop1 <- function(x){
    apply(x, 2, function(y){
      length(na.omit(y))/2
    })
  }
  indtyppop <- indtyppop1(pa1)
  #calculate allele frequencies
  afCalcpop <- sapply(1:length(actab), function(i){
    actab[[i]]/(indtyppop[i] * 2)
  })
  # calculate heterozygosities
  Hexp <- sapply(afCalcpop, function(x){
    1 - (sum(x^2))
  })
  # return(Hexp)
  return(data.frame(alls = Alls, hexp = Hexp))
}
################################################################################
# End arHex
################################################################################