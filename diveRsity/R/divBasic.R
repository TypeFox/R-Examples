################################################################################
# Calculate basic stats
################################################################################
#' @export
#' HWE exact testing added 29/10/2014
divBasic <- function (infile = NULL, outfile = NULL, gp = 3, bootstraps = NULL,
                      HWEexact = FALSE, mcRep = 2000) {
  
  #infile <- "Test_data.gen"
  #outfile <- "kk"
  #HWEexact <- TRUE
  #mcRep <- 5000
  #fileReader <- diveRsity:::fileReader
  #hweFun <- diveRsity:::hweFun
  #rgp <- diveRsity:::rgp
  #gp <- 3
  #bootstraps = 10
  on = outfile
  # create a results dir
  if(!is.null(on)){
    suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                           "-[diveRsity]","/",sep="")))
    of = paste(getwd(), "/", on, "-[diveRsity]", "/", sep = "")
  }
  
  data1 <- fileReader(infile)
  data1[data1 == 0] <- NA
  data1[data1 == "999999"] <- NA
  data1[data1 == "000000"] <- NA
  data1[data1 == "9999"] <- NA
  data1[data1 == "0000"] <- NA
  #raw_data<-data1
  npops<-length(which(toupper(data1[,1]) == "POP"))
  pop_pos<- c(which(toupper(data1[,1]) == "POP"),(nrow(data1)+1))
  pop_sizes <- sapply(1:npops, function(i){
    pop_pos[(i+1)] - pop_pos[i]-1
  })
  minSize <- min(pop_sizes) 
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
  if (gp == 3){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "    NA"]<-NA
    }
  } else if (gp == 2){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "  NA"] <-NA
    }
  }
  
  
  # define a function for calculating allelic richness 
  ARfun <- function(pop_list){
    
    bser<-function(x){
      return(matrix(x[sample(nrow(x), minSize, replace = TRUE), ],ncol=ncol(x)))
    }
    
    pop_list<-lapply(pop_list, bser)
    
    if (gp==3){
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
        pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
        return(pl)
      }
    } else {
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
        pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
        return(pl)
      }
    }
    pop_alleles<-lapply(pop_list,pl_ss)
    alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
      res<-list()
      for(i in 1:ncol(x[[1]])){
        res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
      }
      return(res)
    }
    allele_names<-lapply(pop_alleles,alln)
    Alls <- lapply(allele_names, function(x){
      sapply(x, function(y){
        length(y)
      })
    })
    nAlls <- matrix(unlist(Alls), ncol = npops)
    
    return(nAlls)
  }
  ############################# END AR function ###############################
  # Calculate allelic richness
  ARdata <- replicate(1000, ARfun(pop_list))
  
  AR <- apply(ARdata, 2, function(x){
    round(rowMeans(x), 2)
  })
  meanAR <- apply(ARdata, 3, colMeans, na.rm = TRUE)
  arLCI <- apply(meanAR, 1, quantile, probs = 0.025, na.rm = TRUE)
  arUCI <- apply(meanAR, 1, quantile, probs = 0.975, na.rm = TRUE)
  #locSD <- apply(ARdata, c(1,2), sd, na.rm = TRUE)
  locSD <- rbind(arLCI, arUCI)
  colnames(locSD) <- pop_names
  rownames(locSD) <- c("Lower_CI", "Upper_CI")
  ###vectorize loci_pop_sizes#################################################
  lps<-function(x){
    lsp_count <- as.vector(colSums(!is.na(x)))
    return(lsp_count)
  }
  locPopSize <- sapply(pop_list,lps)
  ###vectorize pop_alleles####################################################
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
      pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
      return(pl)
    }
  }
  pop_alleles<-lapply(pop_list,pl_ss)
  #end vectorize pop_alleles##################################################
  #vectorize allele_names#####################################################
  pop_alleles<-lapply(pop_list,pl_ss)
  # calcluate the observed heterozygosity
  ohcFUN <- function(x){
    lapply(1:ncol(x[[1]]), function(y){
      (x[[1]][,y]!=x[[2]][,y])*1 #multiply by 1 to conver logical to numeric
    })
  }
  ohc_data <- lapply(pop_alleles, ohcFUN)
  ohcConvert <- function(x){
    matrix(unlist(x),nrow = length(x[[1]]))
  }
  ohc<-lapply(ohc_data,ohcConvert)
  rm(ohc_data)
  hetObs <- sapply(ohc, function(x){
    apply(x, 2, function(y){
      sum(na.omit(y))/length(na.omit(y))
    })
  })
  # End
  alln <- function(x){ # where x is the object pop_alleles (returned by pl_ss())
    res <- sapply(1:ncol(x[[1]]), function(i){
      list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
    })
  }
  allele_names<-sapply(pop_alleles,alln)
  
  # Count the number of alleles observed in each population sample per locus
  obsAlls <- apply(allele_names, 2, function(x){
    sapply(x, function(y){
      length(y)
    })
  })
  # Calculate expected He
  if(npops == 1){
    loci_combi <- allele_names[,1]
  } else {
    loci_combi <- apply(allele_names, 1, FUN = 'unlist')
  }
  # fix loci_combi for SNP format
  if(is.matrix(loci_combi)){
    loci_combi <- lapply(1:ncol(loci_combi), function(i){
      return(loci_combi[,i])
    })
  }
  aaList<-function(x){
    return(sort(unique(x,decreasing=FALSE)))
  }
  all_alleles<-lapply(loci_combi,aaList)
  # Create allele frequency holders
  allele_freq <- lapply(1:ncol(pop_list[[1]]), function(i){
    Nrow <- length(all_alleles[[i]])
    Ncol <- length(pop_list)
    mat <- matrix(rep(0,(Ncol * Nrow)), ncol = Ncol)
    rownames(mat) <- all_alleles[[i]]
    return(mat)
  })
  # rbind pop_alleles
  pa1 <- lapply(pop_alleles, function(x){
    rbind(x[[1]],x[[2]])
  })
  
  # Count alleles
  actab <- lapply(pa1, function(x){
    lapply(1:ncol(x), function(i){
      table(x[,i])
    })
  })
  # Count the number of individuals typed per locus per pop
  indtyppop <- lapply(pa1, function(x){
    apply(x, 2, function(y){
      length(na.omit(y))/2
    })
  })
  #calculate allele frequencies
  afCalcpop<-sapply(1:length(actab), function(x){
    sapply(1:length(actab[[x]]),function(y){
      list(actab[[x]][[y]]/(indtyppop[[x]][y]*2))
    })
  })
  preFreq <- lapply(1:nrow(afCalcpop), function(i){
    lapply(1:ncol(afCalcpop), function(j){
      afCalcpop[i,j][[1]]
    })
  })
  rm(afCalcpop)  # remove afCalcpop
  # Assign allele frequencies per locus
  for(i in 1:nloci){
    for(j in 1:npops){
      allele_freq[[i]][names(preFreq[[i]][[j]]), j] <- preFreq[[i]][[j]]
    }
  }
  # calculate Heterozygosity exp
  if(npops > 1){
    tsapply <- function(...){t(sapply(...))}
    hetExp <- tsapply(allele_freq, function(x){
      apply(x, 2, function(y){
        1 - (sum(y^2))
      })
    })
  } else {
    hetExp <- as.matrix(sapply(allele_freq, function(x){
      apply(x, 2, function(y){
        1 - (sum(y^2))
      })
    }), ncol = 1)
  }
  
  
  totAlls <- sapply(allele_freq, FUN = "nrow")
  # Calculate the proportion of alleles per sample
  propAlls <- apply(obsAlls, 2, function(x){
    round((x/totAlls)*100, 2)
  })
  # R function to calculate expected and observed genetype 
  # numbers for HWE testing
  
  if(HWEexact){
    hwe <- hweFun(infile, mcRep)
    HWE <- round(do.call("cbind", lapply(hwe$locus, function(x){
      sapply(x, "[[", "p")
    })), 3)
    chiDif <- round(do.call("cbind", lapply(hwe$locus, function(x){
      sapply(x, "[[", "chisq")
    })), 3)
    HWEall <- round(sapply(hwe$multilocus, "[[", "p"), 3)
    chi.glb <- sapply(hwe$multilocus, "[[", "chisq")
  } else {
    # generate all possible genotypes for each locus per population
    posGeno <- apply(allele_names, 2, function(x) {
      lapply(x, function(y) {
        if (length(y) == 0) {
          return(NA)
        } else {
          genos <- expand.grid(y, y)
          genos.sort <- t(apply(genos, 1, sort))
          genos <- unique(genos.sort)
          geno <- paste(genos[, 1], genos[, 2], sep = "")
          return(geno)
        }
      })
    })
    # Count the number of each genotype observed
    # define a genotype counting function
    obsGeno <- lapply(1:npops, function(i){
      lapply(1:nloci, function(j){
        sapply(posGeno[[i]][[j]], function(x){
          if(is.na(x)){
            return(NA)
          } else {
            length(which(pop_list[[i]][,j] == x))
          }
        })
      })
    })
    expGeno <- lapply(1:npops, function(i){
      lapply(1:nloci, function(j){
        sapply(posGeno[[i]][[j]], function(x){
          if(is.na(x)){
            return(NA)
          } else {
            if(gp == 3){
              allele1 <- substr(x, 1, 3)
              allele2 <- substr(x, 4, 6)
            } else {
              allele1 <- substr(x, 1, 2)
              allele2 <- substr(x, 3, 4)
            }
            Freq1 <- allele_freq[[j]][which(rownames(allele_freq[[j]]) == 
                                              allele1), i]
            Freq2 <- allele_freq[[j]][which(rownames(allele_freq[[j]]) == 
                                              allele2), i]
            if(allele1 != allele2){
              expFreq <- 2 * (Freq1 * Freq2)
              return(as.vector(expNum <- expFreq * locPopSize[j, i]))
            } else {
              expFreq <- Freq1^2
              return(expNum <- as.vector(expFreq * locPopSize[j, i]))
            }
          }
        })
      })
    })
    # Calculate chi-sq
    chiDif <- sapply(1:npops, function(i){
      sapply(1:nloci, function(j){
        if(length(obsGeno[[i]][[j]]) == 1){
          return(NA)
        } else {
          top <- (obsGeno[[i]][[j]] - expGeno[[i]][[j]])^2
          chi <- top/expGeno[[i]][[j]]
          return(round(sum(chi), 2))
        }      
      })
    })
    # Calculate degrees of freedom
    df <- apply(allele_names, 2,   function(x){
      sapply(x, function(y){
        k <- length(y)
        if(k == 1){
          return(NA)
        } else {
          return((k*(k-1))/2)
        }
      })
    })
    # Calculate HWE significance
    HWE <- sapply(1:npops, function(i){
      round(pchisq(q = chiDif[,i], df = df[,i], lower.tail = FALSE), 4)
    })
    # Calculate over all HWE significance
    HWEall <- round(pchisq(q = colSums(chiDif, na.rm = TRUE), 
                           df = colSums(df, na.rm = TRUE), 
                           lower.tail = FALSE), 4)
  }
  
  if(!is.null(bootstraps)){
    # write a function to calculate allele freq and  obsHet from pop_alleles
    # object
    # convert pop_alleles into a list of arrays
    pa <- lapply(pop_alleles, function(x){
      return(array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), 2)))
    })
    # fit boot function
    bootPA <- function(pa){
      if(!is.list(pa)){
        idx <- sample(dim(pa)[1], dim(pa)[1], replace = TRUE)
        pa <- pa[idx,,]
        obsHet <- apply(pa, 2, function(x){
          (x[,1] != x[, 2])*1
        })
        obsHet <- apply(obsHet, 2, function(x){
          sum(na.omit(x))/length(na.omit(x))
        })
        # calculate expected
        htExp <- apply(pa, 2, function(x){
          af <- as.vector((table(c(x[,1], x[,2]))/(length(na.omit(x[,1]))*2))^2)
          return(1 - sum(af, na.rm = TRUE))
        })
        ht <- sum(htExp, na.rm=TRUE)
        ho <- sum(obsHet, na.rm=TRUE)
        overall <- (ht-ho)/ht
        return(c((htExp-obsHet)/htExp, overall))
      } else {
        out <- lapply(pa, function(pasub){
          idx <- sample(dim(pasub)[1], dim(pasub)[1], replace = TRUE)
          pasub <- pasub[idx,,]
          obsHet <- apply(pasub, 2, function(x){
            (x[,1] != x[, 2])*1
          })
          obsHet <- apply(obsHet, 2, function(x){
            sum(na.omit(x))/length(na.omit(x))
          })
          # calculate expected
          htExp <- apply(pasub, 2, function(x){
            af <- as.vector((table(c(x[,1], x[,2]))/(length(na.omit(x[,1]))*2))^2)
            return(1 - sum(af, na.rm = TRUE))
          })
          ht <- htExp
          ho <- obsHet
          fis <- (ht-ho)/ht
          fis[is.nan(fis)] <- NA
          overall <- (sum(ht,na.rm=TRUE)-sum(ho,na.rm=TRUE))/sum(ht,na.rm=TRUE)
          fis <- c(fis, overall)
          return(fis)
        })
        return(do.call("rbind", out))
      }
    }
    # calculate base fis
    fisCalc <- function(x, y){
      ho <- colSums(x, na.rm = TRUE)
      he <- colSums(y, na.rm = TRUE)
      return((he-ho)/he)
    }
    fisLoc <- round((hetExp[-(nloci+1),] - hetObs[-(nloci+1),])/
                      hetExp[-(nloci+1),], 4)
    fisLoc[is.nan(fisLoc)] <- NA
    fisAct <- fisCalc(hetObs, hetExp)
    # calculate fis CIs
    fisBS <- replicate(bootstraps, bootPA(pa))
    # convert fisBS into list format
    fisBSloc <- lapply(1:npops, function(i){
      return(t(fisBS[i,-(nloci+1),]))
    })
    fisBSOverall <- lapply(1:npops, function(i){
      return(as.vector((fisBS[i,nloci+1,])))
    })
    # fix the bias
    biasCor <- lapply(1:npops, function(i){
      bs <- fisBSloc[[i]]
      mnBA <- colMeans(bs, na.rm = TRUE) - fisLoc[,i]
      mnBA[is.nan(mnBA)] <- NA
      bcCor <- t(apply(bs, 1, function(x){
        return(x - mnBA)
      }))
    })
    biasCorall <- lapply(1:npops, function(i){
      bs <- fisBSOverall[[i]]
      mnBA <- mean(bs, na.rm = TRUE) - fisAct[i]
      mnBA[is.nan(mnBA)] <- NA
      return(as.vector(bs - mnBA))
    })
    
    # bias cor CIs
    bcCILoc <- lapply(biasCor, function(x){
      apply(x, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    })
    bcCIall <- lapply(biasCorall, quantile, 
                      probs = c(0.025, 0.975), na.rm = TRUE)
    nbcCILoc <- lapply(fisBSloc, function(x){
      apply(x, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    })
    nbcCIall <- lapply(fisBSOverall, quantile, 
                       probs = c(0.025, 0.975), na.rm = TRUE)
    fisLoc <- lapply(1:npops, function(i){
      return(fisLoc[,i])
    })
    # arrange data for outpup
    # define function
    tableMake <- function(fisLoc, fisAct, bcCIall, bcCILoc, nbcCIall, nbcCILoc,
                          loci_names){
      out <- data.frame(fis = c(fisLoc, fisAct),
                        lower_CI = c(nbcCILoc[1,], nbcCIall[1]),
                        upper_CI = c(nbcCILoc[2,], nbcCIall[2]),
                        BC_lower_CI = c(bcCILoc[1,], bcCIall[1]),
                        BC_upper_CI = c(bcCILoc[2,], bcCIall[2]))
      rownames(out) <- c(loci_names, "overall")
      return(round(out, 4))
    }
    
    output <- mapply(tableMake, fisLoc = fisLoc,
                     fisAct = fisAct, bcCIall = bcCIall,
                     bcCILoc = bcCILoc, nbcCIall = nbcCIall,
                     nbcCILoc = nbcCILoc,  SIMPLIFY = FALSE,
                     MoreArgs = list(loci_names = loci_names))
    
  }
  # Add pop means or totals to each stat object
  # allelic richness
  AR <- round(rbind(AR, colMeans(AR)), 2)
  dimnames(AR) <- list(c(loci_names, "overall"), pop_names)
  # Number of individuals typed per locus per pop
  locPopSize <- rbind(locPopSize, round(colMeans(locPopSize), 2))
  dimnames(locPopSize) <- list(c(loci_names, "overall"), pop_names)
  # proportion of alleles per pop
  propAlls <- round(rbind(propAlls, colMeans(propAlls)), 2)
  dimnames(propAlls) <- list(c(loci_names, "overall"), pop_names)
  # Number of alleles observed per pop
  obsAlls <- rbind(obsAlls, colSums(obsAlls))
  dimnames(obsAlls) <- list(c(loci_names, "overall"), pop_names)
  # Observed heterozygosity
  hetObs <- round(rbind(hetObs, colMeans(hetObs)), 2)
  dimnames(hetObs) <- list(c(loci_names, "overall"), pop_names)
  # Expected heterozygosity
  hetExp <- round(rbind(hetExp, colMeans(hetExp)), 2)
  dimnames(hetExp) <- list(c(loci_names, "overall"), pop_names)
  # HWE
  HWE <- rbind(HWE, HWEall)
  # Compile information into writable format
  if(!is.null(bootstraps)){
    statComp <- lapply(1:npops, function(i){
      pop <- rbind(locPopSize[,i], obsAlls[,i],
                   propAlls[,i], AR[,i], hetObs[,i],
                   hetExp[,i], HWE[,i], output[[i]][,1],
                   output[[i]][,4], output[[i]][,5])
      return(pop)
    })
    if(npops > 1){
      writeOut <- cbind(c(pop_names[1], 
                          "N", "A", "%", "Ar", "Ho", "He", "HWE", "Fis",
                          "Fis_Low", "Fis_High", "\t"), 
                        rbind(c(loci_names, "Overall"), 
                              statComp[[1]], rep("\t", nloci+1)))
      for(i in 2:npops){
        writeOut <- rbind(writeOut, 
                          cbind(c(pop_names[i], 
                                  "N", "A", "%", "Ar", "Ho", "He", "HWE", "Fis",
                                  "Fis_Low", "Fis_High", "\t"),
                                rbind(c(loci_names, "Overall"), statComp[[i]], 
                                      rep("\t", nloci+1))))
      }
    } else {
      writeOut <- cbind(c(pop_names[1], 
                          "N", "A", "%", "Ar", "Ho", "He", "HWE", "Fis",
                          "Fis_Low", "Fis_High", "\t"), 
                        rbind(c(loci_names, "Overall"), 
                              statComp[[1]], rep("\t", nloci+1)))
    }
  } else {
    statComp <- lapply(1:npops, function(i){
      pop <- rbind(locPopSize[,i], obsAlls[,i],
                   propAlls[,i], AR[,i], hetObs[,i],
                   hetExp[,i], HWE[,i])
      return(pop)
    })
    if(npops > 1){
      writeOut <- cbind(c(pop_names[1], 
                          "N", "A", "%", "Ar", "Ho", "He", "HWE", "\t"), 
                        rbind(c(loci_names, "Overall"), 
                              statComp[[1]], rep("\t", nloci+1)))
      for(i in 2:npops){
        writeOut <- rbind(writeOut, 
                          cbind(c(pop_names[i], 
                                  "N", "A", "%", "Ar", "Ho", "He", "HWE", "\t"),
                                rbind(c(loci_names, "Overall"), statComp[[i]], 
                                      rep("\t", nloci+1))))
      }
    } else {
      writeOut <- cbind(c(pop_names[1], 
                          "N", "A", "%", "Ar", "Ho", "He", "HWE", "\t"), 
                        rbind(c(loci_names, "Overall"), 
                              statComp[[1]], rep("\t", nloci+1)))
    }
  }
  
  
  if (!is.null(outfile)){
    #write_res<-is.element("xlsx",installed.packages()[,1])
    if ("xlsx" %in% rownames(installed.packages())) {
      xlsx::write.xlsx(writeOut, file = paste(of, "[divBasic].xlsx", sep = ""),
                       sheetName = "Basic stats", col.names = FALSE,
                       row.names = FALSE, append=FALSE)
    } else {
      out<-file(paste(of, "[divBasic].txt", sep = ""), "w")
      #cat(paste(colnames(pw_bs_out),sep=""),"\n",sep="\t",file=pw_bts)
      for(i in 1:nrow(writeOut)){
        cat(writeOut[i,], "\n", file = out, sep="\t")
      }
      close(out)
    }
  }
  mainTab <- lapply(statComp, function(x){
    if(nrow(x) == 10L){
      as.data.frame(cbind(stats = c("N", "A", "%", "Ar", "Ho", "He", 
                                    "HWE", "Fis", "Fis_Low", "Fis_High"), 
                          round(x, 3)))
    } else {
      as.data.frame(cbind(stat = c("N", "A", "%", "Ar", "Ho", "He", "HWE"),
                          round(x, 3)))
    }
  })
  names(mainTab) <- pop_names
  if(!is.null(bootstraps)){
    list(locus_pop_size = locPopSize,
         Allele_number = obsAlls,
         proportion_Alleles = propAlls,
         Allelic_richness = AR,
         Ho = hetObs,
         He = hetExp,
         HWE = HWE,
         fis = output,
         arCIs = round(locSD, 4),
         mainTab = mainTab)
  } else {
    list(locus_pop_size = locPopSize,
         Allele_number = obsAlls,
         proportion_Alleles = propAlls,
         Allelic_richness = AR,
         Ho = hetObs,
         He = hetExp,
         HWE = HWE,
         arCIs = round(locSD, 4),
         mainTab = mainTab)
  }
  
}