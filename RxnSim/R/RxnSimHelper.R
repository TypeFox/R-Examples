.algoCheck <- function (algo) {
  if (!(algo %in% c('msim', 'msim_max', 'rsim', 'rsim2'))) {
    stop('Invalid reaction similarity algorithm specified.\n
         Use on of: \'msim\', and \'rsim\'', call. = F)
  }
}

.fpTypeCheck <- function (fp.type, fp.mode) {
  if (!(fp.mode %in% c('bit', 'count'))) {
    stop('Invalid fingerprint mode specificed.', call. = F)
  }
  
  if (fp.type %in% c('standard', 'extended', 'graph', 'hybridization', 'maccs', 
                 'estate', 'kr', 'shortestpath', 'pubchem')) {
    if (fp.mode != 'bit') {
      stop('Use \'bit\' mode for \'', fp.type, '\' fingerprint type.\n
           \'signature\' and \'circular\' are allowed fingerprint types for count mode.', call. = F)
    }
  } else if (fp.type == 'signature') {
    if (fp.mode != 'count') {
      stop('Use \'count\' mode for \'', fp.type, '\' fingerprint type.\n
           \'standard\', \'extended\', \'graph\', \'hybridization\', \'maccs\', \'estate\',
           \'kr\', \'circular\', \'pubchem\' and \'shortestpath\' are allowed fingerprint types for bit mode.', 
           call. = F)
    }
  } else if (fp.type != 'circular') {
    stop('Invalid fingerprint type specificed.', call. = F)
  }
}

.simTypeCheck <- function (sim.method, fp.mode) {
  if (fp.mode == 'bit') {
    if (!(sim.method[[1]] %in% c('simple', 'jaccard', 'tanimoto', 'russelrao', 'dice', 'rodgerstanimoto',
                             'achiai', 'cosine', 'kulczynski2', 'mt', 'baroniurbanibuser', 'tversky',
                             'robust', 'hamann', 'pearson', 'yule', 'mcconnaughey', 'simpson'))) {
      stop('Invalid similarity metric: \'', sim.method[[1]], '\' specified.\n
           Use one of the following metric: \'simple\', \'jaccard\', \'tanimoto\', \'russelrao\',
           \'dice\', \'rodgerstanimoto\', \'achiai\', \'cosine\', \'kulczynski2\', \'mt\',
           \'baroniurbanibuser\', \'tversky\', \'robust\', \'hamann\', \'pearson\', \'yule\',
           \'mcconnaughey\', \'simpson\'', call. = F)
    } else if (sim.method[[1]] == 'tversky') {
      if (length(sim.method) != 3) {
        stop('For Tversky metric, please specify Tversky coefficients. E.g., sim.method = c(\'tversky\', 1, 2).',
             call. = F)
      }
    }
  } else if (fp.mode == 'count') {
    if (!(sim.method[[1]] %in% c('tanimoto', 'dice', 'robust', 'jaccard-count', 'tanimoto-count'))) {
      stop('Invalid similarity metric: \'', sim.method[[1]], '\' specified.\n
            Use \'tanimoto\', \'dice\', \'robust\', \'jaccard-count\' or \'tanimoto-count\' metric 
           for feature-vectors.', 
           call. = F)
    }
  }
}

.smilesParser <- function (smiles, standardize, explicitH) {
  tryCatch({
    mol <- parse.smiles(smiles)[[1]]
    .jcall(.javaObj.env$acm, "V", "percieveAtomTypesAndConfigureAtoms", mol)
    if (standardize) {
      .jcall(.javaObj.env$acm, 'Lorg/openscience/cdk/interfaces/IAtomContainer;', 
             'suppressHydrogens', mol)
    }
    if (explicitH) {
      .jcall(.javaObj.env$acm, 'V', 'convertImplicitToExplicitHydrogens', mol)
    }
    mol
  }, error = function(err) {
    stop("Failed to parse: ", smiles, call. = F)
  })
}

.molParser <- function (fileName, standardize, explicitH) {
  reader <- NULL
  tryCatch({
    file <- .jnew('java.io.FileInputStream', fileName)
    reader <- .jnew('org.openscience.cdk.io.MDLV2000Reader', .jcast(file, 'java.io.InputStream'))
    mol <- .jnew('org.openscience.cdk.AtomContainer')
    objMol <- .jcall(reader, 'Lorg/openscience/cdk/interfaces/IChemObject;', 'read',
                     .jcast(mol,'org.openscience.cdk.interfaces/IChemObject'))
    .jcall(reader, 'V', 'close')
    mol <- .jcast(objMol, 'org.openscience.cdk.interfaces.IAtomContainer')
    .jcall(.javaObj.env$acm, "V", "percieveAtomTypesAndConfigureAtoms", mol)
    if (standardize) {
      .jcall(.javaObj.env$acm, 'Lorg/openscience/cdk/interfaces/IAtomContainer;', 
             'suppressHydrogens', mol)
    }
    if (explicitH) {
      .jcall(.javaObj.env$acm, 'V', 'convertImplicitToExplicitHydrogens', mol)
    }
    mol
  }, error = function(err) {
      stop("Failed to parse file: ", fileName, call. = F)
  }, finally = {
      if(!is.null(reader)) {
        .jcall(reader, 'V', 'close')
      }
  })  
}

.rsmiParser <- function (rsmi, standardize, explicitH) {
  tryCatch({
    rsmi <- gsub (" ", "", rsmi)
    objRxn <- .jcall(.javaObj.env$rs_parser, 'Lorg/openscience/cdk/interfaces/IReaction;',
                     'parseReactionSmiles', rsmi)
    
    rxn <- .jrxnParser (objRxn, standardize, explicitH)
    rsmi <- .jcall(.javaObj.env$smilesGen, 'S', 'createReactionSMILES', objRxn)
    rxn <- c(RSMI = rsmi, rxn)
  }, error = function(err) {
    stop("Failed to parse: ", rsmi, call. = F)
  })  
}

.mdlParser <- function (fileName, standardize, explicitH) {
  reader <- NULL
  tryCatch({
    file <- .jnew('java.io.FileInputStream', fileName)
    reader <- .jnew('org.openscience.cdk.io.MDLRXNV2000Reader', .jcast(file, 'java.io.InputStream'))
    rct <- .jnew('org.openscience.cdk.Reaction')
    objRxn <- .jcall(reader, 'Lorg/openscience/cdk/interfaces/IChemObject;', 'read',
                     .jcast(rct,'org.openscience.cdk.interfaces/IChemObject'))
    objRxn <- .jcast(objRxn, 'org.openscience.cdk.interfaces.IReaction')
    .jcall(reader, 'V', 'close')
    
    rxn <- .jrxnParser (objRxn, standardize, explicitH)
    rsmi <- .jcall(.javaObj.env$smilesGen, 'S', 'createReactionSMILES', objRxn)
    rxn <- c(RSMI = rsmi, rxn)
  }, error = function(err) {
    stop("Failed to parse file: ", fileName, call. = F)
  }, finally = {
    if(!is.null(reader)) {
      .jcall(reader, 'V', 'close')
    }
  })
}

.jrxnParser <- function (objRxn, standardize, explicitH) {
  tryCatch({
    objReacts <- .jcall(objRxn, 'Lorg/openscience/cdk/interfaces/IAtomContainerSet;',
                        'getReactants')
    Reacts <- as.list(.jcall(objReacts, 'Ljava/lang/Iterable;', 'atomContainers'))
    Reacts <- lapply(Reacts, .jcast, 'org/openscience/cdk/interfaces/IAtomContainer')
    
    objProds <- .jcall(objRxn, 'Lorg/openscience/cdk/interfaces/IAtomContainerSet;',
                       'getProducts')
    Prods <- as.list(.jcall(objProds, 'Ljava/lang/Iterable;', 'atomContainers'))
    Prods <- lapply(Prods, .jcast, 'org/openscience/cdk/interfaces/IAtomContainer')
    
    for (mol in Reacts) {
      .jcall(.javaObj.env$acm, "V", "percieveAtomTypesAndConfigureAtoms", mol)
    }
    for (mol in Prods) {
      .jcall(.javaObj.env$acm, "V", "percieveAtomTypesAndConfigureAtoms", mol)
    }
    
    if (standardize) {
      for (mol in Reacts) {
        .jcall(.javaObj.env$acm, 'Lorg/openscience/cdk/interfaces/IAtomContainer;', 
               'suppressHydrogens', mol)
      }
      for (mol in Prods) {
        .jcall(.javaObj.env$acm, 'Lorg/openscience/cdk/interfaces/IAtomContainer;',
               'suppressHydrogens', mol)
      }
    }
    
    if (explicitH) {
      for (mol in Reacts) {
        .jcall(.javaObj.env$acm, 'V', 'convertImplicitToExplicitHydrogens', mol)
      }
      for (mol in Prods) {
        .jcall(.javaObj.env$acm, 'V', 'convertImplicitToExplicitHydrogens', mol)
      }
    }
    
    rxn <- list(Reactants = Reacts, Products = Prods)
  }, error = function(err) {
    stop('.jrxnParser: ', err, '\n', call. = F)
  })
}

.similarity <- function (rxnA, rxnB, reversible, algo, sim.method, fp.type, fp.mode, 
                         fp.depth, fp.size, verbose = F, cached = F) {
  if (cached) {
    cache <- .fp.env$fp_map
  } else {
    cache <- NULL
  }
  
  fpA_r <- lapply (rxnA$Reactants, .makeFP, fp.type = fp.type, fp.mode = fp.mode,
                   fp.depth = fp.depth, fp.size = fp.size, cache)
  fpA_p <- lapply (rxnA$Products, .makeFP, fp.type = fp.type, fp.mode = fp.mode,
                   fp.depth = fp.depth, fp.size = fp.size, cache)
  fpB_r <- lapply (rxnB$Reactants, .makeFP, fp.type = fp.type, fp.mode = fp.mode,
                   fp.depth = fp.depth, fp.size = fp.size, cache)
  fpB_p <- lapply (rxnB$Products, .makeFP, fp.type = fp.type, fp.mode = fp.mode,
                   fp.depth = fp.depth, fp.size = fp.size, cache)
  
  if (sum(sapply(c(fpA_r, fpA_p, fpB_r, fpB_p), function(x) {if(is.null(x)){1} else {0}}))) {
    warning('Could not generate fingerprints. Try reducing fingerprint search depth.', call. = F)
    NaN
  }
  else {
    smi <- new.env(parent = emptyenv())
    if (verbose) {
      smi$A_r <- lapply(rxnA$Reactants, get.smiles)
      smi$A_p <- lapply(rxnA$Products, get.smiles)
      smi$B_r <- lapply(rxnA$Reactants, get.smiles)
      smi$B_p <- lapply(rxnA$Products, get.smiles)
    }
    
    .calcSimilarity (fpA_r, fpA_p, fpB_r, fpB_p, reversible, algo, sim.method, verbose, smi)
  }
}

.makeFP <- function (mol, fp.type, fp.mode, fp.depth, fp.size, cache) {
  if (!missing(cache) && !is.null(cache)) {
    smi <- get.smiles(mol, type = 'unique')
    fp <- cache[[smi]]
    if (is.null(fp)) {
      fp <- get.fingerprint(mol, type = fp.type, fp.mode = fp.mode, depth = fp.depth, 
                            size = fp.size, verbose = T)
      cache[[smi]] <- fp
    }
  } else {
    fp <- get.fingerprint(mol, type = fp.type, fp.mode = fp.mode, depth = fp.depth,
                          size = fp.size, verbose = T)
  }
  fp
}

.calcSimilarity <- function (fpA_r, fpA_p, fpB_r, fpB_p, reversible, algo, sim.method, verbose = F, smiles = NULL) {
  maxRS <- FALSE
  if (algo == 'msim_max') {
    algo <- 'msim'
    maxRS <- TRUE
  }
  if (algo == 'msim') {    
    lenA_r <- length(fpA_r)
    lenB_r <- length(fpB_r)
    lenA_p <- length(fpA_p)
    lenB_p <- length(fpB_p)
    
    dfrr <- .calcSimMapping(fpA_r, fpB_r, sim.method)
    dfpp <- .calcSimMapping(fpA_p, fpB_p, sim.method)
    lenDFrr <- nrow(dfrr)
    lenDFpp <- nrow(dfpp)
    divFac <- ifelse(maxRS == T, (lenDFrr + lenDFpp),
                     ((lenA_r + lenB_r - lenDFrr)+(lenA_p + lenB_p - lenDFpp)))    
    straight <- (sum(dfrr[3]) + sum(dfpp[3]))/divFac
    
    if(reversible) {
      dfrp <- .calcSimMapping(fpA_r, fpB_p, sim.method)
      dfpr <- .calcSimMapping(fpA_p, fpB_r, sim.method)
      lenDFrp <- nrow(dfrp)
      lenDFpr <- nrow(dfpr)
      divFac <- ifelse(maxRS == T, (lenDFrp + lenDFpr),
                       ((lenA_r + lenB_p - lenDFrp)+(lenA_p + lenB_r - lenDFpr)))
      cross <- (sum(dfrp[3]) + sum(dfpr[3]))/divFac
    } else {
      cross <- 0
    }
    
    if (verbose) {
      cat('RCT1\n')
      cat('- React(s)\n')
      i <- 1
      for (mol in smiles$A_r) {
        cat('\t\t', paste(i, mol, sep='\t'), '\n')
        i <- i+1
      }
      cat('- Prod(s)\n')
      i <- 1
      for (mol in smiles$A_p) {
        cat('\t\t', paste(i, mol, sep='\t'), '\n')
        i <- i+1
      }
      cat('RCT2\n')
      cat('- React(s)\n')
      i <- 1
      for (mol in smiles$B_r) {
        cat('\t\t', paste(i, mol, sep='\t'), '\n')
        i <- i+1
      }
      cat('- Prod(s)\n')
      i <- 1
      for (mol in smiles$B_p) {
        cat('\t\t', paste(i, mol, sep='\t'), '\n')
        i <- i+1
      }
      cat('\n')
      if(straight >= cross) {
        if (!maxRS) {
          if (lenA_r > lenDFrr) {
            ids <- c(1:lenA_r)
            ids <- ids[!ids %in% dfrr[,1]]
            i <- 1
            for (id in ids) {
              dfrr[lenDFrr + i, ] <- c(id, '-', 0)
              i <- i + 1
            }
          } else if (lenB_r > lenDFrr) {
            ids <- c(1:lenB_r)
            ids <- ids[!ids %in% dfrr[,2]]
            i <- 1
            for (id in ids) {
              dfrr[lenDFrr + i, ] <- c('-', id, 0)
              i <- i + 1
            }
          }
          if (lenA_p > lenDFpp) {
            ids <- c(1:lenA_p)
            ids <- ids[!ids %in% dfpp[,1]]
            i <- 1
            for (id in ids) {
              dfpp[lenDFpp + i, ] <- c(id, '-', 0)
              i <- i + 1
            }
          } else if (lenB_p > lenDFpp) {
            ids <- c(1:lenB_p)
            ids <- ids[!ids %in% dfpp[,2]]
            i <- 1
            for (id in ids) {
              dfpp[lenDFpp + i, ] <- c('-', id, 0)
              i <- i + 1
            }
          }
        }
        colnames(dfrr) <- c('RCT1-React', 'RCT2-React', 'Similarity')
        row.names(dfrr) <- c(1:length(dfrr[[1]]))
        colnames(dfpp) <- c('RCT1-Prod', 'RCT2-Prod', 'Similarity')
        row.names(dfpp) <- c(1:length(dfpp[[1]]))
        print(dfrr, row.names = F)
        print(dfpp, row.names = F)
      } else {
        if (!maxRS) {
          if (lenA_r > lenDFrp) {
            ids <- c(1:lenA_r)
            ids <- ids[!ids %in% dfrp[,1]]
            i <- 1
            for (id in ids) {
              dfrp[lenDFrp + i, ] <- c(id, '-', 0)
              i <- i + 1
            }
          } else if (lenB_p > lenDFrp) {
            ids <- c(1:lenB_p)
            ids <- ids[!ids %in% dfrp[,2]]
            i <- 1
            for (id in ids) {
              dfrp[lenDFrp + i, ] <- c('-', id, 0)
              i <- i + 1
            }
          }
          if (lenA_p > lenDFpr) {
            ids <- c(1:lenA_p)
            ids <- ids[!ids %in% dfpr[,1]]
            i <- 1
            for (id in ids) {
              dfpr[lenDFpr + i, ] <- c(id, '-', 0)
              i <- i + 1
            }
          } else if (lenB_r > lenDFpr) {
            ids <- c(1:lenB_r)
            ids <- ids[!ids %in% dfpr[,2]]
            i <- 1
            for (id in ids) {
              dfpr[lenDFpr + i, ] <- c('-', id, 0)
              i <- i + 1
            }
          }
        }
        colnames(dfrp) <- c('RCT1-React', 'RCT2-Prod', 'Similarity')
        row.names(dfrp) <- c(1:length(dfrp[[1]]))
        colnames(dfpr) <- c('RCT1-Prod', 'RCT2-React', 'Similarity')
        row.names(dfpr) <- c(1:length(dfpr[[1]]))
        
        print(dfrp, row.names = F)
        print(dfpr, row.names = F)
      }
      cat('\n')
      cat('Reaction similarity:\t', ifelse(straight >= cross, straight, cross), '\n\n')
    }
    
    ifelse(straight >= cross, straight, cross)
  } else if (algo == 'rsim') {
    fpAR <- .addFP(fpA_r)
    fpAP <- .addFP(fpA_p)
    fpBR <- .addFP(fpB_r)
    fpBP <- .addFP(fpB_p)
    
    simRR <- .calcDistance(fpAR, fpBR, sim.method = sim.method)
    simPP <- .calcDistance(fpAP, fpBP, sim.method = sim.method)
    simRRPP <- (simRR + simPP)/2
    
    if(reversible) {
      simRP <- .calcDistance(fpAR, fpBP, sim.method = sim.method)
      simPR <- .calcDistance(fpAP, fpBR, sim.method = sim.method)
      simRPPR <- (simRP + simPR)/2
    } else {
      simRPPR <- 0
    }
    
    if (verbose) {
      if (simRRPP > simRPPR) {
        cat('Rct1-Reactant(s) | Rct2-Reactant(s):', simRR, '\n')
        cat('Rct1-Product(s)  |  Rct2-Product(s):', simPP, '\n')
      } else {
        cat('Rct1-Reactant(s) |  Rct2-Product(s):', simRP, '\n')
        cat('Rct1-Product(s)  | Rct2-Reactant(s):', simPR, '\n')
      }
      cat('\n')
      cat('Reaction similarity:\t', ifelse (simRRPP > simRPPR, simRRPP, simRPPR), '\n\n')
    }
    ifelse (simRRPP > simRPPR, simRRPP, simRPPR)
  } else if (algo == 'rsim2') {
    fpA <- .addFP(c(fpA_r, fpA_p))
    fpB <- .addFP(c(fpB_r, fpB_p))
    
    simR <- .calcDistance(fpA, fpB, sim.method = sim.method)
    
    if (verbose) {
      cat('\n')
      cat('Reaction similarity:\t', simR, '\n\n')
    }
    
    simR
  }
}

.calcSimMapping <- function (fpA, fpB, sim.method) {
  
  dfSIM <- expand.grid(1:length(fpA), 1:length(fpB))
  SIM <- mapply(function(i, j) {.calcDistance(fpA = fpA[[i]], fpB = fpB[[j]], sim.method)}, dfSIM[,1], dfSIM[,2])
  dfSIM$SIM <- SIM
  colnames(dfSIM) <- c('ID1', 'ID2', 'SIMILARITY')
  
  dfSIM <- dfSIM[order(-dfSIM$SIMILARITY),]
  
  dfSIM$FLAG <- TRUE
  itr <- 1
  while (itr <= nrow(dfSIM)) {
    id1 = dfSIM[itr,1]
    id2 = dfSIM[itr,2]
    
    dfSIM[dfSIM$ID1 == id1 | dfSIM$ID2 == id2, 4] <- FALSE
    dfSIM[itr, 4] <- TRUE
    
    dfSIM <- dfSIM[dfSIM$FLAG == TRUE,]
    itr <- itr + 1
  }
  dfSIM <- dfSIM[,1:3]
}

.addFP <- function (fpList) {
  if (class(fpList[[1]]) == 'featvec') {
    featrs_map <- new.env(parent = emptyenv(), hash = T)
    for (fp in fpList) {
      for (featr in fp@features) {
        f <- fingerprint::feature(featr)
        c <- fingerprint::count(featr)
        if(is.null(featrs_map[[f]])) {
          featrs_map[[f]] <- c
        } else {
          featrs_map[[f]] <- featrs_map[[f]] + c
        }
      }
    }
    f <- new('feature')
    ftList <- list()
    for (featr in ls(featrs_map)) {
      fingerprint::feature(f) <- featr
      fingerprint::count(f) <- featrs_map[[featr]]
      ftList[[length(ftList) + 1]] <- f
    }
    fpSum <- new('featvec')
    fpSum@features <- ftList
    fpSum
  } else if (class(fpList[[1]]) == 'fingerprint') {
    fpSum <- fpList[[1]]
    for (fp in fpList[-1]) {
      fpSum <- fpSum | fp
    }
    fpSum
  } else {
    stop('Undefined fingerprint class.', class. = F)
  }
}

.JaccardCount <- function (fpA, fpB) {
  if (class(fpA) != 'featvec' || class(fpB) != 'featvec') {
    stop('Inputs should be of \'featvec\' (S4 class) type.')
  }
  
  featrs_mapA <- new.env(parent = emptyenv(), hash = T)
  for (featr in fpA@features) {
    f <- fingerprint::feature(featr)
    c <- fingerprint::count(featr)
    featrs_mapA[[f]] <- c
  }
  featrs_mapB <- new.env(parent = emptyenv(), hash = T)
  for (featr in fpB@features) {
    f <- fingerprint::feature(featr)
    c <- fingerprint::count(featr)
    featrs_mapB[[f]] <- c
  }
  ftA <- ls(featrs_mapA)
  ftB <- ls(featrs_mapB)
  min <- 0
  max <- 0
  for (f in union(ftA, ftB)) {
    cA <- ifelse(is.null(featrs_mapA[[f]]), 0, featrs_mapA[[f]])
    cB <- ifelse(is.null(featrs_mapB[[f]]), 0, featrs_mapB[[f]])
    min <- min + ifelse (cA > cB, cB, cA)
    max <- max + ifelse (cA > cB, cA, cB)
  }
  min/max
}

.TanimotoCount <- function (fpA, fpB) {
  if (class(fpA) != 'featvec' || class(fpB) != 'featvec') {
    stop('Inputs should be of \'featvec\' (S4 class) type.')
  }
  
  featrs_mapA <- new.env(parent = emptyenv(), hash = T)
  sumA <- 0
  for (featr in fpA@features) {
    f <- fingerprint::feature(featr)
    c <- fingerprint::count(featr)
    featrs_mapA[[f]] <- c
    sumA <- sumA + (c * c)
  }
  sumB <- 0
  featrs_mapB <- new.env(parent = emptyenv(), hash = T)
  for (featr in fpB@features) {
    f <- fingerprint::feature(featr)
    c <- fingerprint::count(featr)
    featrs_mapB[[f]] <- c
    sumB <- sumB + (c * c)
  }
  ftA <- ls(featrs_mapA)
  ftB <- ls(featrs_mapB)
  min <- 0
  max <- 0
  sum_intersect <- 0
  for (f in intersect(ftA, ftB)) {
    cA <- featrs_mapA[[f]]
    cB <- featrs_mapB[[f]]
    sum_intersect <- sum_intersect + (cA * cB)
  }
  sum_intersect/(sumA + sumB - sum_intersect)
}

.calcDistance <- function (fpA, fpB, sim.method) {
  if (sim.method[[1]] == 'jaccard-count') {
    .JaccardCount(fpA, fpB)
  } else if (sim.method[[1]] == 'tanimoto-count') {
    .TanimotoCount(fpA, fpB)
  } else if (length(sim.method) != 3) {
    fingerprint::distance (fpA, fpB, method = sim.method[[1]])
  } else {
    fingerprint::distance (fpA, fpB, method = sim.method[[1]], a = as.numeric(sim.method[[2]]), b = as.numeric(sim.method[[3]]))
  }
}