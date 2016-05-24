.global.env <- new.env(parent = emptyenv())
.javaObj.env <- new.env(parent = emptyenv())
.fp.env <- new.env(parent = emptyenv())

.onAttach <- function(libname, pkgname) {
  .javaObj.env$rs_parser <- rcdk::get.smiles.parser()
  .javaObj.env$smilesGen <- .jnew ('org.openscience.cdk.smiles.SmilesGenerator')
  .javaObj.env$acm <- 'org.openscience.cdk.tools.manipulator.AtomContainerManipulator'
  .fp.env$fp_map <- new.env(parent = emptyenv(), hash = T)
  .global.env$DefaultDB <- paste(libname, pkgname, 'DB/Metadata.txt', sep = '/')
  .global.env$Rheav60 <- paste(libname, pkgname, 'DB/RheaData_v60.txt', sep = '/')
}

rs.clearCache <- function () {
  .fp.env$fp_map <- new.env(parent = emptyenv(), hash = T)
}

ms.compute <- function (molA, molB, format = 'smiles', standardize = T, explicitH = F, sim.method = 'tanimoto',
                        fp.type = 'extended', fp.mode = 'bit', fp.depth = 6, fp.size = 1024,
                        fpCached = F) {
  
  format <- tolower(format)
  sim.method[[1]] <- tolower(sim.method[[1]])
  fp.type <- tolower(fp.type)
  fp.mode <- tolower(fp.mode)
  
  if (missing(molA) || missing(molB)) {
    stop("Two inputs needed to compute similarity", call. = F)
  } else if(length(molA) > 1 || length(molB) > 1) {
    warning("Input(s) has length > 1 and only the first element(s) will be used.")
    molA <- molA[[1]]
    molB <- molB[[1]]
  }
  
  if(!is.character(molA) || !is.character(molB)) {
    stop("Invalid input. Enter in SMILES format or path to MOL file.", call. = F)
  }
  
  .fpTypeCheck(fp.type, fp.mode)
  .simTypeCheck(sim.method, fp.mode)
  if (length(format) == 1) {
    format[[2]] <- format[[1]]
  }
  
  tryCatch ({
    if (format[[1]] == 'smiles') {
      mA <- .smilesParser(molA, standardize, explicitH)
    } else if (format[[1]] == 'mol') {
      mA <- .molParser(molA, standardize, explicitH)
    } else {
      stop("Invalid input format.", call. = F)
    }
  
    if (format[[2]] == 'smiles') {
      mB <- .smilesParser(molB, standardize, explicitH)
    } else if (format[[2]] == 'mol') {
      mB <- .molParser(molB, standardize, explicitH)
    } else {
      stop("Invalid input format.", call. = F)
    }
    
    if (fpCached) {
      cache <- .fp.env$fp_map
    } else {
      cache <- NULL
    }
    fA <- .makeFP(mA, fp.type, fp.mode, fp.depth, fp.size, cache)
    fB <- .makeFP(mB, fp.type, fp.mode, fp.depth, fp.size, cache)
    result <- .calcDistance (fA, fB, sim.method)
  }, error = function(err) {
    stop(err)
  })
  result
}

ms.compute.sim.matrix <- function (molA, format = 'smiles', standardize = T, explicitH = F, sim.method = 'tanimoto',
                                fp.type = 'extended', fp.mode = 'bit', fp.depth = 6, fp.size = 1024, clearCache = T) {
  
  format <- tolower(format)
  sim.method[[1]] <- tolower(sim.method[[1]])
  fp.type <- tolower(fp.type)
  fp.mode <- tolower(fp.mode)
  
  if (missing(molA) || length(molA) < 2) {
    stop("Pass two or more molecules to compute similarity.", call. = F)
  }
  
  .fpTypeCheck(fp.type, fp.mode)
  .simTypeCheck(sim.method, fp.mode)
  
  result <- tryCatch({
    if (format[[1]] == 'smiles') {
      molSObj <- lapply(molA, .smilesParser, standardize, explicitH)
    } else if (format[[1]] == 'mol') {
      molSObj <- lapply(molA, .molParser, standardize, explicitH)
    } else {
      stop("Invalid input format.", call. = F)
    }
    
    len <- length(molA)
    sim <- matrix(0, nrow = len, ncol = len)
    
    if (clearCache) {
      rs.clearCache()
    }
    fpMolS <- lapply(molSObj, .makeFP, fp.type, fp.mode, fp.depth, fp.size, .fp.env$fp_map)
    for (i in 1:(len - 1)) {
      v <- unlist(lapply(fpMolS[(i+1):len], .calcDistance, fpMolS[[i]], sim.method))
      sim[i, (i+1):len] <- v
      sim[(i+1):len, i] <- v
    }
    if (clearCache) {
      rs.clearCache()
    }
    diag(sim) <- 1
    sim
  }, error = function (err) {
    stop (err)
  })
  return(result)
}

rs.compute <- function (rxnA, rxnB, format = 'rsmi', standardize = T, explicitH = F, reversible = T,
                        algo = 'msim', sim.method = 'tanimoto', fp.type = 'extended', fp.mode = 'bit',
                        fp.depth = 6, fp.size = 1024, verbose = F, fpCached = F) {
  
  format <- tolower(format)
  algo <- tolower(algo)
  sim.method[[1]] <- tolower(sim.method[[1]])
  fp.type <- tolower(fp.type)
  fp.mode <- tolower(fp.mode)
  
  if (missing(rxnA) || missing(rxnB)) {
    stop("Two reactions needed to compute similarity", call. = F)
  } else if(length(rxnA) > 1 || length(rxnB) > 1) {
    warning("Input(s) has length > 1 and only the first element(s) will be used.")
    rxnA <- rxnA[[1]]
    rxnB <- rxnB[[1]]
  }
  if(!is.character(rxnA) || !is.character(rxnB)) {
    stop("Invalid input. Enter in (REACTION) SMILES format or path to RXN file.", call. = F)
  }
  
  .algoCheck(algo)
  .fpTypeCheck(fp.type, fp.mode)
  .simTypeCheck(sim.method, fp.mode)
  if (length(format) == 1) {
    format[[2]] <- format[[1]]
  }
  
  result <- tryCatch({
    if (format[[1]] == 'rsmi') {
      rA <- .rsmiParser(rxnA, standardize, explicitH)
    } else if (format[[1]] == 'rxn') {
      rA <- .mdlParser(rxnA, standardize, explicitH)
    } else {
      stop("Invalid input format.", call. = F)
    }
    if (format[[2]] == 'rsmi') {
      rB <- .rsmiParser(rxnB, standardize, explicitH)
    } else if (format[[2]] == 'rxn') {
      rB <- .mdlParser(rxnB, standardize, explicitH)
    } else {
      stop("Invalid input format.", call. = F)
    }
    
    .similarity(rA, rB, reversible = reversible, algo = algo, sim.method = sim.method, fp.type = fp.type,
                fp.mode = fp.mode, fp.depth = fp.depth, fp.size = fp.size, verbose = verbose, 
                cached = fpCached)
  }, error = function(err) {
    stop(err)
  })
  
  result
}

rs.compute.list <- function (rxnA, rxnB, format = 'rsmi', standardize = T, explicitH = F, reversible = T,
                             algo = 'msim', sim.method = 'tanimoto', fp.type = 'extended',
                             fp.mode = 'bit', fp.depth = 6, fp.size = 1024, clearCache = T) {
  
  format <- tolower(format)
  algo <- tolower(algo)
  sim.method[[1]] <- tolower(sim.method[[1]])
  fp.type <- tolower(fp.type)
  fp.mode <- tolower(fp.mode)
  
  if (missing(rxnA) || missing(rxnB)) {
    stop("Pass two reaction lists to compute similarity.", call. = F)
  }
  
  .algoCheck(algo)
  .fpTypeCheck(fp.type, fp.mode)
  .simTypeCheck(sim.method, fp.mode)
  if (length(format) == 1) {
    format[[2]] <- format[[1]]
  }
  
  result <- tryCatch({
    if (format[[1]] == 'rsmi') {
      rxnAObjs <- lapply(rxnA, .rsmiParser, standardize, explicitH)
    } else if (format[[1]] == 'rxn') {
      rxnAObjs <- lapply(rxnA, .mdlParser, standardize, explicitH)
    } else {
      stop("Invalid input format.", call. = F)
    }
    if (format[[2]] == 'rsmi') {
      rxnBObjs <- lapply(rxnB, .rsmiParser, standardize, explicitH)
    } else if (format[[2]] == 'rxn') {
      rxnBObjs <- lapply(rxnB, .mdlParser, standardize, explicitH)
    } else {
      stop("Invalid input format.", call. = F)
    }
    
    len <- length(rxnA)
    sim <- matrix(0, nrow = len, ncol = length(rxnB))
    
    if (clearCache) {
      rs.clearCache()
    }
    for (i in 1:len) {
      v <- unlist(lapply(rxnBObjs, .similarity, rxnA = rxnAObjs[[i]], reversible = reversible, algo = algo,
                         sim.method = sim.method, fp.type = fp.type, fp.mode = fp.mode, fp.depth = fp.depth,
                         fp.size = fp.size, cached = T))
      sim[i,] <- v
    }
    if (clearCache) {
      rs.clearCache()
    }
    sim
  }, error = function (err) {
    stop (err)
  })
  return(result)
}

rs.compute.sim.matrix <- function (rxnA, format = 'rsmi', standardize = T, explicitH = F, reversible = T, 
                                   algo = 'msim', sim.method = 'tanimoto', fp.type = 'extended',
                                   fp.mode = 'bit', fp.depth = 6, fp.size = 1024, clearCache = T) {
  
  format <- tolower(format)
  algo <- tolower(algo)
  sim.method[[1]] <- tolower(sim.method[[1]])
  fp.type <- tolower(fp.type)
  fp.mode <- tolower(fp.mode)
  
  if (missing(rxnA) || length(rxnA) < 2) {
    stop("Pass two or more reactions to compute similarity.", call. = F)
  }
  
  .algoCheck(algo)
  .fpTypeCheck(fp.type, fp.mode)
  .simTypeCheck(sim.method, fp.mode)
  
  result <- tryCatch({
    if (format[[1]] == 'rsmi') {
      rxnSObj <- lapply(rxnA, .rsmiParser, standardize, explicitH)
    } else if (format[[1]] == 'rxn') {
      rxnSObj <- lapply(rxnA, .mdlParser, standardize, explicitH)
    } else {
      stop("Invalid input format.", call. = F)
    }
    
    len <- length(rxnA)
    sim <- matrix(0, nrow = len, ncol = len)
    
    if (clearCache) {
      rs.clearCache()
    }
    for (i in 1:(len - 1)) {
      v <- unlist(lapply(rxnSObj[(i+1):len], .similarity, rxnB = rxnSObj[[i]], reversible = reversible, 
                         algo = algo, sim.method = sim.method, fp.type = fp.type, fp.mode = fp.mode,
                         fp.depth = fp.depth, fp.size = fp.size, cached = T))
      sim[i, (i+1):len] <- v
      sim[(i+1):len, i] <- v
    }
    if (clearCache) {
      rs.clearCache()
    }
    diag(sim) <- 1
    sim
  }, error = function (err) {
    stop (err)
  })
  return(result)
}

rs.makeDB <- function (txtFile, header = F, sep = '\t', standardize = T, explicitH = F, fp.type = 
                      'extended', fp.mode = 'bit', fp.depth = 6, fp.size = 1024, useMask = F,
                       maskStructure, mask, recursive = F) {
  
  fp.type <- tolower(fp.type)
  fp.mode <- tolower(fp.mode)
  
  if (missing(txtFile)) {
    msg <- paste ("DB text file not provided, using sample DB (extracted from Rhea v.60).\n",
                  "For complete dataset use: ", .global.env$Rheav60, sep = '')
    warning(msg, call. = F, immediate. = T)
    txtFile <- .global.env$DefaultDB
  }
  
  .fpTypeCheck(fp.type, fp.mode)
  
  if (useMask == T) {
    if (missing(maskStructure) || maskStructure == '') {
      stop('Enter a structure to mask in form of a SMILES or SMARTS.', call. = F)
    }
  }
  
  DB <- NULL
  tryCatch({
    #DB <- read.delim(txtFile, header=header, sep=sep, strip.white=T)
    DB <- data.table::fread(txtFile, header=header, sep=sep, data.table=F)
    colnames(DB) <- c('EC', 'ID', 'RSMI')
    
    rxnObjList <- lapply(as.character(DB$RSMI), .rsmiParser, standardize, explicitH)
    if (useMask == T) {
      rxnObjList <- lapply(rxnObjList, .rct.mask, substructure = maskStructure, mask = mask, recursive = recursive)
      MaskedRSMI <- unlist(lapply(rxnObjList, function(obj) {obj[[1]]}))
      DB <- cbind(DB, MaskedRSMI)
    }
    
    fp_map_cache <- new.env(parent = emptyenv(), hash = T)
    fpList <- list()
    len <- 0
    for (obj in rxnObjList) {
      fp_r <- lapply (obj$Reactants, .makeFP, fp.type = fp.type, fp.mode = fp.mode,
                       fp.depth = fp.depth, fp.size = fp.size, fp_map_cache)
      fp_p <- lapply (obj$Products, .makeFP, fp.type = fp.type, fp.mode = fp.mode,
                       fp.depth = fp.depth, fp.size = fp.size, fp_map_cache)
      
      if (sum(sapply(c(fp_r, fp_p), function(x) {if(is.null(x)){1} else {0}}))) {
        message(paste('Skipping reaction -', DB$ID[[len+1]], '- could not generate fingerprints.'))
        DB <- DB[-(len+1),]
      } else {
        len <- len + 1
        fpList[[len]] <- list(FPR = fp_r, FPP = fp_p)
      }
    }
    list(Data = DB, FP = fpList, standardize = standardize, explicitH = explicitH, fp.type = fp.type,
         fp.mode = fp.mode, fp.depth = fp.depth, fp.size = fp.size)
  }, error = function(err) {
    stop(err)
  })
}

rs.compute.DB <- function (rxnA, DB, format = 'rsmi', ecrange = '*', reversible = T, algo = 'msim',
                           sim.method = 'tanimoto', sort = T, fpCached = F) {
  
  format <- tolower(format)
  sim.method[[1]] <- tolower(sim.method[[1]])
  algo <- tolower(algo)

  if (missing(rxnA)) {
    stop("Reaction needed to compute similarity", call. = F)
  } else if(length(rxnA) > 1) {
    warning("Input(s) has length > 1 and only the first element(s) will be used.")
    rxnA <- rxnA[[1]]
  }
  if(!is.character(rxnA)) {
    stop("Invalid input. Enter in SMILES format or path to RXN file.", call. = F)
  }
  
  if (missing(DB)) {
    stop("DB not specified. Use \'rs.makeDB\' to create DB.")
  }
  
  .algoCheck(algo)
  .simTypeCheck(sim.method, DB$fp.mode)
  
  result <- tryCatch({
    if (format[[1]] == 'rsmi') {
      rct <- .rsmiParser(rxnA, DB$standardize, DB$explicitH)
    } else if (format[[1]] == 'rxn') {
      rct <- .mdlParser(rxnA, DB$standardize, DB$explicitH)
    } else {
      stop("Invalid input format.", call. = F)
    }
    
    if (fpCached) {
      cache <- .fp.env$fp_map
    } else {
      cache <- NULL
    }
    fp_r <- lapply (rct$Reactants, .makeFP, fp.type = DB$fp.type, fp.mode = DB$fp.mode,
                     fp.depth = DB$fp.depth, fp.size = DB$fp.size, cache)
    fp_p <- lapply (rct$Products, .makeFP, fp.type = DB$fp.type, fp.mode = DB$fp.mode,
                     fp.depth = DB$fp.depth, fp.size = DB$fp.size, cache)
        
    ecrange <- gsub('\\.', '\\\\\\.', ecrange)
    ecrange <- gsub('\\*', '.*', ecrange)
    
    #resDF <- as.data.frame(setNames(replicate(length(DB$Data)+1,numeric(0), simplify = F), c(colnames(DB$Data),'SIMILARITY')))
    #for (itr in grep(paste('^', ecrange, sep = ''), DB$Data$EC)) {
    #  sim <- .calcSimilarity (fp_r, fp_p, DB$FP[[itr]]$FPR, DB$FP[[itr]]$FPP, reversible, algo, sim.method)
    #  res <- data.frame(DB$Data[itr,], SIMILARITY = sim)
    #  resDF <- rbind(resDF, res)
    #}
    
    ECs <- grep(paste('^', ecrange, sep = ''), DB$Data$EC)
    resDF <- DB$Data[ECs, ]
    FPs <- DB$FP[ECs]
    SIM <- list()
    for (itr in (1:nrow(resDF))) {
      SIM[itr] <- .calcSimilarity (fp_r, fp_p, FPs[[itr]]$FPR, FPs[[itr]]$FPP, reversible, algo, sim.method)
    }
    #SIM <- lapply(FPs, function(FP, fp_r, fp_p, reversible, algo, sim.method) {.calcSimilarity (fp_r, fp_p, FP$FPR, FP$FPP, reversible, algo, sim.method)}, fp_r, fp_p, reversible, algo, sim.method)
    
    resDF$'SIMILARITY' <- as.numeric(SIM)
    
    
    if (sort) {
      resDF <- resDF[order(-resDF$SIMILARITY),]
    }
    resDF
  }, error = function(err) {
    stop(err)
  })
  result
}