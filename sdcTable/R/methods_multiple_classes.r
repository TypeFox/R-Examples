#' @aliases calc.multiple,character,list-method
#' @rdname calc.multiple-method
setMethod(f='calc.multiple', signature=c('character', 'list'),
  definition=function(type, input) {
    .SD <- ID <- NULL
    if (!type %in% c('makePartitions', 'genMatMFull',
        'makeAttackerProblem', 'calcFullProblem') ) {
      stop("calc.multiple:: argument 'type' is not valid!\n")
    }

    if ( type == 'makePartitions' ) {
      return(c_make_partitions(input))
    }
    if ( type == 'genMatMFull' ) {
      return(c_gen_mat_m(input))
    }
    if ( type == 'makeAttackerProblem' ) {
      return(c_make_att_prob(input))
    }
    if ( type == 'calcFullProblem' ) {
      return(c_calc_full_prob(input))
    }
  }
)


setMethod("c_make_partitions", signature=c("list"), definition=function(input) {
  pI <- input$objectA
  dimInfoObj <- input$objectB
  dimInfo <- g_dim_info(dimInfoObj)
  strIDs <- g_strID(pI)

  ## create classes and groups
  tmpDat <- expand.grid(lapply(1:length(dimInfo), function(x) { 1:g_nr_levels(dimInfo[[x]]) } ))
  groups <- apply(tmpDat, 1, function(x) { paste(x, collapse="-")})
  classes <- apply(tmpDat, 1, sum)
  sortOrder <- order(classes)
  classes <- classes[sortOrder]
  classesUnique <- unique(classes)
  groups <- groups[sortOrder]
  splitGroups <- split(groups, classes)

  ## create tables for all classes and groups
  final <- list()
  final$groups <- as.list(groups)
  final$indices <- list()

  # default_codes and levels
  default_codes <- lapply(1:length(dimInfo), function(x) {
    g_default_codes(dimInfo[[x]])
  })
  dim_levels <- lapply(1:length(dimInfo), function(x) {
    g_levels(dimInfo[[x]])
  })

  # data.table to merge on
  df <- data.table(N=1:length(strIDs), strIDs=strIDs)
  setkey(df, strIDs)

  for ( i in 1:length(groups) ) {
    final$indices[[i]] <- list()
    levs <- as.integer(unlist(sapply(groups[[i]], strsplit, "-")))

    res <- list()
    for ( z in 1:length(dimInfo) ) {
      res[[z]] <- list()
      index <- which(g_levels(dimInfo[[z]]) %in% c(levs[z], levs[z]-1))
      codesDefault <- default_codes[[z]][index]
      if ( levs[z] == 1 ) {
        res[[z]] <- codesDefault
      } else {
        levOrig <- dim_levels[[z]][index]
        diffs <- c(0,diff(levOrig))
        checkInd <- which(diffs == 1)-1
        out <- data.frame(index=index, levOrig=levOrig, codesDefault=codesDefault, ind=NA)
        out$ind[checkInd] <- 1

        checkInd <- c(checkInd, length(index))
        splitVec <- rep(0, length(index))
        for ( j in 2:length(checkInd) ) {
          if ( j < length(checkInd) ) {
            splitVec[checkInd[j-1]:(checkInd[j]-1)] <- j-1
          } else {
            splitVec[checkInd[j-1]:(checkInd[j])] <- j-1
          }
        }
        spl <- split(index, splitVec)
        counter <- 1
        for ( k in 1:length(spl) ) {
          rowInd <- match(spl[[k]], out$index)
          tmp <- out[rowInd,]
          if ( any(tmp[,"levOrig"]==levs[z]) ) {
            tmp <- tmp[1:(max(which(tmp$levOrig==levs[z]))),]
            res[[z]][[length(res[[z]])+1]] <- sort(unique(as.character(tmp$codesDefault)))
          }
        }
      }
    }
    final$indices[[i]] <- list()
    combs <- expand.grid(lapply(1:length(res), function(x) {
      1:length(res[[x]])
    }))

    final$indices[[i]] <- list();
    length(final$indices[[i]]) <- nrow(combs)
    for ( m in 1:nrow(combs) ) {
      final.strIDs <- pasteStrVec(expand(lapply(1:ncol(combs), function(x) { res[[x]][[combs[m,x]]] })), ncol(combs))
      df2 <- data.table(strIDs=final.strIDs)
      setkey(df2, strIDs)
      final$indices[[i]][[m]] <- merge(df, df2)$N
    }
  }
  final$nrGroups <- length(groups)
  final$nrTables <- sum(sapply(1:final$nrGroups, function(x) { length(final$indices[[x]]) } ))
  return(final)
})

setMethod("c_gen_mat_m", signature=c("list"), definition=function(input) {
  x <- input$objectA
  y <- input$objectB

  levelObj <- g_dim_info(y)
  strID <- g_strID(x)
  nrVars <- length(levelObj)
  nrCells <- g_nrVars(x)
  freqs <- g_freq(x)

  constraintM <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, nrow=0, ncol=nrCells)))
  for ( i in 1:nrVars ) {
    lO <- levelObj[[i]]
    keepList <- lapply(g_str_info(y)[-i], function(k) {
      seq(k[1], k[2])
    })
    keepList2 <- lapply(g_str_info(y)[i], function(k) {
      seq(k[1], k[2])
    })
    f1 <- f2 <- mySplitIndicesList(strID, keepList2)

    if ( nrVars > 1 ) {
      f1 <- mySplitIndicesList(strID, keepList)
    }

    dimlO <- g_dims(lO)
    if ( length(unique(f2)) != 1 ) {
      dimInd <- sapply(1:length(dimlO), function(x) { identical( sort(unique(f2)), dimlO[[x]]) } )
      if ( sum(dimInd) == 0 ) {
        for ( j in 1:length(g_dims(lO)) ) {
          splitInd <- which(f2 %in% g_dims(lO)[[j]])
          spl <- split(splitInd, f1[splitInd])
          for ( z in 1:length(spl) ) {
            ind <- rep(1,length(spl[[z]]))
            ind[which.max(freqs[spl[[z]]])] <- -1
            if ( !is.zero(sum(freqs[spl[[z]]]*ind)) ) {
              stop("something went wrong!\n")
            }
            constraintM <- c_add_row(constraintM, input=list(index=spl[[z]], values=ind))
          }
        }
      } else {
        splitInd <- which(f2 %in% g_dims(lO)[[which(dimInd==TRUE)]])
        ## only 1 dimension
        if ( nrVars > 1 ) {
          spl <- split(splitInd, f1[splitInd])
        } else {
          spl <- split(splitInd, rep(1, length(splitInd)))
        }

        for ( z in 1:length(spl) ) {
          ind <- rep(1,length(spl[[z]]))
          ind[which.max(freqs[spl[[z]]])] <- -1
          if ( !is.zero(sum(freqs[spl[[z]]]*ind)) ) {
            stop("something went wrong! (z=",z," und names(spl)[z]='",names(spl)[z],")\n")
          }
          constraintM <- c_add_row(constraintM, input=list(index=spl[[z]], value=ind))
        }
      }
    }
  }
  return(constraintM)
})

setMethod("c_make_att_prob", signature=c("list"), definition=function(input) {
  x <- input$objectA
  y <- input$objectB
  nrVars <- g_nrVars(x)
  A <- c_gen_mat_m(input=list(objectA=x, objectB=y))

  ## calculating (logical) constraints for the master problem ##
  # idea: for each constraint at least 2 suppressions must
  # exist if one xi != 0! (http://www.eia.doe.gov/ices2/missing_papers.pdf)
  newCutsMaster <- init.cutList(type='empty', input=list(nrCols=nrVars))
  #xx <- lapply(1:g_nr_rows(A), function(x) {
  # cols <- g_col_ind(g_row(A, input=list(x)))
  # v <- rep(0, nrVars)
  # v[cols] <- c(1, rep(-1, length(cols)))
  # s_add_complete_constraint(newCutsMaster) <<- list(init.cutList(type='singleCut', input=list(vals=v, dir="<=", rhs=0)))
  #})
  ################################################################

  nrConstraints <- g_nr_rows(A)
  objective <- rep(0, length=2*nrVars+nrConstraints)
  z1 <- init.simpleTriplet(type='simpleTripletDiag', input=list(nrRows=nrVars, negative=FALSE))
  z2 <- init.simpleTriplet(type='simpleTripletDiag', input=list(nrRows=nrVars, negative=TRUE))
  z <- c_bind(object=z1, input=list(z2, bindRow=FALSE))
  A <- c_bind(object=z, input=list(g_transpose(A), bindRow=FALSE))
  direction <- rep("==", g_nr_rows(A))
  rhs <- rep(0, g_nr_rows(A))

  types <- rep("C", g_nr_cols(A))
  boundsLower <- list(ind=1:g_nr_cols(A), val=c(rep(0, 2*nrVars), rep(-Inf, nrConstraints)))
  boundsUpper <- list(ind=1:g_nr_cols(A), val=c(rep(Inf, 2*nrVars), rep(Inf,  nrConstraints)))

  aProb <- new("linProb",
    objective=objective,
    constraints=A,
    direction=direction,
    rhs=rhs,
    boundsLower=boundsLower,
    boundsUpper=boundsUpper,
    types=types)
  return(list(aProb=aProb, newCutsMaster=newCutsMaster))
})

setMethod("c_calc_full_prob", signature=c("list"), definition=function(input) {
  .SD <- ID <- NULL
  x <- input$objectA
  y <- input$objectB
  time.start <- proc.time()
  datO <- g_raw_data(x)
  dimObj <- g_dim_info(y)

  # we have to aggregate if we are dealing with microdata
  if ( g_is_microdata(x) ) {
    rawData <- datO[, lapply(.SD, sum, na.rm=TRUE), by=key(datO), .SDcols=setdiff(colnames(datO), key(datO))]
  } else {
    rawData <- copy(datO)
  }
  ind.dimvars <- g_dimvar_ind(x)
  ind.freq <- g_freqvar_ind(x)

  codes <- list(); length(codes) <- length(ind.dimvars)
  for ( i in 1:length(codes) ) {
    codes[[i]] <- rawData[[ind.dimvars[i]]]
    cDefault <- g_default_codes(dimObj[[i]])
    cOriginal <- g_original_codes(dimObj[[i]])
    cOriginalDups <- g_dups(dimObj[[i]])
    cOriginalDupsUp <- g_dups_up(dimObj[[i]])
    if ( all(codes[[i]] %in% c(cOriginal, cOriginalDups)) ) {
      mInd1 <- match(codes[[i]], cOriginalDups)
      mInd2 <- which(!is.na(mInd1))
      if ( length(mInd2) > 0 ) {
        codes[[i]][mInd2] <- cOriginalDupsUp[mInd1[mInd2]]
      }
      codes[[i]] <- c_match_default_codes(object=dimObj[[i]], input=rawData[[ind.dimvars[i]]])
    } else if ( all(codes[[i]] %in% cDefault) ) {
      # cat("no recoding necessary!\n")
    } else {
      stop("c_calc_full_prob:: recoding not possible!\n")
    }
  }

  ## calculate all possible combinations within the lowest levels of dim-vars
  ## if any combinations are missing (missing.codes), we have to set them to 0 later
  strID <- as.character(pasteStrVec(unlist(codes), length(codes)))
  exDims <- pasteStrVec(unlist(codes), length(codes))
  possDims <- sort(pasteStrVec(as.character(expand(lapply(dimObj, function(x) { 
    g_minimal_default_codes(x) 
  }), vector=TRUE)), length(dimObj)))
  missing.codes <- setdiff(possDims, exDims)

  ## fill the table
  nrIndexvars <- length(ind.dimvars)
  fullDims <- lapply(dimObj, g_dims)

  allCodes <- expand(lapply(dimObj, g_default_codes), vector=FALSE)
  fullTabObj <- data.table(ID=1:length(allCodes[[1]]))
  for ( i in 1:length(allCodes)) {
    fullTabObj[,colnames(rawData)[ind.dimvars][i]:=allCodes[[i]]]
  }
  setkeyv(fullTabObj, colnames(rawData)[ind.dimvars])
  fullTabObj[,ID:=NULL]

  ## revert rawData codes to default codes
  for ( j in seq_along(ind.dimvars) ) {
    v <- c_match_default_codes(object=dimObj[[j]], input=rawData[,get(names(dimObj)[j])])
    set(rawData, NULL, names(dimObj)[j], v)
  }
  setkeyv(rawData, colnames(rawData)[ind.dimvars])

  ## replace NAs in rawData by 0 (required for aggregation)
  cols <- colnames(rawData)[(length(dimObj)+1):ncol(rawData)]
  ind.na <- list(); length(ind.na) <- length(cols); k <- 1
  for ( j in cols ) {
    ind.na[[k]] <- which(is.na(rawData[[j]]))
    set(rawData, ind.na[[k]], j, 0)
    k <- k+1
  }; rm(k)

  ## merge minDat to fullDat
  fullTabObj <- merge(fullTabObj, rawData, all.x=TRUE)

  ## set missing combinations of lowest levels to 0
  ## problematic are all levels that should exist, but do not exist
  ## they are filled with 0 so that we can aggregate
  dim.vars <- colnames(fullTabObj)[ind.dimvars]
  strID <- apply(fullTabObj[,dim.vars,with=FALSE],1,paste0, collapse="")

  if ( length(missing.codes) > 0 ) {
    index <- which(strID%in%missing.codes)
    for ( i in 1:length(cols) ) {
      set(fullTabObj, index, cols[i], 0)
    }
  }

  ## fill up missing dimensions
  not.finished <- TRUE
  
  # which indexvars have any hierarchy (not just the total?)
  # these indiecs specify the dim-variables we loop over
  useInds <- which(sapply(y@dimInfo, function(x) {
    length(x@codesOriginal)>1
  }))  
  while ( not.finished ) {
    cols <- (nrIndexvars+1):ncol(fullTabObj)
    col.names <- colnames(fullTabObj)[cols]
    for ( i in useInds ) {
      if ( length(dim.vars) > 1 ) {
        setkeyv(fullTabObj, dim.vars[-i])
      } else {
        setkeyv(fullTabObj, dim.vars[1])
      }

      cur.dim <- dimObj[[i]]@dims
      for ( j in length(cur.dim):1 ) {
        cur.levs <-  cur.dim[[j]]
        out <- fullTabObj[fullTabObj[[ind.dimvars[i]]] %in% cur.levs[-1],]
        if ( length(dim.vars)==1 ) {
          out <- out[,lapply(.SD,sum), .SDcols=col.names]
        } else {
          out <- out[,lapply(.SD,sum), .SDcols=col.names, by=key(out)]
        }
        row.ind <- which(fullTabObj[[ind.dimvars[i]]] == cur.levs[1])
        for ( j in col.names ) {
          set(fullTabObj, row.ind, j, out[[j]])
        }
      }
    }
    if ( !is.na(fullTabObj[1,ind.freq,with=FALSE]) ) {
      not.finished <- FALSE
    }
  }

  nrV <- nrow(fullTabObj)
  f <- fullTabObj[[ind.freq]]
  strID <- apply(fullTabObj[,dim.vars,with=FALSE],1,paste0, collapse="")
  w <- numVarsList <- NULL
  w.ind <- g_weightvar_ind(x)
  if ( !is.null(w.ind) ) {
    w <- fullTabObj[[w.ind]]
  }
  n.ind <- g_numvar_ind(x)
  if ( !is.null(n.ind) ) {
    numVarsList <- list(); length(numVarsList) <- length(n.ind)
    for ( n in 1:length(n.ind) ) {
      numVarsList[[n]] <- fullTabObj[[n.ind[n]]]
    }
  }

  if ( length(n.ind) > 0 ) {
    names(numVarsList) <- colnames(g_raw_data(x))[n.ind]
  }

  ## replace 0 in rawData by NA if they have been replaced earlier
  for ( i in 1:length(ind.na) ) {
    if ( length(ind.na[[i]]) > 0 ) {
      set(rawData, ind.na[[i]], cols[i], NA)
    }
  }
  s_raw_data(x) <- list(datO)
  problemInstance <- new("problemInstance",
    strID=strID,
    Freq=f,
    w=w,
    numVars=numVarsList,
    lb=rep(0, nrV),
    ub=pmax(2*f, 5),
    LPL=rep(1, nrV),
    UPL=rep(1, nrV),
    SPL=rep(0, nrV),
    sdcStatus=rep("s", nrV)
  )
  partition <- c_make_partitions(input=list(objectA=problemInstance, objectB=y))
  sdcProblem <- new("sdcProblem",
    dataObj=x,
    dimInfo=y,
    problemInstance=problemInstance,
    partition=partition,
    startI=1,
    startJ=1,
    indicesDealtWith=NULL,
    elapsedTime=(proc.time()-time.start)[3]
  )
  return(sdcProblem)
})
