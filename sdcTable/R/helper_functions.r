## Wrapper function for pasting key-Variables
pasteStrVec <- function(strVec, nrVars, coll=NULL) {
  if(length(strVec) %% nrVars != 0)
    stop("Wrong Dimensions!\n")
  if ( is.null(coll) ) {
    return(cpp_myPaste(as.character(strVec), as.integer(nrVars)[1], NA))
  } else {
    return(cpp_myPaste(as.character(strVec), as.integer(nrVars)[1], as.character(coll[1])))
  }
}

# alternative to expand.grid (used for pasteStrVec!)
expand <- function(inputList, vector=TRUE) {
  uniques <- sapply(inputList, length)
  nrPoss <- prod(uniques)
  if ( vector == TRUE ) {
    out <- NULL
    for ( i in 1:length(inputList) ) {
      if ( i == 1 )
        out <- rep(inputList[[i]], nrPoss/length(inputList[[i]]))
      else
        out <- c(out, rep(inputList[[i]], each=prod(uniques[1:(i-1)]), nrPoss/length(rep(inputList[[i]], each=prod(uniques[1:(i-1)])))))
    }
  }
  else {
    out <- list()
    for ( i in 1:length(inputList) ) {
      if ( i == 1 )
        out[[i]] <- rep(inputList[[i]], nrPoss/length(inputList[[i]]))
      else
        out[[i]] <- rep(inputList[[i]], each=prod(uniques[1:(i-1)]), nrPoss/length(rep(inputList[[i]], each=prod(uniques[1:(i-1)]))))
    }
  }
  out
}

# returns a vector original size or str
mySplit <- function(strVec, keepIndices) {
  if ( min(keepIndices) < 1 | max(keepIndices) > nchar(strVec[1]) ) {
    stop("indices must be in 1:",nchar(strVec[1]),"!\n")
  }
  keepIndices <- unique(keepIndices)
  return(cpp_mySplit(as.character(strVec), as.integer(keepIndices)))
}

#strs <- rep(paste(LETTERS[1:6],collapse=""), 10000)
#system.time({
#  sapply(strs, mySplit, c(1,6))
#})

mySplitIndicesList <- function(strVec, keepList, coll="-") {
  u <- unlist(keepList)
  if ( min(u) < 1 | max(u) > nchar(strVec[1]) ) {
    stop("indices must be in 1:",nchar(strVec[1]),"!\n")
  }
  out <- list()
  for ( i in 1:length(keepList) ) {
    out[[i]] <- mySplit(strVec, keepList[[i]])
  }
  out <- cpp_myPaste(as.character(unlist(out)), as.integer(length(out)), coll)
}
# mySplitIndicesList("112233444", list(1:3, 5:6, 7:8))

# check ob 'x' ganzzahlig ist
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
  abs(x - round(x)) < tol
}

is.zero <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - 0) < tol
}

is.one <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - 1) < tol
}

# welche Variable soll als Branching_Variable verwendet werden?
getBranchingVariable <- function(sol, alreadyBranched, primSupps) {
  ind <- setdiff(1:length(sol), c(alreadyBranched, primSupps))
  branchVar <- ind[which.min(0.5 - sol[ind])]
  branchVar
}

my.Rglpk_solve_LP <- function(obj, mat, dir, rhs, types = NULL, max = FALSE, bounds = NULL, verbose = FALSE) {
  get.simpleTriplet <- sdcTable::get.simpleTriplet
  if ( !identical(max, TRUE) && !identical(max, FALSE) ) {
    stop("'Argument 'max' must be either TRUE or FALSE!\n")
  }
  direction_of_optimization <- as.integer(max)
  if ( !identical(verbose, TRUE) && !identical(verbose, FALSE) ) {
    stop("'Argument 'verbose' must be either TRUE or FALSE.\n")
  }
  if ( !class(mat) == "simpleTriplet" ) {
    stop("argument 'mat' must be of class 'simpleTriplet'\n")
  }

  verb <- as.integer(verbose)
  n_of_constraints <- length(dir)
  direction_of_constraints <- match(dir, c("<", "<=", ">", ">=", "=="))
  if ( any(is.na(direction_of_constraints)) ) {
    stop("Argument 'dir' must be either '<', '<=', '>', '>=' or '=='.\n")
  }
  n_of_objective_vars <- length(obj)
  if ( is.null(types) ) {
    types <- "C"
  }
  if ( any(is.na(match(types, c("I", "B", "C"), nomatch = NA))) ) {
    stop("'types' must be either 'B', 'C' or 'I'.\n")
  }
  types <- rep(types, length.out = n_of_objective_vars)
  integers <- types == "I"
  binaries <- types == "B"
  is_integer <- any(binaries | integers)
  bounds <- as.glp_bounds(as.list(bounds), n_of_objective_vars)
  x <- glp_call_interface(
    obj,
    n_of_objective_vars,
    get.simpleTriplet(mat, type="rowInd", input=list()),
    get.simpleTriplet(mat, type="colInd", input=list()),
    get.simpleTriplet(mat, type="values", input=list()),
    length(get.simpleTriplet(mat, type="values", input=list())),
    rhs, direction_of_constraints, n_of_constraints, is_integer,
    integers, binaries,
    direction_of_optimization,
    bounds[,1L],
    bounds[,2L],
    bounds[,3L],
    verb
  )
  solution <- x$lp_objective_vars_values
  solution[integers | binaries] <- round(solution[integers | binaries])
  status <- as.integer(x$lp_status != 5L)
  list(optimum = sum(solution * obj), solution = solution, status = status)
}
environment(my.Rglpk_solve_LP) <- environment(Rglpk_solve_LP)

# calculates years, weeks, days, hours, minutes and seconds from integer number
# secs: count of elapsed seconds (proc.time()[3])
# returns also a formatted string
formatTime <- function(secs){
  time.vec <- rep(NA, 6)
  names(time.vec) <- c('seconds', 'minutes','hours', 'days','weeks','years')

  secs <- ceiling(secs)

  time.vec['years'] <- floor(secs / (60*60*24*7*52))
  if ( time.vec['years'] > 0 ) {
    secs <- secs - (time.vec['years'] * (60*60*24*7*52))
  }

  time.vec['weeks'] <- floor(secs / (60*60*24*7))
  if ( time.vec['weeks'] > 0 ) {
    secs <- secs - (time.vec['weeks'] * (60*60*24*7))
  }

  time.vec['days'] <- floor(secs / (60*60*24))
  if ( time.vec['days'] > 0 ) {
    secs <- secs - time.vec['days']*(60*60*24)
  }

  time.vec['hours'] <- floor(secs / (60*60))
  if ( time.vec['hours'] > 0 ) {
    secs <- secs - time.vec['hours']*(60*60)
  }

  time.vec['minutes'] <- floor(secs / (60))
  if ( time.vec['minutes'] > 0 ) {
    secs <- secs - time.vec['minutes']*(60)
  }

  time.vec['seconds'] <- secs
  time.vec <- rev(time.vec)

  # time str #
  x <- time.vec[time.vec!=0]
  shortNames <- sapply(1:length(x), function(y) { substr(names(x)[y], 1, nchar(names(x)[y])-1)  } )

  time.str <- NULL
  for ( i in seq_along(names(x))) {

    if ( length(x) == 1 ) {
      if ( x[i] > 1 ) {
        time.str <- paste(time.str, x[i], " ", names(x[i]), sep="")
      } else {
        time.str <- paste(time.str, x[i], " ", shortNames[i], sep="")
      }
    }
    else {

      if ( names(x)[i]=="seconds") {
        if ( x[i] > 1 ) {
          time.str <- paste(time.str, "and", x[i], names(x[i]), sep=" ")
        } else {
          time.str <- paste(time.str, "and", x[i], shortNames[i], sep=" ")
        }

      } else {
        if ( x[i] > 1 ) {
          time.str <- paste(time.str, x[i], " ", names(x[i]), sep="")
        } else {
          time.str <- paste(time.str, x[i], " ", shortNames[i], sep="")
        }

        if ( i != length(x)-1 ) {
          time.str <- paste(time.str,", ", sep="")
        }
      }
    }
  }
  return(list(time.vec=time.vec, time.str=time.str))
}

# create default parameter objects suitable for primary|secondary suppression
# if selection == 'control.primary': set arguments suitable for primary suppression
# if selection == 'control.secondary': set arguments suitable for secondary suppression
genParaObj <- function(selection, ...) {
  controlPrimary <- function(...) {
    ### setDefaults ###
    paraObj <- list()

    # freq.rule
    paraObj$maxN <- 3
    paraObj$allowZeros <- FALSE

    # p-percent rule
    paraObj$p <- 80

    # n,k rule
    paraObj$n <- 2
    paraObj$k <- 85

    # pq-rule
    paraObj$pq <- c(25, 50)

    paraObj$numVarInd <- NA

    newPara <- list(...)

    indexNumVarIndices <- which(names(newPara) == "numVarIndices")
    if ( length(indexNumVarIndices) == 0 ) {
      stop("genPara (type=='control.primary'): parameter 'numVarIndices' must be specified\n")
    } else {
      numVarIndices <- newPara[[indexNumVarIndices]]
    }

    for ( i in seq_along(newPara) ) {
      m <- match(names(newPara)[i], names(paraObj))
      if ( !is.na(m) ) {
        paraObj[[m]] <- newPara[[i]]
      }
    }

    #if ( any(sapply(paraObj, length)!=1) ) {
    # stop("genPara (type=='control.primary'): arguments for primary suppression are not valid!\n")
    #}
    if ( !is.logical(paraObj$allowZeros) ) {
      stop("genPara (type=='control.primary'): argument 'allowZeros' must be logical!\n")
    }
    if ( !all(c(is.numeric(paraObj$maxN), is.numeric(paraObj$p), is.numeric(paraObj$n), is.numeric(paraObj$k))) ) {
      stop("genPara (type=='control.primary'): arguments 'maxN', 'p', 'n' and 'k' must be numeric!\n")
    }
    if ( length(paraObj$pq) != 2 ) {
      stop("genPara (type=='control.primary'): length of argument 'pq' must equal 2!\n")
    }
    if ( paraObj$k < 1 | paraObj$k >= 100) {
      stop("genPara (type=='control.primary'): argument 'k' must be >= 1 and < 100!\n")
    }
    if ( paraObj$p < 1 | paraObj$p >= 100) {
      stop("genPara (type=='control.primary'): argument p must be >= 1 and < 100!\n")
    }
    if ( paraObj$pq[1] < 1 | paraObj$pq[1] >= 100) {
      stop("genPara (type=='control.primary'): argument 'p' of 'pq' must be >= 1 and < 100!\n")
    }
    if ( paraObj$pq[2] < 1 | paraObj$pq[2] >= 100) {
      stop("genPara (type=='control.primary'): argument 'q' of 'pq' must be >= 1 and < 100!\n")
    }
    if ( paraObj$pq[1] >= paraObj$pq[2] ) {
      stop("genPara (type=='control.primary'): argument 'p' of 'pq' must be < argument 'q' of 'pq'\n")
    }
    if ( !is.na(paraObj$numVarInd) ) {
      if ( !paraObj$numVarInd %in% 1:length(numVarIndices) ) {
        stop("genPara (type=='control.primary'): argument 'numVarInd' must be >= 1 and <=",length(numVarIndices),"!\n")
      }
    }
    return(paraObj)
  }

  ### create a parameter list with (...) changing the default-values -> used in protectTable()
  controlSecondary <- function(...) {
    ### setDefaults ###
    paraObj <- list()

    # general parameter
    paraObj$method <- NA
    paraObj$verbose <- FALSE
    paraObj$save <- FALSE
    paraObj$solver <- "glpk"

    # HITAS|OPT - parameter
    paraObj$maxIter <- 10
    paraObj$timeLimit <- NULL
    paraObj$maxVars <- NULL
    paraObj$fastSolution <- FALSE
    paraObj$fixVariables <- TRUE
    paraObj$approxPerc <- 10
    paraObj$useC <- FALSE

    # HYPERCUBE - parameter
    paraObj$protectionLevel <- 80
    paraObj$suppMethod <- "minSupps"
    paraObj$suppAdditionalQuader <- FALSE

    # protectLinkedTables
    paraObj$maxIter <- 5

    newPara <- list(...)
    for ( i in seq_along(newPara) ) {
      m <- match(names(newPara)[i], names(paraObj))
      if ( !is.na(m) ) {
        paraObj[[m]] <- newPara[[i]]
      }
    }

    ### checks
    if ( any(sapply(paraObj, length)!=1) ) {
      stop("genPara (type=='control.secondary'): arguments controlObj for sdc-procedure are not valid!\n")
    }
    if ( !all(c(is.numeric(paraObj$maxIter), is.numeric(paraObj$approxPerc), is.numeric(paraObj$protectionLevel), is.numeric(paraObj$maxIter))) ) {
      stop("genPara (type=='control.secondary'): arguments 'maxIter', 'maxIter', 'protectionLevel' and 'maxIter' must be numeric!\n")
    }
    if ( !all(c(is.logical(paraObj$verbose), is.logical(paraObj$save), is.logical(paraObj$fastSolution), is.logical(paraObj$fixVariables), is.logical(paraObj$suppAdditionalQuader))) ) {
      stop("genPara (type=='control.secondary'): arguments 'verbose', 'save', 'fastSolution' 'fixVariables' and 'suppAdditionalQuader' must be numeric!\n")
    }
    if ( !is.null(paraObj$timeLimit) && !paraObj$timeLimit %in% 1:3000 ) {
      stop("genPara (type=='control.secondary'): argument 'timeLimit' must be >= 1 and <= 3000 minutes!\n")
    }
    if ( !length(paraObj$approxPerc) & !paraObj$approxPerc %in% 1:100 ) {
      stop("genPara (type=='control.secondary'): argument 'approxPerc' must be >= 1 and <= 100!\n")
    }
    if ( !paraObj$method %in% c('SIMPLEHEURISTIC', 'HITAS', 'HYPERCUBE', 'OPT') ) {
      stop("genPara (type=='control.secondary'): 'method' must be either 'SIMPLEHEURISTIC', 'HITAS', 'HYPERCUBE' or 'OPT'!\n")
    }
    if ( !paraObj$suppMethod %in% c('minSupps', 'minSum', 'minSumLogs') ) {
      stop("genPara (type=='control.secondary'): 'suppMethod' must be either 'minSupps', 'minSum' or 'minSumLogs'!\n")
    }
    return(paraObj)
  }

  if ( !selection %in% c('control.primary', 'control.secondary') ) {
    stop("genPara:: argument 'selection' must be either 'control.primary' or 'control.secondary'!\n")
  }

  if ( selection == 'control.primary' ) {
    paraObj <- controlPrimary(...)
  }
  if ( selection == 'control.secondary' ) {
    paraObj <- controlSecondary(...)
  }
  return(paraObj)
}

# convert simple triplet to matrix
st_to_mat <- function(x) {
  n.rows <- g_nr_rows(x)
  n.cols <- g_nr_cols(x)
  M <- matrix(0, nrow=n.rows, ncol=n.cols)

  i.x <- g_row_ind(x)
  j.x <- g_col_ind(x)
  v.x <- g_values(x)
  for ( i in 1:g_nr_cells(x) ) {
    M[i.x[i], j.x[i]] <- v.x[i]
  }
  # matrizen from attackers problem are transposed -> switch!
  return(t(M))
}

csp_cpp <- function(sdcProblem, attackonly=FALSE, verbose) {
  pI <- g_problemInstance(sdcProblem)
  dimInfo <- g_dimInfo(sdcProblem)
  aProb <- c_make_att_prob(input=list(objectA=pI, objectB=dimInfo))$aProb

  # already suppressed cells
  ind_prim <- as.integer(sort(c(g_primSupps(pI), g_secondSupps(pI))))
  len_prim <- as.integer(length(ind_prim))
  bounds_min <- bounds_max <- rep(0, len_prim)

  ind_fixed <- as.integer(g_forcedCells(pI))
  len_fixed <- as.integer(length(ind_fixed))

  aProb <- c_make_att_prob(input=list(objectA=pI, objectB=dimInfo))$aProb
  attProbM <- init.simpleTriplet("simpleTriplet", input=list(mat=st_to_mat(aProb@constraints)))

  ia <- as.integer(c(0, g_row_ind(attProbM)))
  ja <- as.integer(c(0, g_col_ind(attProbM)))
  ar <- as.double(c(0, g_values(attProbM)))

  cells_mat <- as.integer(length(ia))
  nr_vars <- as.integer(g_nr_cols(attProbM))
  nr_rows <- as.integer(g_nr_rows(attProbM))

  vals <- as.integer(g_freq(pI))

  lb <- as.double(g_lb(pI))
  ub <- as.double(g_ub(pI))

  LPL <- as.integer(g_LPL(pI))
  UPL <- as.integer(g_UPL(pI))
  SPL <- as.integer(g_SPL(pI))

  if ( attackonly == TRUE ) {
    attackonly <- as.integer(1)
  } else {
    attackonly <- as.integer(0)
  }
  final_pattern <- as.integer(rep(0, length(vals)))
  time.start <- proc.time()
  res <- .C("csp",
    ind_prim=ind_prim,
    len_prim=len_prim,
    bounds_min=bounds_min,
    bounds_max=bounds_max,
    ind_fixed=ind_fixed,
    len_fixed=len_fixed,
    ia=ia,
    ja=ja,
    ar=ar,
    cells_mat=cells_mat,
    nr_vars=nr_vars,
    nr_rows=nr_rows,
    vals=vals,
    lb=lb, ub=ub,
    LPL=LPL,
    UPL=UPL,
    SPL=SPL,
    final_pattern=final_pattern,
    attackonly=attackonly,
    verbose=as.integer(verbose),
    is_ok=0L
  )

  if ( attackonly ) {
    df <- data.frame(prim_supps=res$ind_prim, val=res$vals[res$ind_prim], bounds_low=res$bounds_min, bounds_up=res$bounds_max)
    df$protected <- df$bounds_low <= df$val - LPL[df$prim_supps]  &
    df$bounds_up >=  df$val + UPL[df$prim_supps] &
      df$bounds_up - df$bounds_low >= SPL[df$prim_supps]

    if ( length(g_secondSupps(pI)) > 0 ) {
      index <- g_primSupps(pI)
      df <- df[which(df$prim_supps %in% index),]
    }
    return(df)
  } else {
    if  ( res$is_ok != 0 ) {
      return(NULL)
    } else {
      nr_vars <- g_nrVars(g_problemInstance(sdcProblem))
      status_new <- rep("s", nr_vars)
      status_new[res$final_pattern!=0] <- "x"
      status_new[ind_prim] <- "u"
      if ( length(g_secondSupps(pI)) > 0 ) {
        status_new[g_secondSupps(pI)] <- "x"
      }
      if ( length(ind_fixed) > 0 ) {
        status_new[ind_fixed] <- "z"
      }

      pI <- g_problemInstance(sdcProblem)
      s_sdcStatus(pI) <- list(index=1:nr_vars, vals=status_new)
      s_problemInstance(sdcProblem) <- pI

      time.el <- g_elapsedTime(sdcProblem)+(proc.time()-time.start)[3]
      s_elapsedTime(sdcProblem) <- time.el
      s_indicesDealtWith(sdcProblem) <- 1:nr_vars
      return(sdcProblem)
    }
  }
}

