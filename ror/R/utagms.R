MINEPS <- 1E-10

utagms.strong.necessary <- function(performances, strongPrefs=NULL, weakPrefs=NULL, indifPrefs=NULL, strictVF=FALSE) {
  res <- !t(utagms(performances, strongPrefs, weakPrefs, indifPrefs, strictVF, necessary=FALSE))
  class(res) <- "binary.relation"
  return(res)
}

plot.binary.relation <- function(x, layout=igraph::layout.auto, ...) {
  igraph::plot.igraph(graph.adjacency(x), layout=layout, ...)  
}

utagms <- function(performances, strongPrefs = NULL, weakPrefs = NULL, indifPrefs = NULL, necessary=TRUE, strictVF=FALSE) {
  rel <- matrix(nrow=nrow(performances), ncol=nrow(performances))
  
  if (!checkConsistency(performances, necessary, strictVF, strongPrefs, weakPrefs, indifPrefs)) {
    stop("Model infeasible")
  }

  for (i in 1:nrow(rel)) {
    for(j in 1:nrow(rel)) {
      rel[i,j] = checkRelation(performances, i, j, necessary=necessary, strictVF=strictVF,
           strongPrefs=strongPrefs, weakPrefs=weakPrefs, indifPrefs=indifPrefs)
    }
  }
  if (!is.null(rownames(performances))) {
    rownames(rel) <- rownames(performances)
    colnames(rel) <- rownames(performances)
  }

  class(rel) <- "binary.relation"
  return(rel)
}

checkConsistency <- function(perf, necessary, strictVF, strongPrefs, weakPrefs, indifPrefs) {
  ## check vars
  stopifnot(is.logical(necessary))
  stopifnot(is.logical(strictVF))
  
  altVars <- buildAltVariableMatrix(perf)
  baseModel <- buildBaseLPModel(perf, strictVF=strictVF, strongPrefs=strongPrefs,
                                weakPrefs=weakPrefs, indifPrefs=indifPrefs)

  ret <- solveModel(perf, baseModel)

  return(ret$status$code == 0 && ret$objval >= MINEPS)
}

solveModel <- function(perf, model) {
  obj <- L_objective(buildObjectiveFunction(perf))
  roiConst <- L_constraint(model$lhs, model$dir, model$rhs)
  lp <- OP(objective=obj, constraints=roiConst, maximum=TRUE)
  ROI_solve(lp, .solver)
}

checkRelation <- function(perf, a, b, necessary, strictVF, strongPrefs, weakPrefs, indifPrefs) {
  ## check vars
  stopifnot(is.logical(necessary))
  stopifnot(is.logical(strictVF))
  if (a == b) {
    return(TRUE)
  }
  altVars <- buildAltVariableMatrix(perf)  
  baseModel <- buildBaseLPModel(perf, strictVF=strictVF, strongPrefs=strongPrefs, weakPrefs=weakPrefs, indifPrefs=indifPrefs)

  addConst <- c()
  if (necessary == TRUE) {
    addConst <- buildStrongPreferenceConstraint(b, a, altVars)
  } else { ## possible
    addConst <- buildWeakPreferenceConstraint(a, b, altVars)
  }
  allConst <- combineConstraints(baseModel, addConst)

  ret <- solveModel(perf, allConst)
#  cat("a", a, "b", b, "code", ret$status$code, "objval", ret$objval, "\n")

  if (necessary == TRUE) {
    return(ret$status$code != 0 || ret$objval < MINEPS)
  } else { # possible
    return(ret$status$code == 0 && ret$objval >= MINEPS)
  }
}

buildObjectiveFunction <- function(perf) {
  levels <- getLevels(perf)
  nrVars <- getNrVars(levels)

  lhs <- rep(0, nrVars)
  lhs[length(lhs)] = 1
  return(lhs)
}


## perf: the performance matrix
## strictVF = TRUE -> value functions strictly increasing (instead of monotonous increasing)
## *prefs: an n x 2 matrix, where each row (a, b) means
## that a is [strongly or weakly preferred, or indifferent] to b.
buildBaseLPModel <- function(perf, strictVF, strongPrefs, weakPrefs, indifPrefs) {
  altVars <- buildAltVariableMatrix(perf)

  c1 <- buildMonotonousConstraints(perf, strictVF=strictVF)
  c2 <- buildFirstLevelZeroConstraints(perf)
  c3 <- buildBestLevelsAddToUnityConstraint(perf)
  c4 <- buildAllVariablesLessThan1Constraint(perf)
  c5 <- buildEpsilonStrictlyPositiveConstraint(perf)

  allConst <- combineConstraints(c1, c2, c3, c4, c5)

  if (is.matrix(strongPrefs)) {
    for (i in 1:nrow(strongPrefs)) {
      prefConst <- buildStrongPreferenceConstraint(strongPrefs[i,1], strongPrefs[i,2], altVars)
      allConst <- combineConstraints(allConst, prefConst);
    }
  }
  if (is.matrix(weakPrefs)) {
    for (i in 1:nrow(weakPrefs)) {
      prefConst <- buildWeakPreferenceConstraint(weakPrefs[i,1], weakPrefs[i,2], altVars)
      allConst <- combineConstraints(allConst, prefConst);
    }
  }
  if (is.matrix(indifPrefs)) {
    for (i in 1:nrow(indifPrefs)) {
      prefConst <- buildIndifPreferenceConstraint(indifPrefs[i,1], indifPrefs[i,2], altVars)
      allConst <- combineConstraints(allConst, prefConst);
    }
  }

  return(allConst)
}

buildStrongPreferenceConstraint <- function(a, b, altVars) {
  nrVars <- dim(altVars)[2]

  lhs <- altVars[a,]
  lhs[length(lhs)] = -1
  lhs <- lhs - altVars[b,]
  
  return(list(lhs=lhs, dir=">=", rhs=0))
}

buildIndifPreferenceConstraint <- function(a, b, altVars) {
  lhs <- altVars[a,]
  lhs <- lhs - altVars[b,]
  
  return(list(lhs=lhs, dir="==", rhs=0))
}

buildWeakPreferenceConstraint <- function(a, b, altVars) {
  lhs <- altVars[a,]
  lhs <- lhs - altVars[b,]
  
  return(list(lhs=lhs, dir=">=", rhs=0))
}

buildEpsilonStrictlyPositiveConstraint <- function(perf) {
  levels <- getLevels(perf)
  nrVars <- getNrVars(levels)

  lhs <- rep(0, nrVars)
  lhs[length(lhs)] = 1
  return(list(lhs=lhs, dir=">=", rhs=MINEPS))
}

buildAllVariablesLessThan1Constraint <- function(perf) {
  levels <- getLevels(perf)
  nrVars <- getNrVars(levels)

  lhs <- diag(nrVars)

  return(list(lhs=lhs, dir=rep("<=", nrVars), rhs=rep(1, nrVars)))
}

buildBestLevelsAddToUnityConstraint <- function(perf) {
  levels <- getLevels(perf)
  offsets <- getOffsets(levels)
  nrVars <- getNrVars(levels)

  lhs <- rep(0, nrVars)
  ind <- c((offsets-1)[-1], nrVars-1)
  lhs[ind] = 1
  return(list(lhs=lhs, dir="==", rhs=1))
}

buildFirstLevelZeroConstraints <- function(perf) {
  levels <- getLevels(perf)
  offsets <- getOffsets(levels)
  nrVars <- getNrVars(levels)
  
  res <- matrix(0, nrow=length(offsets),ncol=nrVars)

  for (i in seq(1:length(offsets))) {
    res[i,offsets[i]] = 1
  }

  return(list(lhs=res,dir=rep("==", length(offsets)),rhs=rep(0,length(offsets))))
}

buildMonotonousConstraints <- function(perf, strictVF=FALSE) {
  
  stopifnot(is.logical(strictVF))
  
  levels <- getLevels(perf)
  offsets <- getOffsets(levels)
  nrVars <- getNrVars(levels)

  res <- c()
  
  for (i in seq(1:length(levels))) {
    for (j in seq(1:(length(levels[[i]])-1))) {
      index <- offsets[i] + j - 1
      lhs <- array(0, dim=nrVars)
      lhs[index] <- 1
      lhs[index+1] <- -1
      if (strictVF == TRUE) {
        lhs[length(lhs)] = 1
      }
      res <- rbind(res, lhs)
    }
  }

  return(list(lhs=res, dir=rep("<=", nrow(res)), rhs=rep(0, nrow(res))))
}

getLevels <- function(perf) {
  res <- list()
  for (i in 1:ncol(perf)) {
    res[[i]] <- sort(unique(perf[,i]))
  }
  return(res)
}

getNrVars <- function(levels) {
  return(sum(as.numeric(lapply(levels, length))) + 1)
}

buildAltVariableMatrix <- function(perf) {
  levels <- getLevels(perf)

  offsets <- getOffsets(levels)
  nrAlts <- nrow(perf)
  nrCrit <- ncol(perf)

  nrVars <- getNrVars(levels)

  resMat = matrix(nrow=nrAlts,ncol=nrVars)
  
  for (i in seq(1:nrAlts)) {
    vec <- array(0, dim=nrVars)
    indices <- sapply(seq(1:nrCrit), function(x) {which(levels[[x]] == perf[i,x])})
    vec[indices + offsets - 1] = 1
    resMat[i,] = vec
  }
  return(resMat)
}

getOffsets <- function(levels) {
  x <- cumsum(lapply(levels, length))
  return(c(1, x[1:length(x)-1] + 1))
}
