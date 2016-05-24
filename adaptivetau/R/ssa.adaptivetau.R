ssa.adaptivetau <-
function(init.values, transitions, rateFunc, params, tf,
         jacobianFunc = NULL, maxTauFunc = NULL,
         deterministic = NULL, halting = NULL,
         relratechange=rep(1, length(init.values)),
         tl.params = NULL) {
  return(.Call('simAdaptiveTau', PACKAGE='adaptivetau',
               init.values, transitions,
               rateFunc, jacobianFunc,
               params, tf, deterministic, halting,
               relratechange, tl.params, maxTauFunc))
}

ssa.exact <-
function(init.values, transitions, rateFunc, params, tf) {
  return(.Call('simExact', PACKAGE='adaptivetau',
               init.values, transitions, rateFunc, params, tf))
}

ssa.maketrans <- function(variables, ...) {
  userTrans = list(...)
  if (length(userTrans) == 0) {
    stop("no transitions passed into ssa.maketrans!")
  }
  if (length(variables) == 1  &&  is.numeric(variables)) {
    numVariables = variables;
  } else if (is.character(variables)) {
    numVariables = length(variables)
  } else {
    stop("Cannot deduce number of variables -- ssa.maketrans requires ",
         "either a vector of variable names or the number of variables")
  }

  allTrans = vector("list", length = sum(sapply(userTrans, function(x) (
                                if (is.matrix(x)) ncol(x) else 1))));

  trI = 0
  for (i in 1:length(userTrans)) {
    x = userTrans[[i]]
    if (length(x) == 1  &&  is.na(x)) {
        trI = trI + 1
        allTrans[[trI]] = integer(0)
        next
    }
    idx = (1:(nrow(x) %/% 2))*2-1 #indices of variables
    mag = idx + 1                 #indices of magnitudes
    for (j in seq_len(ncol(x))) {
      trI = trI + 1
      if (is.numeric(x)) {
        if (any(x[idx,j] < 1  |  x[idx,j] > numVariables)) {
          stop("variable index outside valid range (1:numVariables)")
        }
        allTrans[[trI]] = structure(as.integer(x[mag,j]), names = x[idx,j])
      } else if (is.character(x)) {
        if (any(!(x[idx,j] %in% variables))) {
          stop("unknown variable(s): ",
               paste(x[idx,j][!(x[idx,j] %in% variables)], collapse=", "))
        }
        allTrans[[trI]] = structure(as.integer(x[mag,j]), names = x[idx,j])
      } else {
        stop("transitions passed to ssa.maketrans must be integer or character")
      }
    }
  }

  allTrans
}
