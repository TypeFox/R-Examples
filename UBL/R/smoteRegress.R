## ===================================================
## Creating a SMOTEd training sample for regression problems
# 
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae), ]
#   C.perc = list(0.1, 8) 
#   mysmote.alg <- SmoteRegress(a7~., clean.algae, dist = "HEOM",
#                               C.perc = C.perc)
#   smoteBal.alg <- SmoteRegress(a7~., clean.algae, dist = "HEOM",
#                               C.perc = "balance")
#   smoteExt.alg <- SmoteRegress(a7~., clean.algae, dist = "HEOM",
#                               C.perc = "extreme")
# 
#   ir<- iris[-c(95:130), ]
#   mysmote.iris <- SmoteRegress(Sepal.Width~., ir, dist = "HEOM",
#                                C.perc = list(0.5,2.5))
#   mysmote.iris <- SmoteRegress(Sepal.Width~., ir, dist = "HEOM",
#                                C.perc = list(0.2,4), thr.rel = 0.8)
#   smoteBalan.iris <- SmoteRegress(Sepal.Width~., ir, dist = "HEOM",
#                                C.perc = "balance")
#   smoteExtre.iris <- SmoteRegress(Sepal.Width~., ir, dist = "HEOM",
#                                C.perc = "extreme")
# 
#   rel <- matrix(0, ncol = 3, nrow = 0)
#   rel <- rbind(rel, c(2, 1, 0))
#   rel <- rbind(rel, c(3, 0, 0))
#   rel <- rbind(rel, c(4, 1, 0))
#
#   sP.ir <- SmoteRegress(Sepal.Width~., ir, dist = "HEOM", rel = rel,
#                        C.perc = list(4,0.5,4))
# 
# L. Torgo, Jun 2008
# P. Branco, Mar, Apr 2015 Apr 2016
# ---------------------------------------------------
SmoteRegress <- function(form, dat, rel = "auto", thr.rel = 0.5,
                         C.perc = "balance", k = 5, repl = FALSE,
                         dist = "Euclidean", p = 2)
  
  # INPUTS:
  # form    a model formula
  # dat    the original training set (with the unbalanced distribution)
  # rel     is the relevance determined automatically (default: "auto") 
  #         or provided by the user through a matrix. See examples.
  # thr.rel is the relevance threshold above which a case is considered
  #         as belonging to the rare "class"
  # C.perc is a list containing the percentage of under- or/and 
  #         over-sampling to apply to each "class" obtained with the threshold.
  #         The over-sampling percentage means that the examples above the 
  #         threshold are increased by this percentage. The under sampling
  #         percentage means that the normal cases (cases below the threshold)
  #         are under-sampled by this percentage. Alternatively it may be
  #         "balance" or "extreme", cases where the sampling percentages
  #         are automatically estimated.
  # k       is the number of neighbors to consider as the pool from where
  #         the new synthetic examples are generated
  # repl    is it allowed to perform sampling with replacement
  # dist    is the distance measure to be used (defaults to "Euclidean")
  # p       is a parameter used when a p-norm is computed
{
 
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }
  
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  
  if (tgt < ncol(dat)) {
    orig.order <- colnames(dat)
    cols <- 1:ncol(dat)
    cols[c(tgt, ncol(dat))] <- cols[c(ncol(dat), tgt)]
    dat <- dat[, cols]
  }
  if (is.na(thr.rel)) {
    stop("Future work!")
  }
  

  y <- dat[, ncol(dat)]
  attr(y, "names") <- rownames(dat)
  s.y <- sort(y)
  
  if (is.matrix(rel)) { 
    pc <- phi.control(y, method = "range", control.pts = rel)
  } else if (is.list(rel)) { 
    pc <- rel
  } else if (rel == "auto") {
    pc <- phi.control(y, method = "extremes")
  } else {# handle other relevance functions and not using the threshold!
    stop("future work!")
  }
  
  temp <- y.relev <- phi(s.y, pc)
  if (!length(which(temp < 1))) {
    stop("All the points have relevance 1. 
         Please, redefine your relevance function!")
  }
  if (!length(which(temp > 0))) {
    stop("All the points have relevance 0. 
         Please, redefine your relevance function!")
  }
  temp[which(y.relev > thr.rel)] <- -temp[which(y.relev > thr.rel)]
  bumps <- c()
  for (i in 1:(length(y) - 1)) { 
    if (temp[i] * temp[i + 1] < 0) bumps <- c(bumps, i) 
  }
  nbump <- length(bumps) + 1 # number of different "classes"
  
  # collect the indexes in each "class"
    obs.ind <- as.list(rep(NA, nbump))
    last <- 1
    for (i in 1:length(bumps)) {
      obs.ind[[i]] <- s.y[last:bumps[i]]
      last <- bumps[i] + 1
    }
    obs.ind[[nbump]] <- s.y[last:length(s.y)]

  newdata <- data.frame()
  
  if (is.list(C.perc)) {
    if (length(C.perc) != nbump){
      stop("The percentages provided must be the same length as the number
           of bumps!")
    }
  } else if (C.perc == "balance") {
    # estimate the percentages of over/under sampling
    B <- round(nrow(dat)/nbump, 0)
    C.perc <- B/sapply(obs.ind, length)        
  } else if (C.perc == "extreme") {
    B <- round(nrow(dat)/nbump, 0)
    rescale <- nbump * B/sum(B^2/sapply(obs.ind, length))
    obj <- round((B^2/sapply(obs.ind, length)) * rescale, 2)
    C.perc <- round(obj/sapply(obs.ind, length), 1)
  }
  
  for (i in 1:nbump) {
    if (C.perc[[i]] == 1) {
      newdata <- rbind(newdata, dat[names(obs.ind[[i]]), ])
    } else if (C.perc[[i]] > 1) {
      newExs <- Smote.exsRegress(dat[names(obs.ind[[i]]), ],
                                 ncol(dat),
                                 C.perc[[i]],
                                 k,
                                 dist,
                                 p)
      # add original rare examples and synthetic generated examples
      newdata <- rbind(newdata, newExs, dat[names(obs.ind[[i]]), ])
    } else if (C.perc[[i]] < 1) {
      sel.maj <- sample(1:length(obs.ind[[i]]),
                        as.integer(C.perc[[i]] * length(obs.ind[[i]])),
                        replace = repl)
      newdata <- rbind(newdata, dat[names(obs.ind[[i]][sel.maj]), ])
    }
  }
  
  if (tgt < ncol(dat)) {
    newdata <- newdata[, cols]
    dat <- dat[, cols]
  }
  
  newdata
}



# ===================================================
# Obtain a set of smoted examples for a set of rare cases.
#
# L. Torgo, Jun 2008
# P.Branco, Mar 2015 Apr 2016
# ---------------------------------------------------
Smote.exsRegress <- function(dat, tgt, N, k, dist, p)
  # INPUTS:
  # dat are the rare cases (the minority "class" cases)
  # tgt the column nr of the target variable
  # N is the percentage of over-sampling to carry out;
  # and k is the number of nearest neighours
  # dist is the distance function used for the neighours computation
  # p is an integer used when a "p-norm" distance is selected
  # OUTPUTS:
  # The result of the function is a (N-1)*nrow(dat) set of generate
  # examples with rare values on the target
{
  
  nomatr <- c()
  T <- matrix(nrow = dim(dat)[1], ncol = dim(dat)[2])
  for (col in seq.int(dim(T)[2])){
    if (class(dat[, col]) %in% c('factor', 'character')) {
      T[, col] <- as.integer(dat[, col])
      nomatr <- c(nomatr, col)
    } else {
      T[, col] <- dat[, col]
    }
  }
  nC <- dim(T)[2]
  nT <- dim(T)[1]
  

  ranges <- rep(1, nC)
  if (length(nomatr)) {
    for (x in (1:nC)[-c(nomatr)]) {
      ranges[x] <- max(T[, x]) - min(T[, x])
    }
  } else {
    for(x in (1:nC)) {
      ranges[x] <- max(T[, x]) - min(T[, x])
    }
  }

  kNNs <- neighbours(tgt, dat, dist, p, k)
    
  nexs <- as.integer(N - 1) # nr of examples to generate for each rare case
  extra <- as.integer(nT * (N - 1 - nexs)) # the extra examples to generate
  idx <- sample(1:nT, extra)
  newM <- matrix(nrow = nexs * nT + extra, ncol = nC)    # the new cases
 
  if (nexs) {
    for (i in 1:nT) {
      for (n in 1:nexs) {
        # select randomly one of the k NNs
        neig <- sample(1:k, 1)
        # the attribute values of the generated case
        difs <- T[kNNs[i, neig], -tgt] - T[i, -tgt]
        newM[(i - 1) * nexs + n, -tgt] <- T[i, -tgt] + runif(1) * difs
        for (a in nomatr) {
          # nominal attributes are randomly selected among the existing
          # values of seed and the selected neighbour 
          newM[(i - 1) * nexs + n, a] <- c(T[kNNs[i, neig], a],
                                          T[i, a])[1 + round(runif(1), 0)]
        }
        # now the target value (weighted (by inverse distance) average)
        d1 <- d2 <- 0
        for (x in (1:nC)[-c(nomatr, tgt)]) {
          d1 <- abs(T[i, x] - newM[(i - 1) * nexs + n, x])/ranges[x]
          d2 <- abs(T[kNNs[i, neig], x] - newM[(i - 1) * nexs + n, x])/ranges[x]
        }
        if (length(nomatr)) {
          d1 <- d1 + sum(T[i, nomatr] != newM[(i - 1) * nexs + n, nomatr])
          d2 <- d2 + 
               sum(T[kNNs[i, neig], nomatr] != newM[(i - 1) * nexs + n, nomatr])
        }
        # (d2+d1-d1 = d2 and d2+d1-d2 = d1) the more distant the less weight
        if (d1 == d2) {
          newM[(i - 1) * nexs + n, tgt] <- (T[i, tgt] + T[kNNs[i, neig], tgt])/2
        } else {
          newM[(i - 1) * nexs + n, tgt] <- (d2 * T[i, tgt] + 
                                           d1 * T[kNNs[i, neig], tgt])/(d1 + d2)
        }
      }
    }
  }
  
  if (extra) {
    count <- 1
    for (i in idx) {
      # select randomly one of the k NNs
      neig <- sample(1:k, 1) 
      
      # the attribute values of the generated case
      difs <- T[kNNs[i, neig], -tgt] - T[i, -tgt]
      newM[nexs * nT + count, -tgt] <- T[i, -tgt] + runif(1) * difs
      for (a in nomatr) {
        newM[nexs * nT + count, a] <- c(T[kNNs[i,neig], a], 
                                       T[i, a])[1 + round(runif(1), 0)]
      }
      
      # now the target value (weighted (by inverse distance) average)
      d1 <- d2 <- 0
      for (x in (1:nC)[-c(nomatr,tgt)]) {
        d1 <- abs(T[i, x] - newM[nexs * nT + count, x])/ranges[x]
        d2 <- abs(T[kNNs[i, neig], x] - newM[nexs * nT + count, x])/ranges[x]
      }
      if (length(nomatr)) {
        d1 <- d1 + sum(T[i,nomatr] != newM[(i - 1) * nexs + n, nomatr])
        d2 <- d2 + 
              sum(T[kNNs[i, neig], nomatr] != newM[(i - 1) * nexs + n, nomatr])
      }
      # (d2+d1-d1 = d2 and d2+d1-d2 = d1) the more distant the less weight
      if (d1 == d2) {
        newM[nexs * nT + count, tgt] <- (T[i, tgt] + T[kNNs[i, neig], tgt])/2
      } else {
        newM[nexs * nT + count, tgt] <- (d2 * T[i, tgt] + 
                                        d1 * T[kNNs[i, neig],tgt])/(d1 + d2)
      }
      
      count <- count + 1
    }
  }

  newCases <- data.frame(newM)

  for (a in nomatr) {
    newCases[, a] <- factor(newCases[, a],
                            levels = 1:nlevels(dat[, a]),
                            labels = levels(dat[, a]))
  }
  
  colnames(newCases) <- colnames(dat)
  newCases
  
}
