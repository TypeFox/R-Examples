## ===================================================
## Creating a SMOTEd training sample for classification problems
# 
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae), ]
#   C.perc = list(autumn = 2, summer = 1.5, winter = 0.9) 
#   class spring remains unchanged
#   mysmote.algae <- SmoteClassif(season~., clean.algae, C.perc)
#   smoteBalan.algae <- SmoteClassif(season~., clean.algae, "balance")
#   smoteExtre.algae <- SmoteClassif(season~., clean.algae, "extreme")
# 
# data(iris)
# dat <- iris[, c(1, 2, 5)]
# dat$Species <- factor(ifelse(dat$Species == "setosa", "rare", "common")) 
# ## checking the class distribution of this artificial data set
# table(dat$Species)
# 
# ## now using SMOTE to create a more "balanced problem"
# newData <- SmoteClassif(Species ~ ., dat,
#                         C.perc = list(common = 1,rare = 6))
# table(newData$Species)
# 
# ## Checking visually the created data
#   par(mfrow = c(1, 2))
#   plot(dat[, 1], dat[, 2], pch = 19 + as.integer(dat[, 3]),
#        main = "Original Data")
#   plot(newData[, 1], newData[, 2], pch = 19 + as.integer(newData[, 3]),
#        main = "SMOTE'd Data")
# ## Another example with iris data 
#   ir <- iris[-c(95:130), ]
#   mysmote.iris <- SmoteClassif(Species~., ir, 
#                                list(setosa = 0.6, virginica = 1.5))
#   smoteBalan.iris <- SmoteClassif(Species~., ir, "balance")
#   smoteExtre.iris <- SmoteClassif(Species~., ir, "extreme")
# 
# 
#   library(MASS)
#   data(cats)
#   mysmote.cats <- SmoteClassif(Sex~., cats, list(M = 0.8, F = 1.8))
#   smoteBalan.cats <- SmoteClassif(Sex~., cats, "balance")
#   smoteExtre.cats <- SmoteClassif(Sex~., cats, "extreme")
#
## L. Torgo, Feb 2010, Nov 2014
## P.Branco, Mar,Apr 2015 Apr 2016
## ---------------------------------------------------
SmoteClassif <- function(form, dat, C.perc = "balance",
                         k = 5, repl = FALSE, dist = "Euclidean", p = 2)
  # Args:
  # form    a model formula
  # dat    the original training set (with the unbalanced distribution)
  # C.perc  named list containing each class percentage of under- or 
  #         over-sampling to apply between 0 and 1. The user may provide
  #         only a subset of the existing classes where sampling is to
  #         be applied. Alternatively it may be "balance" (the default) or
  #         "extreme", cases where the sampling percentages are automatically
  #         estimated.
  # k       is the number of neighbors to consider as the pool from where
  #         the new examples are generated
  # repl    is it allowed to perform sampling with replacement, when 
  #         performing under-sampling
  # dist    is the distance measure to be used (defaults to "Euclidean")
  # p       is a parameter used when a p-norm is computed
  #
  # Returns: a new data frame modified through the smote algorithm

{
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  names <- sort(unique(dat[, tgt]))
  li <- class.freq(dat, tgt)
  if (tgt < ncol(dat)) {
    orig.order <- colnames(dat)
    cols <- 1:ncol(dat)
    cols[c(tgt, ncol(dat))] <- cols[c(ncol(dat), tgt)]
    dat <- dat[, cols]
  }
  
  if (is.list(C.perc)) {
    names.und <- names(which(C.perc < 1))
    names.ove <- names(which(C.perc > 1))
    names.same <- setdiff(names, union(names.und, names.ove))
    
    # include examples from classes unchanged
    newdata <- dat[which(dat[, ncol(dat)] %in% names.same), ]

    if (length(names.und)) {  # perform under-sampling
      for (i in 1:length(names.und)) {
        Exs <- which(dat[, ncol(dat)] == names.und[i])
        sel <- sample(Exs,
                      as.integer(C.perc[[names.und[i]]] * length(Exs)),
                      replace = repl)
        newdata <- rbind(newdata, dat[sel, ])
      }
    }
    if (length(names.ove)) { # perform over-sampling
      for (i in 1:length(names.ove)) {
        newExs <- Smote.exsClassif(dat[which(dat[, ncol(dat)] == names.ove[i]), ],
                                   ncol(dat),
                                   C.perc[[names.ove[i]]],
                                   k,
                                   dist,
                                   p)
        # add original rare examples and synthetic generated examples
        newdata <- rbind(newdata,
                         newExs,
                         dat[which(dat[, ncol(dat)] == names.ove[i]), ])
      }
    }
  } else {
    if (C.perc == "balance") {  
      li[[3]] <- round(sum(li[[2]])/length(li[[2]]), 0) - li[[2]]
    } else if (C.perc == "extreme") {
      med <- sum(li[[2]])/length(li[[2]])
      li[[3]] <- round(med^2/li[[2]] * sum(li[[2]])/sum(med^2/li[[2]]), 0) - li[[2]]
    } else {
      stop("Please provide a list with classes to under-/over-sample
           or alternatively indicate 'balance' or 'extreme'.")
    }
    und <- which(li[[3]] < 0) # classes to under-sample
    ove <- which(li[[3]] > 0) #classes to over-sample
    same <- which(li[[3]] == 0) # unchanged classes
    
    # include examples from classes unchanged
    newdata <- dat[which(dat[, ncol(dat)] %in% li[[1]][same]), ]
    
    if (length(und)) { #perform under-sampling
      for (i in 1:length(und)) { 
        Exs <- which(dat[, ncol(dat)] == li[[1]][und[i]])
        sel <- sample(Exs,
                      as.integer(li[[2]][und[i]] + li[[3]][und[i]]),
                      replace = repl)
        newdata <- rbind(newdata, dat[sel, ])
      }
    }
    
    if (length(ove)) { #perform over-sampling
      for (i in 1:length(ove)) {
        newExs <- Smote.exsClassif(dat[which(dat[, ncol(dat)] == li[[1]][ove[i]]), ],
                                   ncol(dat),
                                   li[[3]][ove[i]]/li[[2]][ove[i]] + 1,
                                   k,
                                   dist,
                                   p)
        # add original rare examples and synthetic generated examples
        dc <- ncol(dat)
        newdata <- rbind(newdata, newExs, dat[which(dat[,dc] == li[[1]][ove[i]]),])
      } 
    }

  }
  
  if (tgt < ncol(dat)) {
    newdata <- newdata[,cols]
    dat <- dat[,cols]
  }
  
  newdata
}


# ===================================================
# Obtain a set of smoted examples for a set of rare cases.
# L. Torgo, Feb 2010
# P.Branco, Mar,Apr 2015
# ---------------------------------------------------
Smote.exsClassif <- function(dat, tgt, N, k, dist, p)
  # INPUTS:
  # dat   are the rare cases (the minority class cases)
  # tgt    is the name of the target variable
  # N      is the percentage of over-sampling to carry out;
  # k      is the number of nearest neighors to use for the generation
  # dist   is the distance function to use for the neighbors computation
  # p      is an integer used when a "p-norm" distance is selected
  # OUTPUTS:
  # The result of the function is a (N-1)*nrow(dat) set of generated
  # examples with rare class on the target
{
  nomatr <- c()
  T <- matrix(nrow = dim(dat)[1], ncol = dim(dat)[2] - 1)
  for (col in seq.int(dim(T)[2])) { 
    if (class(dat[, col]) %in% c('factor', 'character')) {
      T[, col] <- as.integer(dat[, col])
      nomatr <- c(nomatr, col)
    } else {
      T[, col] <- dat[, col]
    }
  }
  nC <- dim(T)[2]
  nT <- dim(T)[1]
  
  # check if there is enough data to determine the k neighbors
  if (nT <= k) {
    stop("Trying to determine ",k, " neighbors for a subset with only ",
         nT, " examples")
  }

  kNNs <- neighbours(tgt, dat, dist, p, k)
  
  nexs <-  as.integer(N - 1) # nr of examples to generate for each rare case
  extra <- as.integer(nT * (N - 1 - nexs)) # the extra examples to generate
  idx <- sample(1:nT, extra)
  newM <- matrix(nrow = nexs * nT + extra, ncol = nC)    # the new cases
  if (nexs) {
    for (i in 1:nT) {
      for (n in 1:nexs) {
        # select randomly one of the k NNs
        neig <- sample(1:k, 1)
      
        # the attribute values of the generated case
        difs <- T[kNNs[i, neig], ] - T[i, ]
        newM[(i - 1) * nexs + n, ] <- T[i, ] + runif(1) * difs
        for (a in nomatr) {
          # nominal attributes are randomly selected among the existing values
          # of seed and the selected neighbor 
          newM[(i - 1) * nexs + n, a] <- c(T[kNNs[i, neig], a], 
                                          T[i, a])[1 + round(runif(1), 0)]
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
      difs <- T[kNNs[i, neig], ] - T[i, ]
      newM[nexs * nT + count, ] <- T[i, ] + runif(1) * difs
      for (a in nomatr) {
        newM[nexs * nT + count, a] <- c(T[kNNs[i, neig], a], 
                                       T[i, a])[1 + round(runif(1), 0)]
      }
      count <- count + 1
    }
  }
  newCases <- data.frame(newM)
  
  for (a in nomatr){
    newCases[, a] <- factor(newCases[, a],
                            levels = 1:nlevels(dat[, a]),
                            labels = levels(dat[, a]))
  }
  newCases[, tgt] <- factor(rep(dat[1, tgt], nrow(newCases)),
                            levels = levels(dat[, tgt]))
  colnames(newCases) <- colnames(dat)
  newCases
}
