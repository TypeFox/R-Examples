## ===================================================
## Creating a new training sample generated with the introduction
## of Gaussian Noise for classification problems
## Examples:
## library(DMwR)
## data(algae)
## clean.algae <- algae[complete.cases(algae), ]
## C.perc = list(autumn = 2, summer = 1.5, winter = 0.9)
## gn1 <- GaussNoiseClassif(season~., clean.algae)
## gn2 <- GaussNoiseClassif(season~., clean.algae, C.perc)
## P.Branco, May 2015 April 2016
## ---------------------------------------------------
GaussNoiseClassif <- function(form, dat, C.perc = "balance", pert = 0.1,
                              repl = FALSE)
  # Args:
  # form    a model formula
  # dat     the original training set (with the unbalanced distribution)
  # C.perc  named list containing each class percentage of under- or 
  #         over-sampling to apply. The user may provide
  #         only a subset of the existing classes where sampling is to
  #         be applied. Alternatively it may be "balance" or "extreme",
  #         cases where the sampling percentages are automatically estimated.
  # pert    the level of perturbation to introduce when generating synthetic 
  #         examples. Assuming as center the base example, this parameter 
  #         defines the radius (based on the standard deviation) where the
  #         new example is generated. 
  # repl    is it allowed to perform sampling with replacement (when 
  #         performing under-sampling)
  #
  # Returns: a data frame with the new modified data set through the 
  #         Gaussian Noise strategy
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
                      replace=repl)
        newdata <- rbind(newdata,dat[sel, ])
      }
    }
    
    if (length(names.ove)) { # perform over-sampling
      for (i in 1:length(names.ove)) {
        newExs <- Gn.exsClassif(dat[which(dat[, ncol(dat)] == names.ove[i]), ],
                                ncol(dat),
                                C.perc[[names.ove[i]]],
                                pert)
        # add original rare examples and synthetic generated examples
        newdata <- rbind(newdata, newExs,
                         dat[which(dat[, ncol(dat)] == names.ove[i]), ])
      }
    }
  } else {
    if (C.perc == "balance") {
      li[[3]] <- round(sum(li[[2]]) / length(li[[2]]), 0) - li[[2]]
    } else if (C.perc =="extreme") {
      med <- sum(li[[2]])/length(li[[2]])
      li[[3]] <- round(med^2/li[[2]] * sum(li[[2]])/sum(med^2/li[[2]]), 0) - 
                 li[[2]]
    } else {
      stop("Please provide a list with classes to under/over-sample or 
           'balance' or 'extreme'.")
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
        newExs <- Gn.exsClassif(dat[which(dat[, ncol(dat)] == li[[1]][ove[i]]),],
                                ncol(dat),
                                li[[3]][ove[i]]/li[[2]][ove[i]] + 1,
                                pert)
        # add original rare examples and synthetic generated examples
        newdata <- rbind(newdata, newExs, 
                         dat[which(dat[, ncol(dat)] == li[[1]][ove[i]]), ])
      } 
    }
  }
  
  if (tgt < ncol(dat)) {
    newdata <- newdata[, cols]
    dat <- dat[, cols]
  }
  
  newdata
}


# ===================================================
# Obtain a set of synthetic examples generated with Gaussian Noise 
# perturbance for a set of rare cases.
# 
# P.Branco, May 2015
# ---------------------------------------------------
Gn.exsClassif <- function(dat, tgt, N, pert)
  # Args:
  # dat are the minority class cases
  # tgt the column nr of the target variable
  # N is the percentage of over-sampling to carry out;
  # pert is the amount of disturbance between 0 and 1 of standard deviation 
  # Returns:
  # The result of the function is a (N-1)*nrow(dat) set of generated
  # examples with rare class on the target
{
  nC <- dim(dat)[2]
  nL <- dim(dat)[1]
  
  nomatr <- c()
  T <- matrix(nrow = nL,ncol = nC - 1)
  for(col in seq.int(nC - 1))
    if (class(dat[, col]) %in% c('factor','character')) {
      T[, col] <- as.integer(dat[, col])
      nomatr <- c(nomatr, col)
    } else {
      T[, col] <- dat[, col]
    }
  
  numatr <- (1:nC)[-c(nomatr, tgt)]
    
  # number of artificial exs to generate for each rare case
  nexs <- as.integer(N - 1) 
  # the extra examples to generate
  extra <- as.integer(nL * (N - 1 - nexs))
  id.ex <- sample(1:nL, extra)
  

  newdata <- matrix(nrow = nexs * nL + extra, ncol = nC)

  if (nexs) {
    for (i in 1:nL) {
      for (n in 1:nexs) {
        idx <- (i - 1) * nexs + n 
        for (num in 1:(nC - 1)) {
          newdata[idx, num] <- T[i, num] + rnorm(1, 0, sd(T[, num]) * pert)
          if (num %in% nomatr) {
            probs <- c()
            for (u in 1:length(unique(T[, num]))) {
              probs <- c(probs, length(which(T[, num] == unique(T[, num])[u])))
            }
            newdata[idx, num] <- sample(unique(T[, num]), 1, prob = probs)
          }
        }
      }
    }
  }
  
  if (extra) {
    count <- 1
    for (i in id.ex) {
      for (num in 1:(nC-1)) {
        newdata[nexs * nL + count, num] <- T[i, num] + 
                                           rnorm(1, 0, sd(T[, num]) * pert)
        if (num %in% nomatr) {
          probs <- c()
          for (u in 1:length(unique(T[, num]))) {
            probs <- c(probs,length(which(T[, num] == unique(T[, num])[u])))
          }
          newdata[nexs * nL + count, num] <- sample(unique(T[, num]),
                                                    1, prob = probs)
        }
      }
      count <- count + 1
    }
  }
  
  newCases <- data.frame(newdata)
  
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

