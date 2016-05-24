## ===================================================
## Creating a new training sample generated with the introduction
## of Gaussian Noise for regression problems
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae), ]
#   C.perc = list(0.5, 3) 
#   mygn.alg <- GaussNoiseRegress(a7~., clean.algae, C.perc = C.perc)
#   gnB.alg <- GaussNoiseRegress(a7~., clean.algae, C.perc = "balance", 
#                               pert = 0.1)
#   gnE.alg <- GaussNoiseRegress(a7~., clean.algae, C.perc = "extreme")
#   
#   plot(density(clean.algae$a7))
#   lines(density(gnE.alg$a7), col = 2)
#   lines(density(gnB.alg$a7), col = 3)
#   lines(density(mygn.alg$a7), col = 4)
# 
# 
# 
#   ir<- iris[-c(95:130), ]
#   mygn1.iris <- GaussNoiseRegress(Sepal.Width~., ir, C.perc = list(0.5, 2.5))
#   mygn2.iris <- GaussNoiseRegress(Sepal.Width~., ir, C.perc = list(0.2, 4),
#                                   thr.rel = 0.8)
#   gnB.iris <- GaussNoiseRegress(Sepal.Width~., ir, C.perc = "balance")
#   gnE.iris <- GaussNoiseRegress(Sepal.Width~., ir, C.perc = "extreme")
# 
#   rel <- matrix(0, ncol = 3, nrow = 0)
#   rel <- rbind(rel, c(2, 1, 0))
#   rel <- rbind(rel, c(3, 0, 0))
#   rel <- rbind(rel, c(4, 1, 0))
# 
#   gn.rel <- GaussNoiseRegress(Sepal.Width~., ir, rel = rel,
#                               C.perc = list(5, 0.2, 5))
# 
# plot(density(ir$Sepal.Width), ylim = c(0, 1))
# lines(density(gn.rel$Sepal.Width), col = 2)
# lines(density(gnB.iris$Sepal.Width), col = 3)
# lines(density(gnE.iris$Sepal.Width, bw = 0.3), col = 4)
# 
## P.Branco, May 2015 Apr 2016
## ---------------------------------------------------
GaussNoiseRegress <- function(form, dat, rel = "auto", thr.rel = 0.5,
                              C.perc = "balance", pert = 0.1, repl = FALSE)
  
  # Args:
  # form    a model formula
  # dat    the original training set (with the unbalanced distribution)
  # C.perc  named list containing each class percentage of under- or 
  #         over-sampling to apply between 0 and 1. The user may provide
  #         only a subset of the existing classes where sampling is to
  #         be applied. Alternatively it may be "balance" or "extreme",
  #         cases where the sampling percentages are automatically estimated.
  # pert    the level of perturbation to introduce when generating synthetic
  #         examples
  # repl    is it allowed to perform sampling with replacement (when 
  #         performing under-sampling)
  #
  # Returns: a data frame modified by the Gaussian Noise strategy

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
    stop("All the points have relevance 1. Please, redefine your 
         relevance function!")
  }
  if (!length(which(temp > 0))) {
    stop("All the points have relevance 0. Please, redefine your
         relevance function!")
  }

  temp[which(y.relev > thr.rel)] <- -temp[which(y.relev > thr.rel)]
  bumps <- c()
  for (i in 1:(length(y) - 1)) {
    if (temp[i] * temp[i + 1] < 0) {
      bumps <- c(bumps, i)
    }
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
    if (length(C.perc) != nbump) {
      stop("The percentages provided must be the same length as the number
           of bumps!")
    }
  } else if (C.perc == "balance") { 
    # estimate the percentages of over-/under-sampling
    B <- round(nrow(dat)/nbump, 0)
    C.perc <- B/sapply(obs.ind, length)        
  } else if (C.perc == "extreme") {
    B <- round(nrow(dat)/nbump, 0)
    rescale <- nbump * B/sum(B^2/sapply(obs.ind, length))
    obj <- round((B^2/sapply(obs.ind, length)) * rescale, 2)
    C.perc <- round(obj/sapply(obs.ind, length), 4)
  }
  
  for (i in 1:nbump) {
    if (C.perc[[i]] == 1) {
      newdata <- rbind(newdata, dat[names(obs.ind[[i]]), ])
    } else if (C.perc[[i]] > 1) {
      newExs <- Gn.exsRegress(dat[names(obs.ind[[i]]), ],
                              ncol(dat),
                              C.perc[[i]],
                              pert)
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
# Obtain a set of synthetic examples generated with Gaussian Noise 
# perturbance for a set of rare cases.
#
# 
# P.Branco, May 2015 Apr 2016
# ---------------------------------------------------
Gn.exsRegress <- function(dat, tgt, N, pert)
  # INPUTS:
  # dat are the rare cases (the minority "class" cases)
  # tgt the column nr of the target variable
  # N is the percentage of over-sampling to carry out;
  # pert is the amount of disturbance between 0 and 1 of standard deviation
  # OUTPUTS:
  # The result of the function is a (N-1)*nrow(dat) set of generate
  # examples with rare values on the target
{
  nC <- dim(dat)[2]
  nL <- dim(dat)[1]
  nomatr <- c()
  T <- matrix(nrow = nL, ncol = nC)
  for (col in seq.int(nC)) {
    if (class(dat[, col]) %in% c('factor', 'character')) {
      T[, col] <- as.integer(dat[, col])
      nomatr <- c(nomatr, col)
    }else {
      T[, col] <- dat[, col]
    }
    if (length(nomatr)) {
      numatr <- (1:nC)[-nomatr]    
    } else {
      numatr <- (1:nC)
    }
  }

  nexs <-  as.integer(N - 1) # number of artificial exs to generate for each rare case
  extra <- as.integer(nL * (N - 1 - nexs)) # the extra examples to generate
  id.ex <- sample(1:nL, extra)
  
  
  newdata <- matrix(nrow = nexs * nL + extra, ncol = nC)
  
  if (nexs) {
    for (i in 1:nL) {
      for (n in 1:nexs) {
        # the attribute values of the generated case
        idx <- (i - 1) * nexs + n 
        for (num in 1:nC) {
          if (is.na(T[i, num])) {
            newdata[idx, num] <- NA
          } else {
            newdata[idx, num] <- T[i, num] + rnorm(1,
                                                   0,
                                                   sd(T[, num], 
                                                      na.rm = TRUE) * pert)
            if (num %in% nomatr) {
              probs <- c()
              if (length(unique(T[, num])) == 1) {
                newdata[idx, num] <- T[1, num]
              } else {
                for (u in 1:length(unique(T[, num]))) {
                  probs <- c(probs, 
                             length(which(T[, num] == unique(T[, num])[u])))
                }
                newdata[idx, num] <- sample(unique(T[, num]), 1, prob = probs)
              }
            }
          }
        }
      }
    }
  }
  
  
  if (extra) {
    count <- 1
    for (i in id.ex) {    
      for (num in 1:nC) {
        if (is.na(T[i, num])) {
          newdata[nexs * nL + count, num] <- NA
        } else {
          newdata[nexs * nL + count, num] <- T[i, num] + 
                                             rnorm(1, 0, 
                                                   sd(T[, num],
                                                      na.rm = TRUE) * pert)
          if (num %in% nomatr) {
            probs <- c()
            if (length(unique(T[, num])) == 1) {
              newdata[nexs * nL + count, num] <- T[1, num]
            } else {
              for (u in 1:length(unique(T[, num]))) {
                probs <- c(probs,
                           length(which(T[, num] == unique(T[, num])[u])))
              }
              newdata[nexs * nL + count, num] <- sample(unique(T[, num]),
                                                        1, prob = probs)
            }
          }
        }
      }
      count <- count + 1
    }
  }
  
  newCases <- data.frame(newdata)
  
  for(a in nomatr){
    newCases[, a] <- factor(newCases[, a],
                            levels = 1:nlevels(dat[, a]),
                            labels = levels(dat[, a]))
  }
  
  colnames(newCases) <- colnames(dat)
  newCases
    
}
