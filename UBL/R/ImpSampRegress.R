## ===================================================
## Creating a new training sample for regression problems
# based on the relevance function to perform both over and under-sampling
# 
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae), ]
#   IS.ext <-ImpSampRegress(a7~., clean.algae, rel = "auto", thr.rel = 0.7,
#                           C.perc = "extreme")
#   IS.bal <-ImpSampRegress(a7~., clean.algae, rel = "auto", thr.rel = 0.7,
#                           C.perc = "balance")
#   myIS <-ImpSampRegress(a7~., clean.algae, rel = "auto", thr.rel = 0.7,
#                           C.perc = list(0.2, 6))
##  everything automatic
#   IS.auto <- ImpSampRegress(a7~., clean.algae)
##  select the importance given to phi (over-sampling) 
##  and to 1-phi (under-sampling)
#   IS.auto2 <- ImpSampRegress(a7~., clean.algae, O = 0.8, U = 0.2)
# 
# P. Branco, May 2015 Apr 2016
# ---------------------------------------------------
ImpSampRegress <- function(form, dat, rel = "auto", thr.rel = NA, 
                           C.perc = "balance", O = 0.5, U = 0.5)
  
  # Args:
  # form    a model formula
  # dat    the original training set (with the unbalanced distribution)
  # rel     is the relevance determined automatically (default: "auto") 
  #         or provided by the user through a matrix.
  # thr.rel is the relevance threshold above which a case is considered
  #         as belonging to the rare "class". If thr.rel is NA then both 
  #         under- and over-sampling are applied sampling examples
  #         accordingly to the relevance of the target variable (over-sampling)
  #         or accordingly to 1-phi of the target (uder-sampling)
  #       
  # C.perc is a list containing the percentage of under- or/and 
  #        over-sampling to apply to each "class" obtained with the threshold.
  #        To use this list, a thr.rel must be provided otherwise this 
  #        parameter is ignored. The over-sampling percentage means that the
  #        examples above the threshold are increased by this percentage. 
  #        The under-sampling percentage means that the normal cases (cases
  #        below the threshold) are under-sampled by this percentage.
  #        Alternatively it may be "balance" or "extreme", cases where the
  #        sampling percentages are automatically estimated.
  # O      is a number expressing the importance given to over-sampling when the
  #        thr.rel parameter is NA. When O increases the number of examples to
  #        include during the over-sampling step also increase. Defaults to 0.5.
  # U      is a number expressing the importance given to under-sampling when 
  #        the thr.rel parameter is NA. When U increases, the number of
  #        examples removed during the under-sampling step also increase.
  #        The default is 0.5.
  #
  # Returns: a new data frame modified through the Importance Sampling strategy

{
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  
  y <- dat[, tgt]
  attr(y, "names") <- rownames(dat)

  if (is.matrix(rel)) {
    pc <- phi.control(y, method = "range", control.pts = rel)
  } else if (is.list(rel)) { 
    pc <- rel
  } else if (rel == "auto") {
    pc <- phi.control(y, method = "extremes")
  } else {# handle other relevance functions and not using the threshold!
    stop("future work!")
  }

  if (!is.na(thr.rel)) {
    s.y <- sort(y)
    temp <- y.relev <- phi(s.y, pc)
    if(!length(which(temp < 1))) {
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
      # estimate the percentages of over/under sampling
      B <- round(nrow(dat)/nbump, 0)
      C.perc <- B/sapply(obs.ind, length)        
    } else if (C.perc == "extreme") {
      B <- round(nrow(dat)/nbump,0)
      rescale <- nbump * B/sum(B^2/sapply(obs.ind, length))
      obj <- round((B^2/sapply(obs.ind, length)) * rescale, 2)
      C.perc <- round(obj/sapply(obs.ind, length), 1)
    }
    
    for (i in 1:nbump) {
      if (C.perc[[i]] == 1) {
        newdata <- rbind(newdata, dat[names(obs.ind[[i]]), ])
      } else if (C.perc[[i]] > 1) {
        s <- sample(names(obs.ind[[i]]), 
                    as.integer(C.perc[[i]] * length(obs.ind[[i]])), 
                    replace = TRUE, 
                    prob = y.relev[which(s.y %in% obs.ind[[i]])])
        newdata <- rbind(newdata, dat[s, ])
        
      }else if (C.perc[[i]] < 1) {
        s <- sample(names(obs.ind[[i]]), 
                    as.integer(C.perc[[i]] * length(obs.ind[[i]])),
                    replace = TRUE, 
                    prob = y.relev[which(s.y %in% obs.ind[[i]])])
        newdata <- rbind(newdata, dat[s, ])
        
      }
    }
  }
  
  if (is.na(thr.rel)) {
    y.relev <- phi(y, pc)
    if (!length(which(y.relev < 1))) {
      stop("All the points have relevance 1. 
           Please, redefine your relevance function!")
    }
    zero <- which(y.relev == 0)
    if (length(zero)) {
      s.ove <- sample(setdiff(1:nrow(dat), zero),
                      as.integer(O * nrow(dat)),
                      replace = TRUE, 
                      prob = y.relev[-zero])
    } else {
      s.ove <- sample(1:nrow(dat),
                      as.integer(O * nrow(dat)), 
                      replace = TRUE,
                      prob = y.relev)
    }
    newdata <- rbind(dat, dat[s.ove, ])
    one <- which(y.relev == 1)
    s.und <- sample(setdiff(1:nrow(dat), one), 
                    as.integer(U * nrow(dat)), 
                    replace = TRUE,
                    prob = 1 - y.relev[-one])
    newdata <- newdata[-s.und,]
  }

  newdata
}
