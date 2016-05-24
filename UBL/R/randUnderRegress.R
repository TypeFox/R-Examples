# ===================================================
# Performs a random under-sampling strategy for regression problems.
# Basically randomly removes a percentage of cases of the "class(es)"
# (bumps below a relevance threshold) selected by the user. 
# Alternatively, it can either balance all the 
# existing classes or it can "smoothly invert" the frequency
# of the examples in each "class".
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae), ]
#   C.perc = list(0.5) 
#   alg.myUnd <- RandUnderRegress(a7~., clean.algae, C.perc = C.perc)
#   alg.Bal <- RandUnderRegress(a7~., clean.algae, C.perc = "balance")
#   alg.Ext <- RandUnderRegress(a7~., clean.algae, C.perc = "extreme")
# 
# L. Torgo, Jun 2008
# P. Branco, May 2015 Apr 2016
# ---------------------------------------------------
RandUnderRegress <- function(form, dat, rel = "auto", thr.rel = 0.5,
                             C.perc = "balance", repl = FALSE)
  # Args:
  # form a model formula
  # dat the original training set (with the unbalanced distribution)
  # rel is relevance determined automatically (default) 
  #     or provided by the user
  # thr.rel is the relevance threshold above which a case is considered
  #         as belonging to the rare "class"
  # C.perc is a list containing the under-sampling percentage/s to apply 
  #        to all/each "class" obtained with the relevance threshold. This
  #        percentage represents the percentage of examples that is maintained
  #        in each "class". Examples are randomly removed in each "class".
  #        Moreover, different percentages may be provided for each "class". 
  #        Alternatively, it may be "balance" or "extreme", cases 
  #        where the under-sampling percentages are automatically estimated.
  # repl is it allowed to perform sampling with replacement.
  #
  # Returns: a data frame modified by the random under-sampling strategy

{
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  
  y <- dat[, tgt]
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


  imp <- sapply(obs.ind, function(x) mean(phi(x, pc)))
  
  und <- which(imp < thr.rel)
  ove <- which(imp > thr.rel)
  
  newdata <- NULL
  for (j in 1:length(ove)) {
  newdata <- rbind(newdata, dat[names(obs.ind[[ove[j]]]), ])
  # start with the examples from the minority "classes"
  }
  
  # set the under-sampling percentages
  if (is.list(C.perc)) {
    if (length(und) > 1 & length(C.perc) == 1) { 
      # the same under-sampling percentage is applied to all the "classes"
      C.perc <- rep(C.perc[1], length(und))
    } else if (length(und) > length(C.perc) & length(C.perc) > 1) {
      stop("The number of under-sampling percentages must be equal 
           to the number of bumps below the threshold defined!")      
    }else if (length(und) < length(C.perc)) {
      stop("the number of under-sampling percentages must be at most
           the number of bumps below the threshold defined!")
    }
  } else if (C.perc == "balance") {
    B <- sum(sapply(obs.ind[ove], length))
    obj <- B/length(und)
    C.perc <- as.list(round(obj/sapply(obs.ind[und], length), 5))
  } else if (C.perc == "extreme") {
    Bove <- sum(sapply(obs.ind[ove], length))/length(ove)
    obj <- Bove^2/sapply(obs.ind[und], length)
    C.perc <- as.list(round(obj/sapply(obs.ind[und], length), 5))
  }
  
  for (j in 1:length(und)) {
    sel <- sample(names(obs.ind[[und[j]]]),
                  C.perc[[j]] * length(obs.ind[[und[j]]]),
                  replace = repl)
    newdata <- rbind(newdata, dat[sel, ])
  }
  
  newdata
}

