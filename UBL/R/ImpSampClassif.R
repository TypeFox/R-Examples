## ===================================================
## Creating a new training sample for classification problems
# based on the relevance assigned to each class to perform 
# both random over and under-sampling
# 
# Examples:
# ir <- iris[-c(51:70,111:150), ]
# IS.ext <- ImpSampClassif(Species~., ir, C.perc = "extreme")
# IS.bal <- ImpSampClassif(Species~., ir, C.perc = "balance")
# myIS <- ImpSampClassif(Species~., ir, 
#                        C.perc = list(setosa = 0.2,
#                                      versicolor = 2,
#                                      virginica = 6))
# P. Branco, July 2015 Apr 2016
# ---------------------------------------------------
ImpSampClassif <- function(form, dat, C.perc = "balance")
  
  # Args:
  # form   a model formula
  # dat    the original training set (with the unbalanced distribution)
  # C.perc is a list containing the percentage of under- or/and 
  #        over-sampling to apply to each "class" obtained with the 
  #        threshold. To use this list, a thr.rel must be provided otherwise
  #        this parameter is ignored.
  #        The over-sampling percentage means that the examples above the 
  #        threshold are increased by this percentage. The under-sampling
  #        percentage means that the normal cases (cases below the threshold)
  #        are under-sampled by this percentage. Alternatively it may be
  #        "balance" (the default) or "extreme", cases where the sampling 
  #        percentages are automatically estimated.
  #
  # Returns: a data frame with the data modified through the Importance 
  #        Sampling strategy.

{
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  names <- sort(unique(dat[, tgt]))
  li <- class.freq(dat, tgt)
  

  if (is.list(C.perc)) {
    
    names.und <- names(which(C.perc < 1))
    names.ove <- names(which(C.perc > 1))
    names.same <- setdiff(names, union(names.und, names.ove))
    
    # include examples from classes unchanged
    newdata <- dat[which(dat[, tgt] %in% names.same), ]
    
    if (length(names.und)) {  # perform under-sampling
      for (i in 1:length(names.und)) { 
        Exs <- which(dat[, tgt] == names.und[i])
        sel <- sample(Exs,
                      as.integer(C.perc[[names.und[i]]] * length(Exs)),
                      replace = FALSE)
        newdata <- rbind(newdata, dat[sel, ])
      }
    }
    
    if (length(names.ove)) {  # perform over-sampling
      for (i in 1:length(names.ove)) { 
        Exs <- which(dat[, tgt] == names.ove[i])
        if (length(Exs) == 1) {
          sel <- rep(Exs, C.perc[[names.ove[i]]])
        } else {
          sel <- sample(Exs,
                        as.integer(C.perc[[names.ove[i]]] * length(Exs)),
                        replace = TRUE)
        }
        newdata <- rbind(newdata, dat[sel, ])
      }
    }

  } else {
    
    if (C.perc == "balance") {  
      li[[3]] <- round(sum(li[[2]])/length(li[[2]]), 0) - li[[2]]
    } else if (C.perc == "extreme") {
      med <- sum(li[[2]])/length(li[[2]])
      li[[3]] <- round(med^2/li[[2]] * sum(li[[2]])/sum(med^2/li[[2]]), 
                       0) - li[[2]]
    } else {
      stop("Please provide a list with classes to under/over-sample 
           or alternatively provide 'balance' or 'extreme'.")
    }
    
    und <-which(li[[3]] < 0) # classes to under-sample
    ove <- which(li[[3]] > 0) #classes to over-sample
    same <- which(li[[3]] == 0) # unchanged classes
    
    # include examples from classes unchanged
    newdata <- dat[which(dat[, tgt] %in% li[[1]][same]), ]
    
    if (length(und)) { #perform under-sampling
      for (i in 1:length(und)) { 
        Exs <- which(dat[, tgt] == li[[1]][und[i]])
        sel <- sample(Exs,
                      as.integer(li[[2]][und[i]] + li[[3]][und[i]]),
                      replace = FALSE)
        newdata <- rbind(newdata, dat[sel, ])
      }
    }
    
    if (length(ove)) { #perform over-sampling
      for (i in 1:length(ove)) {
        Exs <- which(dat[, tgt] == li[[1]][ove[i]])
        if (length(Exs) == 1) {
          sel <- rep(Exs, as.integer(li[[2]][ove[i]] + li[[3]][ove[i]]))
        } else {
          sel <- sample(Exs,
                        as.integer(li[[2]][ove[i]] + li[[3]][ove[i]]),
                        replace = TRUE)
        }
        newdata <- rbind(newdata, dat[sel, ])
      } 
    }
    
  }

  newdata
}
