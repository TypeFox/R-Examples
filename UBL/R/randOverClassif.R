# ===================================================
# Random over-sampling strategy for multiclass problems.
# Basically randomly copies a percentage of cases of the class(es) 
# selected by the user. Alternatively, it can either balance all the 
# existing classes or it can "smoothly invert" the frequency
# of the examples in each class
# Examples:
#   ir <- iris[-c(95:130), ]
#   myover.iris <- RandOverClassif(Species~., ir, list(versicolor = 1.2, 
#                                                     virginica = 2.3))
#   oveBalan.iris <- RandOverClassif(Species~., ir, "balance")
#   oveInvert.iris <- RandOverClassif(Species~., ir, "extreme")
# 
#   library(DMwR)
#   data(algae)
#   classes spring and winter remain unchanged:
#   C.perc=list(autumn = 2, summer = 1.5, spring = 1) 
#   myover.algae <- RandOverClassif(season~., algae, C.perc)
#   oveBalan.algae <- RandOverClassif(season~., algae, "balance")
#   oveInvert.algae <- RandOverClassif(season~., algae, "extreme")
#   
#   library(MASS)
#   data(cats)
#   myover.cats <- RandOverClassif(Sex~., cats, list(M = 1.5))
#   oveBalan.cats <- RandOverClassif(Sex~., cats, "balance")
#   oveInvert.cats <- RandOverClassif(Sex~., cats, "extreme")
#
# P. Branco, Mar 2015
# ---------------------------------------------------
RandOverClassif <- function(form, dat, C.perc = "balance", repl = TRUE)
  # INPUTS:
  # form a model formula
  # dat the original training set (with the unbalanced distribution)
  # C.perc is a named list containing each class over-sampling percentage(>=1).
  #       The user may only provide the classes where he wants to 
  #       apply the random over-sampling strategy.
  #       Alternatively it may be "balance" or "extreme", cases where 
  #       the over-sampling percentages are automatically estimated
  # repl is it allowed or not to perform sampling with replacement
  #       defaults to TRUE because if the over-sampling percentage is 
  #       >2 this is necessary.
{
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  names <- sort(unique(dat[, tgt]))
  li <-class.freq(dat, tgt)
  
  # include base examples (i.e., the original data set)
  newdata <- dat
  
  if (is.list(C.perc)) { # over-sampling percentages are provided by the user
    if (any(C.perc < 1)) {
      stop("percentages provided must be > 1 to perform over-sampling!")
    }
    names.ove <- names(which(C.perc > 1))
    
    for (i in 1:length(names.ove)) { # over-sampling each class provided
      Exs <- which(dat[, tgt] == names.ove[i])
      sel <- sample(Exs,
                    as.integer((C.perc[[names.ove[i]]] - 1) * length(Exs)),
                    replace = repl)
      newdata <- rbind(newdata, dat[sel, ])
    }
  } else if (C.perc == "balance") { # over-sampling percent. will be calculated
    majCl <- names(which(table(dat[, tgt]) == max(table(dat[, tgt]))))
    
    if (length(majCl) == length(names)) {
      stop("Classes are already balanced!")
    }
    
    names.ove <- names[which(!(names %in% majCl))]
    
    # over-sample all the other classes
    for (i in 1:length(names.ove)) {
      Exs <- which(dat[, tgt] == names.ove[i])
      num1 <- li[[2]][as.numeric(match(majCl, names))[1]]
      num2<-  li[[2]][as.numeric(names.ove[i])]
      sel <- sample(Exs,
                    as.integer(num1 - num2),
                    replace = repl)
      newdata <- rbind(newdata, dat[sel, ])
    }
  } else if (C.perc == "extreme") { 
    # "reverse" the classes frequencies (fre.maj^2/freq.each class)
    
    majCl <- names(which(table(dat[, tgt]) == max(table(dat[, tgt]))))
    
    if (length(majCl) == length(names)) {
      stop("Classes are balanced. Unable to reverse the frequencies!")
    }
    
    names.ove <- names[which(!(names %in% majCl))]
    
    # over-sample all the other classes reversing frequencies
    for (i in 1:length(names.ove)) { 
      Exs <- which(dat[, tgt] == names.ove[i])
      mmcl <- as.numeric(match(majCl, names))[1]
      n1 <- (li[[2]][mmcl])^2/li[[2]][as.numeric(match(names.ove[i], names))]
      n2 <- li[[2]][as.numeric(match(names.ove[i], names))]
      sel <- sample(Exs,
                    round(n1 - n2, 0),
                    replace = repl)
      newdata <- rbind(newdata, dat[sel, ])
    }
  } else {
    stop("Please provide a list with classes to over-sample
         or alternatively provide 'balance' or 'extreme'.")
  }
  
  newdata
}
