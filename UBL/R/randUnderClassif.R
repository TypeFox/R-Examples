# ===================================================
# Random under-sampling strategy for multiclass problems.
# Basically randomly removes a percentage of cases of the class(es) 
# selected by the user. Alternatively, it can either balance all the 
# existing classes or it can "smoothly invert" the frequency
# of the examples in each class.
# Examples:
#   ir<- iris[-c(95:130), ]
#   myunder.iris <- RandUnderClassif(Species~., ir, list(setosa = 0.5, 
#                                                        versicolor = 0.8))
#   undBalan.iris <- RandUnderClassif(Species~., ir, "balance")
#   undInvert.iris <- RandUnderClassif(Species~., ir, "extreme")
# 
#   library(DMwR)
#   data(algae)
#   classes autumn and spring remain unchanged
#   C.perc = list(autumn = 1, summer = 0.9, winter = 0.4)
#   myunder.algae <- RandUnderClassif(season~., algae, C.perc)
#   undBalan.algae <- RandUnderClassif(season~., algae, "balance")
#   undInvert.algae <- RandUnderClassif(season~., algae, "extreme")
#   
#   library(MASS)
#   data(cats)
#   myunder.cats <- RandUnderClassif(Sex~., cats, list(M = 0.8))
#   undBalan.cats <- RandUnderClassif(Sex~., cats, "balance")
#   undInvert.cats <- RandUnderClassif(Sex~., cats, "extreme")
# 
# P. Branco, Mar 2015, Apr 2016
# ---------------------------------------------------
RandUnderClassif <- function(form, dat, C.perc = "balance", repl = FALSE)
  # Args:
  # form   a model formula
  # dat   the original training set (with the imbalanced distribution)
  # C.perc is a named list containing each class under-sampling percentage
  #        (between 0 and 1).
  #        The user may only provide the classes where he wants to apply 
  #        under-sampling. Alternatively it may be "balance" (the default)
  #        or "extreme", cases where the under-sampling percentages
  #        are automatically estimated
  # repl   is it allowed or not to perform sampling with replacement
  #
  # Returns: a new data frame modified through the random 
  #          under-sampling strategy
  
{
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  names <- sort(unique(dat[, tgt]))
  li <- class.freq(dat, tgt)
  
  if (is.list(C.perc)) { 
    # the under-sampling percentages are provided by the user
    if (any(C.perc > 1)) {
      stop("percentages provided must be < 1 to perform under-sampling!")
    }
    names.und <- names(which(C.perc < 1))

    # include examples from classes unchanged
    newdata <- dat[which(dat[, tgt] %in% names[which(!(names %in% names.und))]), ]
  
    for (i in 1:length(names.und)) { # under-sampling each class provided
      Exs <- which(dat[,tgt] == names.und[i])
      sel <- sample(Exs,
                    as.integer(C.perc[[names.und[i]]] * length(Exs)),
                    replace = repl)
      newdata <- rbind(newdata, dat[sel, ])
    }
  } else if (C.perc == "balance") { 
  # the under-sampling percentages must be calculated
    minCl <- names(which(table(dat[, tgt]) == min(table(dat[, tgt]))))
    if (length(minCl) == length(names)) {
      stop("Classes are already balanced!")
    }
    # add the cases of the minority classes
    minExs <- which(dat[, tgt] %in% minCl)
    newdata <- dat[minExs, ]
    names.und <- names[which(!(names %in% minCl))]
    
    # under-sample all the other classes
    for (i in 1:length(names.und)) { 
      Exs <- which(dat[, tgt] == names.und[i])
      sel <- sample(Exs,
                    as.integer(li[[2]][as.numeric(match(minCl, names))[1]]),
                    replace = repl)
      newdata <- rbind(newdata, dat[sel, ])
    }      
  }else if (C.perc == "extreme") {
    #"reverse" the classes frequencies (freq.min^2/freq. each class)
    minCl <- names(which(table(dat[, tgt]) == min(table(dat[, tgt]))))
    if (length(minCl) == length(names)) {
      stop("Classes are balanced. Unable to reverse the frequencies!")
    }
    # add the cases of the minority classes
    minExs <- which(dat[, tgt] %in% minCl)
    newdata <- dat[minExs, ]
    names.und <- names[which(!(names %in% minCl))]
    
    # under-sample all the other classes reversing frequencies 
    for(i in 1:length(names.und)){ 
      Exs <- which(dat[, tgt] == names.und[i])
      mmcl <- as.numeric(match(minCl, names))
      num1 <- (li[[2]][mmcl[1]])^2/li[[2]][as.numeric(match(names.und[i], 
                                                            names))]
      sel <- sample(Exs,
                    as.integer(num1),
                    replace = repl)
      newdata <- rbind(newdata, dat[sel, ])
    }      
  } else {
    stop("Please provide a list with classes to under-sample 
         or alternative specify 'balance' or 'extreme'.")
  }
  
  newdata
}

# ===================================================
# Auxiliar function which returns a list with the classes names
# and frequency of a data set
# P.Branco, Mar 2015
# ---------------------------------------------------

class.freq <- function(dat, tgt){
  names <- sort(unique(dat[, tgt]))
  li <- list(names, 
             sapply(names, 
                    function(x) length(which(dat[, tgt] == names[x]))))
  li
}

  