# ===================================================
# Neighborhood Cleaning Rule (NCL) is an under-sampling strategy 
# based in ENN which in it's original formulation attempts to clean
# the neighborhood of the class of interest (Laurikalla, 2001). NCLClassif is 
# suitable for multiclass problems.
# Examples:
#
# ir<- iris[-c(95:130), ]
# ir.M1 <- NCLClassif(Species~., ir, k = 3, dist = "Manhattan", Cl = "smaller")
# ir.Def <- NCLClassif(Species~., ir)
# ir.Ch <- NCLClassif(Species~., ir, k = 7, dist = "Chebyshev",
#                     Cl = "virginica")
# ir.Eu <- NCLClassif(Species~., ir, k = 5, Cl = c("versicolor", "virginica"))
# summary(ir$Species)
# summary(ir.M1$Species)
# summary(ir.Def$Species)
# summary(ir.Ch$Species)
# summary(ir.Eu$Species)
#
# P. Branco, April 2015 April 2016
# ---------------------------------------------------
NCLClassif <- function(form, dat, k = 3, dist = "Euclidean", 
                       p = 2, Cl = "smaller")
  # Args:
  # form  a model formula
  # dat the original training set (with the unbalanced distribution)
  # k    is the number of neighbors considered (should be odd to avoid ties)
  # dist represents the distance function to be used for the kNN 
  #      computation
  # p    a parameter used when the dist is set to "p-norm" which represents
  #      the used p.
  # Cl   is a vector indicating the names of the most important classes for 
  #      the user. 
  #      Alternatively, a string can be provided which can be set to "smaller".
  #      The default is "smaller" meaning that the smaller classes 
  #      (classes with #< (#original data)/(#classes)) are the most important.
  #
  # Returns: a data frame with the data modified through the NCL strategy

{
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }
  
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  classes <- levels(dat[, tgt])
  nrCl <- length(classes)

  # When Cl is "all" or "smaller" the Cl and otherCl classes 
  # must be determined. Alternatively the user may choose the
  # Cl class(es) and then only otherCl is determined.
  
  if (Cl[1] == "smaller") {
    Cl <- names(which(table(dat[, tgt]) < nrow(dat)/nrCl))
  }
    otherCl <- setdiff(classes, Cl)

  
  if (!length(Cl)) {
    res <- ENNClassif(form, dat, k, dist, p)
    warning("Only ENN applied because no minority class was detected!")
    return(res)
  }


  # first remove noisy data with ENN on every classes except the important
  # ones (Cl) ENN can not be applied when only one class is in the less 
  # important group
  if (length(otherCl) > 1) {
  A1 <- ENNClassif(form, dat[which(dat[, tgt] %in% otherCl), ],
                   k, dist, p)
  } else {
    A1 <- list(dat[which(dat[, tgt] %in% otherCl), ], c())
  }
  
  
  originalDat.Cl <- subset(dat, dat[, tgt] %in% Cl)
  
  # now clean the neighborhood of the important class(es)

  neig <- neighbours(tgt, dat, dist, p, k)
  
  # construct the classes in otherCl that can be further under-sampled
  minVal <- min(sapply(Cl, function(x)length(which(dat[, tgt] == x))))
  pmv <- 0.5  * minVal
  To.rem <- names(which(sapply(otherCl, 
                               function(x) length(which(dat[, tgt] == x))) > pmv))

  rm.idx <- c()
  for (class in Cl) {
    ind <- which(dat[, tgt] == class)
    for (i in ind) {
# remove the neighbors of index i if all the k neighbors
# have a != tgt class and the neighbors 
# belong to a class with at least 1/2 of the smaller important class
      if (sum(dat[neig[i, ], tgt] != class) == k) {
        for (indNeig in neig[i, ]) {
          if (dat[indNeig, tgt] %in% To.rem) {
            rm.idx <- c(rm.idx, indNeig)
          }
        }
      }
    } 
  }

  rm.idx <- unique(rm.idx)
  if (length(rm.idx)) {
    res <- dat[-rm.idx, ]
  } else {
    res <- dat
  }
  res <- res[setdiff(rownames(res),
                     rownames(dat[which(dat[, tgt] %in% otherCl), ])[A1[[2]]]), ]

  res
}
