# ===================================================
# Edited Nearest Neighbors (ENN) rule is an under-sampling strategy 
# based in a cleaning technique which was developed by Wilson & Martinez, 
# 1972. ENNClassif is suitable for multiclass problems.
# Basically it removes cases (possibily from all classes) 
# which do not agree with the majority of the k neighbors.
# The Cl parameter allows for the user to specify which particular 
# classes should be under-sampled.
# Examples:
# 
#   ir <- iris[-c(95:130), ]
#   ir1norm <- ENNClassif(Species~., ir, k = 5, dist = "p-norm",
#                         p = 1, Cl = "all")
#   irManhat <- ENNClassif(Species~., ir, k = 5, dist = "Manhattan",
#                         Cl = "all") 
#   irEucl <- ENNClassif(Species~., ir)
#   irCheby <- ENNClassif(Species~., ir, k = 7, dist = "Chebyshev", 
#                        Cl = c("virginica", "setosa"))
#   irChebyAll <- ENNClassif(Species~., ir, k = 7, dist = "Chebyshev")
#   irHVDM <- ENNClassif(Species~., ir, k = 3, dist = "HVDM")
# 
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae), ]
#   alg.HVDMWinter <- ENNClassif(season~., clean.algae, k = 3, 
#                               dist = "HVDM", Cl = "winter")
#   alg.HVDM <- ENNClassif(season~., clean.algae, k = 1, dist = "HVDM")
# # the data set constains nominal and numeric attributes, therefore if 
# # the Euclidean distance is used (the default) an error will occur
# #
#   ENNClassif(season~., clean.algae, k = 3)
#   alg.HEOM <- ENNClassif(season~., clean.algae, k = 5, dist = "HEOM",
#                          Cl = c("winter", "summer"))
#   summary(clean.algae$season)
#   summary(alg.HVDMWinter[[1]]$season)
#   summary(alg.HVDM[[1]]$season)
#   summary(alg.HEOM[[1]]$season)

# P. Branco, Mar 2015 Apr 2016
# ---------------------------------------------------
ENNClassif <- function(form, dat, k = 3, dist = "Euclidean", 
                       p = 2, Cl = "all")
  # Args:
  # form   a model formula
  # dat    the original training set (with the unbalanced distribution)
  # k      is the number of neighbors considered (should be odd to avoid ties)
  # dist   represents the distance function to be used for the kNN 
  #        computation
  # p      a parameter used when the dist is set to "p-norm" which 
  #        represents the used p.
  # Cl     is a vector with the names of classes that should be under-sampled. 
  #        The default is "all" meaning that examples from all the
  #        existing classes can be removed
  # Returns: a list with the cleaned dataframe and the indexes of the examples
  #        removed 
  
{
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }
  
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  nom.at <- which(sapply(dat[, -tgt], is.numeric) == FALSE)
  
  if (!length(unique(dat[, tgt]))) {
    stop("The data set contains only one class!")
  }
  

    neig <- neighbours(tgt, dat, dist, p, k)
    rm.idx <- c()
    
    if (Cl[1] == "all") {
      for (i in 1:nrow(dat)) {
        if (sum(dat[neig[i, ], tgt] != dat[i, tgt]) > k/2) {
          rm.idx <- c(rm.idx, i)
        }
      }
    } else if (length(Cl)) {
      for (i in Cl) { 
        for (j in which(dat[, tgt] == i)) {
          if (sum(dat[neig[j, ], tgt] != i) > k/2) {
            rm.idx <- c(rm.idx, j)
          }
        }
      }
    }
    
    if (length(rm.idx)) {
    res <- dat[-rm.idx, ]
    } else {
      warning("There are no examples to remove!")
      res <- dat
    }

  # to ensure that all the classes have at least one example in the
  # final data set  
  if (length(unique(dat[, tgt])) > length(unique(res[, tgt]))) {
    Cdif <- setdiff(unique(dat[, tgt]), unique(res[, tgt]))
    ad.idx <- sapply(Cdif, function(x) sample(which(dat[, tgt] %in% x), 1))
    res <- rbind(res, dat[ad.idx, ])
    rm.idx <- setdiff(rm.idx, ad.idx)
  }

  return(list(res, rm.idx))
}
