# ===================================================
# Tomek Links are can be used as an under-sampling strategy. Two examples form a 
# Tomek link if their class label is different and they are each other nearest 
# neighbor. After find the existing Tomek links two different under-sampling 
# strategies can be applied: to remove only the example of the larger class, 
# or to discard both examples. TomekClassif is suitable for multiclass problems.
# Examples:
# 
#   ir <- iris[-c(95:130), ]
#   ir1norm <- TomekClassif(Species~., ir, dist = "p-norm", p = 1, Cl = "all")
#   irMan <- TomekClassif(Species~., ir, dist = "Manhattan", 
#                         Cl = "all", rem = "maj") 
#   irEuc <- TomekClassif(Species~., ir)
#   irCheb <- TomekClassif(Species~., ir, dist = "Chebyshev",
#                         Cl = c("virginica", "setosa"))
#   irChebAll <- TomekClassif(Species~., ir, dist = "Chebyshev")
#   irChebMaj <- TomekClassif(Species~., ir, dist = "Chebyshev", rem = "maj")
#   irHVDM <- TomekClassif(Species~., ir, dist = "HVDM")
#   summary(ir1norm[[1]]$Species)
#   summary(irMan[[1]]$Species)
#   summary(irEuc[[1]]$Species)
#   summary(irCheb[[1]]$Species)
#   summary(irChebAll[[1]]$Species)
#   summary(irChebMaj[[1]]$Species)
#   summary(irHVDM[[1]]$Species)
# 
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae), ]
#   alg.HVDM1 <- TomekClassif(season~., clean.algae, dist = "HVDM", 
#                            Cl = c("winter", "spring"), rem = "both")
#   alg.HVDM2 <- TomekClassif(season~., clean.algae, 
#                           dist = "HVDM", rem = "maj")
#   # removes only examples from class summer which are the majority class
#   # in the link
#   alg.EuM <- TomekClassif(season~., clean.algae, dist = "HEOM", 
#                          Cl = "summer", rem = "maj") 
#   # removes only examples from class summer in every link they appear
#   alg.EuB <- TomekClassif(season~., clean.algae, dist = "HEOM",
#                          Cl = "summer", rem = "both") 
#   summary(clean.algae$season)
#   summary(alg.HVDM1[[1]]$season)
#   summary(alg.HVDM2[[1]]$season)
#   summary(alg.EuM[[1]]$season)
#   summary(alg.EuB[[1]]$season)
#   P. Branco, April 2015 April 2016
# ---------------------------------------------------
TomekClassif <- function(form, dat, dist = "Euclidean", p = 2, Cl = "all",
                         rem = "both")
  
  # Args:
  # form a model formula
  # dat the original training set (with the unbalanced distribution)
  # dist represents the distance function to be used for the kNN 
  #      computation
  # p    a parameter used when the dist is set to "p-norm" which represents 
  #      the used p.
  # Cl   is a vector with the names of the classes that should be 
  #      under-sampled. Defaults to "all" meaning that all classes
  #      are considered for under-sampling.
  # rem  represents the under-sampling technique applied. When set to "both", 
  #      both examples will be remove; when set to "maj" only the link that
  #      belongs to the majority class (or the class(es) provided in Cl)
  #      is removed. When "both" is selected and only one example in the Tomek
  #      Link belongs to a class selected to be under-sampled, then  
  #      only that example is removed.
  #
  # Returns:
  # a list with the cleaned dataframe and the indexes of the examples removed
  
{
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }
  
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))

  ClCode <- as.integer(as.factor(levels(dat[, tgt])))
  ClChar <- levels(dat[, tgt])
  Clcount <- as.vector(table(dat[, tgt]))
  
  if (Cl[[1]] == "all") {
  rm.code <- ClCode
  } else {
   rm.code <- as.integer(factor(Cl, levels = ClChar)) 
  }

  neig <- neighbours(tgt, dat, dist, p, k = 1)
  # determine the Tomek links
  TL.idx <- matrix(nrow = 0, ncol = 4)
  for (i in 1:nrow(neig)) {
    if (neig[neig[i,],] == i & 
        !(neig[i, ] %in% TL.idx[, 1]) & 
        dat[i, tgt] != dat[neig[i, ], tgt]) {
      TL.idx <- rbind(TL.idx, c(i, neig[i, ], dat[i, tgt], dat[neig[i, ], tgt]))
    }
  }
  
  rm.idx <- c()
    
  # check which examples should be removed to break the link 
  if (nrow(TL.idx)) {
    for (i in 1:nrow(TL.idx)) {
      tr <- which(TL.idx[i, 3:4] %in% rm.code)
      if (rem == "both" & length(tr)) {
        rm.idx <- c(rm.idx, TL.idx[i, tr])
      } else if (length(tr)) { # rem=="maj"
        # only one of the examples is candidate for removal but we must check if it belongs to the larger class
        if (length(tr) == 1) {
          if (which.max(Clcount[TL.idx[i, 3:4]]) == tr) {
            rm.idx <- c(rm.idx, TL.idx[i, tr])
          }
        } else {# two candidate examples but only one can be removed: chech which one belongs to the larger class
          rm.idx <- c(rm.idx, TL.idx[i, which.max(Clcount[TL.idx[i, 3:4]])]) 
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

  return(list(res, rm.idx))
}
