# ===================================================
# One-sided selection strategy for multiclass imbalanced problems.
#
# Examples:
#   ir <- iris[-c(95:130), ]
#   ir1 <- OSSClassif(Species~., ir, dist = "HVDM")
#   ir2 <- OSSClassif(Species~., ir, dist = "p-norm", p = 3, Cl = "virginica")
#   ir3 <- OSSClassif(Species~., ir, start = "Tomek")
#   ir4 <- OSSClassif(Species~., ir)
#   summary(ir1$Species)
#   summary(ir2$Species)
#   summary(ir3$Species)
#   summary(ir4$Species)
# 
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae), ]
#   alg1 <- OSSClassif(season~., clean.algae, dist = "HVDM", 
#                      Cl = c("spring", "summer"))
#   alg2 <- OSSClassif(season~., clean.algae, dist = "HEOM", 
#                      Cl = c("spring", "summer"), start = "Tomek")
#   alg3 <- OSSClassif(season~., clean.algae, dist = "HVDM", start = "CNN")
#   alg4 <- OSSClassif(season~., clean.algae, dist = "HVDM", start = "Tomek")
#   alg5 <- OSSClassif(season~., clean.algae, dist = "HEOM", Cl = "winter")
#   summary(alg1$season)
#   summary(alg2$season)
#   summary(alg3$season)
#   summary(alg4$season)
#   summary(alg5$season)
#
# P. Branco, April 2015 April 2016
# ---------------------------------------------------
OSSClassif <- function(form, dat, dist = "Euclidean",
                       p = 2, Cl = "smaller", start = "CNN") 
  # Args:
  # form    a model formula
  # dat     the original training set (with the unbalanced distribution)
  # dist    represents the distance function to be used for the kNN computation
  # p       a parameter used when the dist is set to "p-norm" which represents
  #         the used p.
  # Cl      is a vector with the names of the more important classes. Defaults 
  #         to "smaller" which automatically decides which are the relevant 
  #         classes. In this case, all the classes that have frequency below 
  #         (nr. examples)/(nr. classes) are considered important.
  # start   is a string which determines which strategy (CNN or Tomek links) 
  #         should be performed first. If set to "CNN" (the default) this 
  #         strategy will be performed first and Tomek links are applied after.
  #         If set to "Tomek" the reverse order is applied.
  #
  # Returns: a new data frame with the data set modified through the One-sided
  #         Selection strategy

{
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }
  
  if (start == "CNN") {
  #  Obtain the reduced data set with CNN
  d1 <- CNNClassif(form, dat, dist = dist, p = p, Cl = Cl)
  
  d2 <- TomekClassif(form, d1[[1]], dist = dist, p = p, Cl = d1[[3]],
                     rem="both")
  return(d2[[1]])
  } else if (start == "Tomek") {
    # the column where the target variable is
    tgt <- which(names(dat) == as.character(form[[2]]))
    classes <- levels(dat[, tgt])
    nrCl <- length(classes)
    
    if (Cl[[1]] == "smaller") { # define which is(are) the important class(es)
      Cl <- names(which(table(dat[, tgt]) < nrow(dat)/nrCl))
    }
    otherCl <- setdiff(levels(dat[, tgt]), Cl)
    d1 <- TomekClassif(form, dat, dist = dist, p = p, Cl = otherCl,
                       rem = "both")
    d2 <- CNNClassif(form, d1[[1]], dist = dist, p = p, Cl = Cl)
    return(d2[[1]])
  } else {
    stop("start parameter must be CNN or Tomek!")
  }
}
