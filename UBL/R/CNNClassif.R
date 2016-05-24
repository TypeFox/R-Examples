## ===================================================
## Condensed Nearest Neighbors strategy for multiclass imbalanced problems.
# 
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae), ]
#   myCNN <- CNNClassif(season~., clean.algae, 
#                       Cl = c("summer", "spring", "winter"), dist = "HEOM")
#   CNN1 <- CNNClassif(season~., clean.algae, 
#                      Cl = "smaller", dist = "HEOM")
#   CNN2<- CNNClassif(season~., clean.algae, Cl = "summer", dist = "HEOM")
#   summary(myCNN[[1]]$season)
#   summary(CNN1[[1]]$season)
#   summary(CNN2[[1]]$season)
#   ir<- iris[-c(95:130),]
#   myCNN.iris <- CNNClassif(Species~., ir, Cl = c("setosa", "virginica"))
#   CNN.iris1 <- CNNClassif(Species~., ir, Cl = "smaller")
#   CNN.iris2 <- CNNClassif(Species~., ir, Cl = "versicolor")
#   summary(ir$Species)
#   summary(myCNN.iris[[1]]$Species)
#   summary(CNN.iris1[[1]]$Species)
#   summary(CNN.iris2[[1]]$Species)
# 
# 
#   library(MASS)
#   data(cats)
#   CNN.catsF <- CNNClassif(Sex~., cats, Cl = "F")
#   CNN.cats <- CNNClassif(Sex~., cats, Cl = "smaller")
# 
## P.Branco, April 2015 April 2016
## ---------------------------------------------------
CNNClassif <- function(form, dat, dist = "Euclidean", p = 2, Cl = "smaller")
#   
#   Args:
#   form a model formula
#   dat  the original training set (with the imbalanced distribution)
#   dist represents the distance function to be used for the kNN computation
#   p    a parameter used when the dist is set to "p-norm" which represents the 
#        used p.
#   Cl   is a vector with the names of the more important classes. Defaults to 
#        "smaller" which automatically decides which are the relevant classes.
#        In this case, all the classes that have frequency below 
#        (nr.examples)/(nr.classes) are considered important.
#
#  Returns: a data frame with the data modified through the CNN strategy.
#

{
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  
  classes <- levels(dat[, tgt])
  nrCl <- length(classes)
  
  if (Cl[[1]] == "smaller") { # define which is(are) the important class(es)
    Cl <- names(which(table(dat[, tgt]) < nrow(dat)/nrCl))
  }
  otherCl <- setdiff(levels(dat[, tgt]), Cl)   
  
  # construct a set with all examples from important classes and one 
  # random example from the other classes
  C.I <- dat[dat[, tgt] %in% Cl, ]
  for (i in 1:length(otherCl)) {
    C.I <- rbind(C.I, dat[sample(which(dat[, tgt] %in% otherCl[i]), 1), ])
  }
  test <- which(!(rownames(dat) %in% rownames(C.I)))
  ad.idx <- c()
  
  for (i in 1:length(test)) {
    neig <- neighbours(tgt, rbind(C.I, dat[test[i], ]), dist, p, k = 1)
    d.n <- c(rownames(C.I), rownames(dat[test[i], ]))
    if (dat[d.n[neig[nrow(neig)]], tgt] != dat[test[i], tgt]) {
      ad.idx <- c(ad.idx, test[i])
    }
  }
  
  C.I <- rbind(C.I, dat[ad.idx, ])
  return(list(C.I, Cl, otherCl))
  
}
