### R code from vignette source 'dawai.Rnw'

###################################################
### code chunk number 1: dawai.Rnw:5-6 (eval = FALSE)
###################################################
## library(dawai)


###################################################
### code chunk number 2: biologicalApplication (eval = FALSE)
###################################################
## library("survival")
## data("pbc")
## data <- pbc[, c("bili", "albumin", "platelet", "stage")]
## data <- data[rowSums(is.na(data)) == 0, ]
## data <- cbind(data[, "stage", drop = FALSE],
##               "logBili" = log(data[["bili"]]),
##               "logAlbumin" = log(data[["albumin"]]),
##               "logPlatelet" = log(data[["platelet"]]))
## data$stage <- as.factor(data$stage)
## levels(data$stage)
## table(data$stage)
## levels(data$stage) <- c(1, 1, 2, 3)
## table(data$stage)
## A <- matrix(0, ncol = 9, nrow = 6)
## A[t(matrix(c(1, 1, 4, 4, 2, 5, 3, 6, 5, 8, 6, 9), nrow = 2))] <- 1
## A[t(matrix(c(1, 4, 4, 7, 2, 2, 3, 3, 5, 5, 6, 6), nrow = 2))] <- -1
## A
## set.seed(-5436)
## values <- runif(dim(data)[1])
## trainsubset <- (values < 0.25)
## testsubset <- (values >= 0.25)
## obj <- rlda(stage ~ logBili + logAlbumin + logPlatelet, data,
##            subset = trainsubset, gamma = c(0, 0.75, 1),
##            resmatrix = A, prior = c(1/3, 1/3, 1/3))
## obj
## pred <- predict(obj, newdata = data[testsubset,],
##               grouping = data[testsubset, "stage"])
## pred$error.rate
## err.est(obj)


###################################################
### code chunk number 3: patternRecognitionApplication (eval = FALSE)
###################################################
## data("Vehicle2")
## data <- Vehicle2[, 1:4]
## grouping <- Vehicle2$Class
## levels(grouping)
## levels(grouping) <- c(4, 2, 1, 3)
## set.seed(-9152)
## values <- runif(dim(data)[1])
## trainsubset <- (values < 0.25)
## obj <- rqda(data, grouping, subset = trainsubset, restext = "s>1,2,3")
## obj
## testsubset <- (values >= 0.25)
## pred <- predict(obj, newdata = data[testsubset,],
##               grouping = grouping[testsubset])
## pred$error.rate
## err.est(obj)


