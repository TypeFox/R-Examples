############################################################################
#     MLwiN User Manual
#
# 4   Random Intercept and Random Slope Models . . . . . . . . . . . . . .47
#
#     Rasbash, J., Steele, F., Browne, W. J. and Goldstein, H. (2012).
#     A User's Guide to MLwiN, v2.26. Centre for Multilevel Modelling,
#     University of Bristol.
############################################################################
#     R script to replicate all analyses using R2MLwiN
#
#     Zhang, Z., Charlton, C., Parker, R, Leckie, G., and Browne, W.J.
#     Centre for Multilevel Modelling, 2012
#     http://www.bristol.ac.uk/cmm/software/R2MLwiN/
############################################################################

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)


# 4.1 Random intercept models . . . . . . . . . . . . . . . . . . . . . . 47

data(tutorial, package = "R2MLwiN")

plot(tutorial$standlrt, tutorial$normexam, asp = 1)

(mymodel1 <- runMLwiN(normexam ~ 1 + standlrt + (1 | student), data = tutorial))

(mymodel2 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(resi.store = TRUE), 
  data = tutorial))

# 4.2 Graphing predicted school lines from a random intercept model . . . 51

xb <- predict(mymodel2)

plot(tutorial$standlrt, xb, type = "l")

u0 <- mymodel2@residual$lev_2_resi_est_Intercept

xbu <- xb + u0[mymodel2@data$school]

head(u0)

plot(tutorial$standlrt, xbu, type = "l")

pred <- as.data.frame(cbind(mymodel2@data$school, mymodel2@data$standlrt, xbu)[order(mymodel2@data$school, mymodel2@data$standlrt), 
  ])

colnames(pred) <- c("school", "standlrt", "xbu")

xyplot(xbu ~ standlrt, type = "l", group = school, data = pred)

# 4.3 The effect of clustering on the standard errors of coeficients . . .58

(mymodel3 <- runMLwiN(normexam ~ 1 + standlrt + schgend + (1 | school) + (1 | student), data = tutorial))

(mymodel4 <- runMLwiN(normexam ~ 1 + standlrt + schgend + (1 | student), data = tutorial))

# 4.4 Does the coeficient of standlrt vary across schools? Introducing a random slope . . . . . . . . . . . . . .
# . . . . . . . . . . . . . .59

(mymodel5 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(resi.store = TRUE), 
  data = tutorial))

# 4.5 Graphing predicted school lines from a random slope model . . . . . 62

xb <- predict(mymodel5)

u <- cbind(mymodel5@residual$lev_2_resi_est_Intercept, mymodel5@residual$lev_2_resi_est_standlrt)

rphat <- rowSums(as.matrix(mymodel5@data[, c("Intercept", "standlrt")]) * as.matrix(u[tutorial$school, ]))

xbu <- xb + rphat

pred <- as.data.frame(cbind(mymodel5@data$school, mymodel5@data$standlrt, xbu)[order(mymodel5@data$school, mymodel5@data$standlrt), 
  ])

colnames(pred) <- c("school", "standlrt", "xbu")

xyplot(xbu ~ standlrt, type = "l", group = school, data = pred)


#     Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . 64

############################################################################
