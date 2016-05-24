############################################################################
#     MLwiN User Manual
#
# 12  Modelling Count Data . . . . . . . . . . . . . . . . . . . . . . . 181
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


# 12.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . .181

data(mmmec, package = "R2MLwiN")
summary(mmmec)

# 12.2 Fitting a simple Poisson model . . . . . . . . . . . . . . . . . .182

(mymodel1 <- runMLwiN(log(obs) ~ 1 + uvbi + offset(log(exp)), D = "Poisson", data = mmmec))

# 12.3 A three-level analysis . . . . . . . . . . . . . . . . . . . . . .184

(mymodel2 <- runMLwiN(log(obs) ~ 1 + offset(log(exp)) + (1 | nation) + (1 | region), D = "Poisson", estoptions = list(Meth = 0), 
  data = mmmec))

(mymodel3 <- runMLwiN(log(obs) ~ 1 + offset(log(exp)) + (1 | nation) + (1 | region), D = "Poisson", estoptions = list(Meth = 0, 
  nonlinear = c(N = 1, M = 2)), data = mmmec))

(mymodel4 <- runMLwiN(log(obs) ~ 1 + uvbi + offset(log(exp)) + (1 | nation) + (1 | region), D = "Poisson", estoptions = list(Meth = 0, 
  nonlinear = c(N = 1, M = 2)), data = mmmec))

# 12.4 A two-level model using separate country terms . . . . . . . . . .186

addmargins(with(mmmec, table(nation)))

contrasts(mmmec$nation, 9) <- diag(9)

(mymodel5 <- runMLwiN(log(obs) ~ 0 + nation + nation:uvbi + offset(log(exp)) + (1 | region), D = "Poisson", estoptions = list(Meth = 0, 
  nonlinear = c(N = 1, M = 2)), data = mmmec))

xb <- predict(mymodel5)

plot(mymodel5@data$uvbi, xb, xlab = "UV B radiation", ylab = "prediction", type = "n")
lines(mymodel5@data$uvbi[mymodel5@data$nationBelgium == 1], xb[mymodel5@data$nationBelgium == 1], col = 1)
lines(mymodel5@data$uvbi[mymodel5@data$nationW_Germany == 1], xb[mymodel5@data$nationW_Germany == 1], col = 2)
lines(mymodel5@data$uvbi[mymodel5@data$nationDenmark == 1], xb[mymodel5@data$nationDenmark == 1], col = 3)
lines(mymodel5@data$uvbi[mymodel5@data$nationFrance == 1], xb[mymodel5@data$nationFrance == 1], col = 4)
lines(mymodel5@data$uvbi[mymodel5@data$nationUK == 1], xb[mymodel5@data$nationUK == 1], col = 5)
lines(mymodel5@data$uvbi[mymodel5@data$nationItaly == 1], xb[mymodel5@data$nationItaly == 1], col = 6)
lines(mymodel5@data$uvbi[mymodel5@data$nationIreland == 1], xb[mymodel5@data$nationIreland == 1], col = 7)
lines(mymodel5@data$uvbi[mymodel5@data$nationLuxembourg == 1], xb[mymodel5@data$nationLuxembourg == 1], col = 8)
lines(mymodel5@data$uvbi[mymodel5@data$nationNetherlands == 1], xb[mymodel5@data$nationNetherlands == 1], col = 9)
legend(7, 0.7, c("belgium", "wgermany", "denmark", "france", "uk", "italy", "ireland", "luxembourg", "netherlands"), 
  lty = 1, col = 1:9)

# 12.5 Some issues and problems for discrete response models . . . . . . 190

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .190

############################################################################
