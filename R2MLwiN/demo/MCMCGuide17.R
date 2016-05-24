############################################################################
#     MLwiN MCMC Manual
#
# 17  Modelling Spatial Data . . . . . . . . . . . . . . . . . . . . . . 247
#
#     Browne, W.J. (2009) MCMC Estimation in MLwiN, v2.13. Centre for
#     Multilevel Modelling, University of Bristol.
############################################################################
#     R script to replicate all analyses using R2MLwiN
#
#     Zhang, Z., Charlton, C., Parker, R, Leckie, G., and Browne, W.J.
#     Centre for Multilevel Modelling, 2012
#     http://www.bristol.ac.uk/cmm/software/R2MLwiN/
############################################################################

# 17.1 Scottish lip cancer dataset . . . . . . . . . . . . . . . . . . . 247

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)

## openbugs executable
if (!exists("openbugs")) openbugs <- "C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe"
while (!file.access(openbugs, mode = 0) == 0 || !file.access(openbugs, mode = 1) == 0 || !file.access(openbugs, mode = 4) ==
  0) {
  cat("Please specify the path for the OpenBUGS executable:\n")
  openbugs <- scan(what = character(0), sep = "\n")
  openbugs <- gsub("\\", "/", openbugs, fixed = TRUE)
}

## winbugs executable
#winbugs="C:/Program Files (x86)/WinBUGS14/WinBUGS14.exe"

## Read lips1 data
data(lips1, package = "R2MLwiN")
summary(lips1)

# 17.2 Fixed effects models . . . . . . . . . . . . . . . . . . . . . . .248

lips1 <- lips1[with(lips1, order(neigh1, area, area)), ]

(mymodel <- runMLwiN(log(obs) ~ 1 + offset(offs), D = "Poisson", estoptions = list(EstM = 1, notation = "class"),
  data = lips1))

(mymodel <- runMLwiN(log(obs) ~ 1 + perc_aff + offset(offs), D = "Poisson", estoptions = list(EstM = 1, notation = "class"),
  data = lips1))

# 17.3 Random effects models . . . . . . . . . . . . . . . . . . . . . . 251

(mymodel <- runMLwiN(log(obs) ~ 1 + perc_aff + offset(offs) + (0 | neigh1) + (1 | area), D = "Poisson", estoptions = list(EstM = 1,
  notation = "class", mcmcMeth = list(iterations = 50000)), data = lips1))

# 17.4 A spatial multiple-membership (MM) model . . . . . . . . . . . . .252

(mymodel <- runMLwiN(log(obs) ~ 1 + perc_aff + offset(offs) + (1 | neigh1) + (1 | area), D = "Poisson", estoptions = list(mm = list(list(mmvar = list("neigh1",
  "neigh2", "neigh3", "neigh4", "neigh5", "neigh6", "neigh7", "neigh8", "neigh9", "neigh10", "neigh11"), weights = list("weight1",
  "weight2", "weight3", "weight4", "weight5", "weight6", "weight7", "weight8", "weight9", "weight10", "weight11")),
  NA, NA), EstM = 1, mcmcMeth = list(iterations = 50000)), data = lips1))
# 17.5 Other spatial models . . . . . . . . . . . . . . . . . . . . . . .255

# 17.6 Fitting a CAR model in MLwiN . . . . . . . . . . . . . . . . . . .255




mymodel <- runMLwiN(log(obs) ~ 0 + perc_aff + offset(offs) + (1 | area) + (0 | area), D = "Poisson", estoptions = list(car = list(list(carvar = list("neigh1",
  "neigh2", "neigh3", "neigh4", "neigh5", "neigh6", "neigh7", "neigh8", "neigh9", "neigh10", "neigh11"), weights = list("wcar1",
  "wcar2", "wcar3", "wcar4", "wcar5", "wcar6", "wcar7", "wcar8", "wcar9", "wcar10", "wcar11")), NA, NA), EstM = 1,
  mcmcMeth = list(iterations = 50000)), BUGO = c(version = 4, n.chains = 1, bugs = openbugs, OpenBugs = TRUE), data = lips1)

summary(mymodel)
# 17.7 Including exchangeable random effects . . . . . . . . . . . . . . 259

(mymodel <- runMLwiN(log(obs) ~ 0 + perc_aff + offset(offs) + (1 | area) + (1 | area), D = "Poisson", estoptions = list(car = list(list(carvar = list("neigh1",
  "neigh2", "neigh3", "neigh4", "neigh5", "neigh6", "neigh7", "neigh8", "neigh9", "neigh10", "neigh11"), weights = list("wcar1",
  "wcar2", "wcar3", "wcar4", "wcar5", "wcar6", "wcar7", "wcar8", "wcar9", "wcar10", "wcar11")), NA, NA), EstM = 1,
  mcmcMeth = list(iterations = 50000)), data = lips1))

# 17.8 Further reading on spatial modelling . . . . . . . . . . . . . . .260

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .128





############################################################################
