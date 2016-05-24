## ---- eval=FALSE, install------------------------------------------------
#  # install.packages("devtools")
#  library(devtools)
#  #devtools::install_github("soodoku/guess")

## ----stndcor-------------------------------------------------------------
# Load library
library(guess)

# Generate some data without DK
pre_test <- data.frame(item1=c(1,0,0,1,0), item2=c(1,NA,0,1,0)) 
pst_test <-  pre_test + cbind(c(0,1,1,0,0), c(0,1,0,0,1))
lucky <- rep(.25, 2)

# Unadjusted Effect
# Treating Don't Know as ignorance
colMeans(nona(pst_test) - nona(pre_test))

# MCAR
colMeans(pst_test - pre_test, na.rm=T)

# Adjusted Effect
stndcor(pre_test, pst_test, lucky)

## ----transmat------------------------------------------------------------
# Without Don't Know
pre_test_var <- c(1,0,0,1,0,1,0) 
pst_test_var <- c(1,0,1,1,0,1,1)
print(transmat(pre_test_var, pst_test_var))

# With Don't Know
pre_test_var <- c(1,0,NA,1,"d","d",0,1,0)
pst_test_var <- c(1,0,1,"d",1,0,1,1,"d")
print(transmat(pre_test_var, pst_test_var))


## ----guesstimate---------------------------------------------------------
# load(system.file("data/alldat.rda", package = "guess"))
load("../data/alldat.rda")

# nitems
nitems <- length(alldat)/400

# Vectors of Names
t1 <- paste0("guess.t1", 1:nitems)
t2 <- paste0("guess.t2", 1:nitems)

transmatrix <- multi_transmat(alldat[,t1], alldat[,t2])

res <- guesstimate(transmatrix)

round(res$param.lca[,1:4], 3)

round(res$est.learning[1:4], 3)

# Guesstimate with DK
# load(system.file("data/alldat_dk.rda", package = "guess"))
load("../data/alldat_dk.rda")
transmatrix <- multi_transmat(alldat_dk[,t1], alldat_dk[,t2], force9=T)
res_dk <- guesstimate(transmatrix)

round(res_dk$param.lca[,1:4], 3)

round(res_dk$est.learning[1:4], 3)

## ----lca_err-------------------------------------------------------------

# Raw 
# Generate some data without DK
pre_test <- data.frame(item1=c(1,0,0,1,0), item2=c(1,NA,0,1,0)) 
pst_test <-  pre_test + cbind(c(0,1,1,0,0), c(0,1,0,0,1))
diff <- pst_test - pre_test
stnd_err <-  sapply(diff, function(x) sqrt(var(x, na.rm=T)/(length(x)-1)))

# Bootstrapped s.e.

# LCA model
lca_stnd_err <- guess_stnderr(alldat[,t1], alldat[,t2], 10)
sapply(lca_stnd_err, function(x) round(head(x, 1),3))

lca_dk_stnd_err <- guess_stnderr(alldat_dk[,t1], alldat_dk[,t2], 10)
sapply(lca_dk_stnd_err, function(x) round(head(x, 1),3))

## ----fit_lca-------------------------------------------------------------
fit <- fit_nodk(alldat[,t1], alldat[,t2], res$param.lca[4,], res$param.lca[1:3,])

print(fit[,1:4])

fit <- fit_dk(alldat_dk[,t1], alldat_dk[,t2], res_dk$param.lca[8,], res_dk$param.lca[1:7,], force9=TRUE)

print(fit[,1:4])

