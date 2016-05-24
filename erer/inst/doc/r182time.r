# A. Source a program
#    Program 15.7 A new function for OLS with the S4 mechanism
setwd("C:/aErer"); source("r157s4.r"); library(erer); data(daIns)

# B. Compare time used by lm2 and lm 
t0 <- proc.time()
  for (i in 1:100) {
    ny <- lm2(daIns, name.y = "Y", name.x = colnames(daIns)[-c(1, 14)])
  }
t1 <- proc.time(); ta <- t1 - t0

tb <- system.time(
  for (i in 1:100) {
    ny <- lm2(daIns, name.y = "Y", name.x = colnames(daIns)[-c(1, 14)])
  }
)

tc <- system.time(
  for (i in 1:100) {
    gg <- lm(formula = Y ~ 1 + Injury + HuntYrs + Nonres + Lspman + Lnong + 
      Gender + Age + Race + Marital + Edu + Inc + TownPop, data = daIns)
  }
)
total <- rbind(ta, tb, tc)
rownames(total) <- c("lm2 trial 1", "lm2 trial 2", "lm default"); total

# C. Recording time used by component in a function
# Turn profiling on and off
Rprof(filename = "r182Timefile.out", memory.profiling = TRUE)
  for (i in 1:100) {
    ny <- lm2(daIns, name.y = "Y", name.x = colnames(daIns)[-c(1, 14)])
  }
Rprof(NULL)

res <- summaryRprof(filename = "r182Timefile.out")
names(res); res[["by.total"]][1:10, ]

# D. Memory
gc()
object.size(ny)
lss(n = 4)
memory.profile(); memory.size(); memory.limit()