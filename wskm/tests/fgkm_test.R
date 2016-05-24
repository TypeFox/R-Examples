# rm(list=ls())

library("wskm")
x_short <- scale(fgkm.sample)
strGroup_short <- "0-9;10-19;20-49"
strGroup_short4 <- "10-19;0-9;20-49"
group_short <- c(rep(0, 10), rep(1, 10), rep(2, 30))
group_short3 <- c(rep("a", 10), rep("b", 10), rep("c", 30))
group_short4 <- c(rep(3, 10), rep(2, 10), rep(4, 30))

f1_short <- fgkm(x_short, 3, strGroup_short, 3, 1, seed = 19)
f2_short <- fgkm(x_short, 3, strGroup_short, 3, 1, seed = 19)
F1_short <- fgkm(x_short, 3, group_short, 3, 1, seed = 19)
F2_short <- fgkm(x_short, 3, group_short, 3, 1, seed = 19)
all.equal(f1_short, f2_short)
all.equal(F1_short, F2_short)
all.equal(F1_short, f1_short)

f1_short3 <- fgkm(x_short, 3, strGroup_short, 3, 1, seed = 19)
f2_short3 <- fgkm(x_short, 3, strGroup_short, 3, 1, seed = 19)
F1_short3 <- fgkm(x_short, 3, group_short3, 3, 1, seed = 19)
F2_short3 <- fgkm(x_short, 3, group_short3, 3, 1, seed = 19)
all.equal(f1_short3, f2_short3)
all.equal(F1_short3, F2_short3)
all.equal(F1_short3, f1_short3)

f1_short4 <- fgkm(x_short, 3, strGroup_short4, 3, 1, seed = 19)
f2_short4 <- fgkm(x_short, 3, strGroup_short4, 3, 1, seed = 19)
F1_short4 <- fgkm(x_short, 3, group_short4, 3, 1, seed = 19)
F2_short4 <- fgkm(x_short, 3, group_short4, 3, 1, seed = 19)
all.equal(f1_short4, f2_short4)
all.equal(F1_short4, F2_short4)
all.equal(F1_short4, f1_short4)
