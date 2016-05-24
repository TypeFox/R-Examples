# rm(list=ls())

library("wskm")
x_short <- scale(fgkm.sample)
strGroup_short <- "0-9;10-19;20-49"
strGroup_short4 <- "10-19;0-9;20-49"
group_short <- c(rep(0, 10), rep(1, 10), rep(2, 30))
group_short3 <- c(rep("a", 10), rep("b", 10), rep("c", 30))
group_short4 <- c(rep(3, 10), rep(2, 10), rep(4, 30))

t1_short <- twkm(x_short, 3, strGroup_short, 3, 1, seed = 19)
t2_short <- twkm(x_short, 3, strGroup_short, 3, 1, seed = 19)
T1_short <- twkm(x_short, 3, group_short, 3, 1, seed = 19)
T2_short <- twkm(x_short, 3, group_short, 3, 1, seed = 19)
all.equal(t1_short, t2_short)
all.equal(T1_short, T2_short)
all.equal(T1_short, t1_short)

t1_short3 <- twkm(x_short, 3, strGroup_short, 3, 1, seed = 19)
t2_short3 <- twkm(x_short, 3, strGroup_short, 3, 1, seed = 19)
T1_short3 <- twkm(x_short, 3, group_short3, 3, 1, seed = 19)
T2_short3 <- twkm(x_short, 3, group_short3, 3, 1, seed = 19)
all.equal(t1_short3, t2_short3)
all.equal(T1_short3, T2_short3)
all.equal(T1_short3, t1_short3)

t1_short4 <- twkm(x_short, 3, strGroup_short4, 3, 1, seed = 19)
t2_short4 <- twkm(x_short, 3, strGroup_short4, 3, 1, seed = 19)
T1_short4 <- twkm(x_short, 3, group_short4, 3, 1, seed = 19)
T2_short4 <- twkm(x_short, 3, group_short4, 3, 1, seed = 19)
all.equal(t1_short4, t2_short4)
all.equal(T1_short4, T2_short4)
all.equal(T1_short4, t1_short4)

