### test file for etmCIF.
### Really simple tests and comparison with etm

require(etm)

data(abortion)

from <- rep(0, nrow(abortion))
to <- abortion$cause
entry <- abortion$entry
exit <- abortion$exit
id <- 1:nrow(abortion)
data <- data.frame(id, from, to, entry, exit, group = abortion$group)

## Computation of the CIFs with etm
tra <- matrix(FALSE, 4, 4)
tra[1, 2:4] <- TRUE

cif.control <- etm(data[data$group == 0, ], c("0", "1", "2", "3"), 
                        tra, NULL, 0)
cif.exposed <- etm(data[data$group == 1, ], c("0", "1", "2", "3"), 
                        tra, NULL, 0)


## Computation of the CIFs with etmCIF

netm <- etmCIF(Surv(entry, exit, cause != 0) ~ group, abortion,
               etype = cause, failcode = 3)

### let's do some comparisons :-)

all.equal(trprob(cif.control, "0 3"), netm[[1]]$est["0", "3", ])
all.equal(trprob(cif.control, "0 2"), netm[[1]]$est["0", "2", ])
all.equal(trprob(cif.control, "0 1"), netm[[1]]$est["0", "1", ])

all.equal(trprob(cif.exposed, "0 3"), netm[[2]]$est["0", "3", ])
all.equal(trprob(cif.exposed, "0 2"), netm[[2]]$est["0", "2", ])
all.equal(trprob(cif.exposed, "0 1"), netm[[2]]$est["0", "1", ])


all.equal(trcov(cif.control, "0 3"), netm[[1]]$cov["0 3", "0 3", ])
all.equal(trcov(cif.control, "0 2"), netm[[1]]$cov["0 2", "0 2", ])
all.equal(trcov(cif.control, "0 1"), netm[[1]]$cov["0 1", "0 1", ])

all.equal(trcov(cif.exposed, "0 3"), netm[[2]]$cov["0 3", "0 3", ])
all.equal(trcov(cif.exposed, "0 2"), netm[[2]]$cov["0 2", "0 2", ])
all.equal(trcov(cif.exposed, "0 1"), netm[[2]]$cov["0 1", "0 1", ])


netm

## test on the summary
snetm <- summary(netm)

snetm

all.equal(unname(trprob(cif.control, "0 3")), snetm[[1]][[3]]$P)
all.equal(unname(trprob(cif.control, "0 2")), snetm[[1]][[2]]$P)
all.equal(unname(trprob(cif.control, "0 1")), snetm[[1]][[1]]$P)

all.equal(unname(trprob(cif.exposed, "0 3")), snetm[[2]][[3]]$P)
all.equal(unname(trprob(cif.exposed, "0 2")), snetm[[2]][[2]]$P)
all.equal(unname(trprob(cif.exposed, "0 1")), snetm[[2]][[1]]$P)

scif.control <- summary(cif.control, ci.fun = "cloglog")
scif.exposed <- summary(cif.exposed, ci.fun = "cloglog")

all.equal(scif.control[[3]]$lower, snetm[[1]][[3]]$lower)
all.equal(scif.control[[3]]$upper, snetm[[1]][[3]]$upper)

all.equal(scif.exposed[[3]]$lower, snetm[[2]][[3]]$lower)
all.equal(scif.exposed[[3]]$upper, snetm[[2]][[3]]$upper)


### test with factors in the input
abortion$status <- with(abortion, ifelse(cause == 2, "life birth",
                                         ifelse(cause == 1, "ETOP", "spontaneous abortion")))

abortion$status <- factor(abortion$status)

netm.factor <- etmCIF(Surv(entry, exit, status != "cens") ~ group, abortion,
                      etype = status, failcode = "spontaneous abortion")


all.equal(trprob(cif.control, "0 3"), netm.factor[[1]]$est["0", "spontaneous abortion", ])
all.equal(trprob(cif.control, "0 2"), netm.factor[[1]]$est["0", "life birth", ])

netm.factor

summary(netm.factor)

### test with group as a character vector
abortion$ttt <- with(abortion, ifelse(group == 0, "control", "exposed"))
abortion$ttt <- factor(abortion$ttt)

netm.ttt <- etmCIF(Surv(entry, exit, status != "cens") ~ ttt, abortion,
                   etype = status, failcode = "spontaneous abortion")

all.equal(trprob(cif.control, "0 3"), netm.ttt[[1]]$est["0", "spontaneous abortion", ])
all.equal(trprob(cif.control, "0 2"), netm.ttt[[1]]$est["0", "life birth", ])

netm.ttt

summary(netm.ttt)


### A couple of comparisons with simulated data
set.seed(1313)
time <- rexp(100)
to <- rbinom(100, 2, prob = c(1/3, 1/3, 1/3))
from <- rep(11, 100)
id <- 1:100
cov <- rbinom(100, 1, 0.5)

dat.s <- data.frame(id, time, from, to, cov)

traa <- matrix(FALSE, 3, 3)
traa[1, 2:3] <- TRUE

aa0 <- etm(dat.s[dat.s$cov == 0, ], c("11", "1", "2"), traa, "0", 0)
aa1 <- etm(dat.s[dat.s$cov == 1, ], c("11", "1", "2"), traa, "0", 0)
aa <- etm(dat.s, c("11", "1", "2"), traa, "0", 0)

test <- etmCIF(Surv(time, to != 0) ~ 1, dat.s, etype = to)

test.c <- etmCIF(Surv(time, to != 0) ~ cov, dat.s, etype = to)

all.equal(trprob(aa, "11 1"), test[[1]]$est["0", "1", ])
all.equal(trprob(aa, "11 2"), test[[1]]$est["0", "2", ])

all.equal(trprob(aa0, "11 1"), test.c[[1]]$est["0", "1", ])
all.equal(trprob(aa0, "11 2"), test.c[[1]]$est["0", "2", ])

all.equal(trprob(aa1, "11 1"), test.c[[2]]$est["0", "1", ])
all.equal(trprob(aa1, "11 2"), test.c[[2]]$est["0", "2", ])

test

test.c

summary(test)
summary(test.c)
