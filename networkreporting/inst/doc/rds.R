## ---- message=FALSE------------------------------------------------------
library(networkreporting)
library(plyr)
library(ggplot2) # we'll use qplot from ggplot2 for plots
theme_set(theme_minimal())

data.dir <- file.path("~", "Dropbox", "dennis_matt", "data", "processed")

# eventually, a sample RDS dataset will be included with the
# package. for now, this will only work on a few machines
# and only after gnsum-curitiba-goc-prep.r has been run
load(file.path(data.dir, "gnsum-curitiba-goc-data.RData"))

## ------------------------------------------------------------------------
## build the nomination chains based on the
## ids in the data
seed.ids <- c("1", "2", "3", "4", "5")
chains <- llply(seed.ids,
                ## NB: make.chain is not currently exported, so we
                ## have to specify the package namespace explicitly
                networkreporting:::make.chain,
                survey.data)

## ------------------------------------------------------------------------
## pick a set of traits and estimate the mixing
## info for them
these.traits <- c("netsize.bss.big")

mm <- networkreporting:::estimate.mixing(nonseed.survey.data, parent.data, these.traits)

## ------------------------------------------------------------------------
dd <- networkreporting:::estimate.degree.distns(survey.data,
                                                d.hat.vals="netsize.bss",
                                                traits=these.traits,
                                                keep.vars=c("total.alters",
                                                            "total.aware"))

## ------------------------------------------------------------------------
num.reps <- 100
boot.mc.dat <- networkreporting:::rds.mc.boot.draws(chains,
                                                    mm,
                                                    dd,
                                                    num.reps=num.reps)

## ------------------------------------------------------------------------
## use the RDS-II estimator on each bootstrap resample

## get estimates of the total number of alters and the total number of
## alters aware that the ego is a heavy drug user for each bootstrap rep
boot.mc.total.ests <- ldply(boot.mc.dat,
                            rdsII.estimator,
                            ## now we use the name of the resampled degree
                            ## (degree) instead of the original one
                            d.hat.vals="degree",
                            y.vals="total.alters",
                            missing="complete.obs")

boot.mc.aware.ests <- ldply(boot.mc.dat,
                            rdsII.estimator,
                            d.hat.vals="degree",
                            y.vals="total.aware",
                            missing="complete.obs")

boot.mc.ests <- data.frame(estimate=boot.mc.aware.ests[,1] / boot.mc.total.ests[,1],
                           estimator="mc")


## ------------------------------------------------------------------------
boot.chain.dat <- networkreporting:::rds.chain.boot.draws(chains,
                                                          mm,
                                                          dd,
                                                          num.reps=num.reps)

## ------------------------------------------------------------------------
## use the RDS-II estimator on each bootstrap resample

## get estimates of the total number of alters and the total number of
## alters aware that the ego is a heavy drug user for each bootstrap rep
boot.chain.total.ests <- ldply(boot.chain.dat,
                               rdsII.estimator,
                               ## now we use the name of the resampled degree
                               ## (degree) instead of the original one
                               d.hat.vals="degree",
                               y.vals="total.alters",
                               missing="complete.obs")

boot.chain.aware.ests <- ldply(boot.chain.dat,
                               rdsII.estimator,
                               d.hat.vals="degree",
                               y.vals="total.aware",
                               missing="complete.obs")

boot.chain.ests <- data.frame(estimate=boot.chain.aware.ests[,1] / boot.chain.total.ests[,1],
                              estimator="chain")


## ------------------------------------------------------------------------
both.ests <- rbind(boot.chain.ests, boot.mc.ests)

comp.plot <- ggplot(both.ests) +
             geom_density(aes(x=estimate, color=estimator, group=estimator))

