### R code from vignette source 'burkina-faso-females.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: sweave-setup
###################################################

options(width = 70, continue = " "#, useFancyQuotes="UTF-8"
        ,SweaveHooks=list(fig=function()
par(mar=c(5.1, 4.1, 1.1, 2.1)))
)
pdf.options(pointsize = 9)

library(reshape)

library(ggplot2)
theme_set(theme_grey(base_size = 9))    # pointsize for plot text

library(popReconstruct)



###################################################
### code chunk number 2: source-functions
###################################################

## source(file.path(.find.package("popReconstruct"), "source"
##                  ,"pop_reconstruction_functions.R")
##        )

## dir(file.path(.find.package("popReconstruct"), "R", "NAMESPACE"))




###################################################
### code chunk number 3: load-burkina-data
###################################################

## load(file.path(.find.package("popReconstruct"), "data", "burkina_faso_females.RData"))
data(burkina_faso_females)



###################################################
### code chunk number 4: fert-rate-init-ests
###################################################

burkina.faso.females$fertility.rates



###################################################
### code chunk number 5: baseline-count-init-ests
###################################################

burkina.faso.females$baseline.pop.counts



###################################################
### code chunk number 6: bkfas-do-the-recon
###################################################

## set the seed random for the random number generator
set.seed(1)
###
### The reconstruction:
###
##
## !!! WARNING: This takes over 24 hours !!!
##
## commented out --->|
## BKFem.Recon.MCMC <-
##     popRecon.sampler(## Size of the MCMC sample and burn in
##                      n.iter = 4E4,
##                      burn.in = 500,
##                      thin.by = 50,

##                      ## initial estimates and census counts
##                      mean.f = burkina.faso.females$fertility.rates,
##                      mean.s = burkina.faso.females$survival.proportions,
##                      mean.g = burkina.faso.females$migration.proportions,
##                      mean.b = burkina.faso.females$baseline.pop.counts,
##                      pop.data = burkina.faso.females$census.pop.counts,

##                      ## Metropolis proposal variances
##                      prop.vars = burkina.faso.prop.vars,
##                      verb=TRUE
##                      )
## save(BKFem.Recon.MCMC, file = "Burkina_Faso_Recon_RESULTS.RData")
## |<--- end comment
load(file = "Burkina_Faso_Recon_RESULTS.RData")



###################################################
### code chunk number 7: results-age-spec-fert
###################################################

apply(BKFem.Recon.MCMC$fert.rate.mcmc, 2, "quantile", c(0.025, 0.5, 0.975))[,1:7]



###################################################
### code chunk number 8: PLOT-age-specific-fert-rate
###################################################
getOption("SweaveHooks")[["fig"]]()

######################################################################
###
### Calculate posterior quantiles for age-specific fertility rate and
### plot
###
######################################################################

require(ggplot2)
require(gdata)
##
## Posterior quantiles
##
vital.chain <- BKFem.Recon.MCMC$fert.rate.mcmc
q.to.plot = c(0.025, 0.5, 0.975)
q.vital <- apply(vital.chain, 2, function(z) quantile(z, probs = q.to.plot))
dimnames(q.vital) <- list(as.character(q.to.plot), colnames(vital.chain))
##
## Age, year labels
##
colspl <- strsplit(colnames(vital.chain), ".", fixed = TRUE)
years <- unique(sapply(colspl, FUN = function(z) z[1]))
fert.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
fert.ages.numeric <- as.numeric(gsub("[^0-9]", "", fert.ages))
##
## Reshape data frame
##
qvit.melt <- melt(q.vital)
qvit.melt.col <- cbind(qvit.melt
                   ,expand.grid(quant = q.to.plot, ages = fert.ages.numeric
                               ,years = years)
                   )
##
## Initial estimates
##
nzfr <- BKFem.Recon.MCMC$alg.params$non.zero.fert.rows
vital.init.est <-
    BKFem.Recon.MCMC$fixed.params$mean.fert.rate[nzfr,]
vital.init.est.melt.col <-
    cbind(value = melt(vital.init.est)$value
          ,expand.grid(ages = fert.ages.numeric
                       ,years = years, quant = 5) # use quant=5 for init.est
        )
##
## Prepare data sets
##
alpha <- BKFem.Recon.MCMC$fixed.params$alpha.fert.rate
beta <- BKFem.Recon.MCMC$fixed.params$beta.fert.rate

qvit.melt.df <- t(q.vital)
colnames(qvit.melt.df) <-
    paste("fert.rate.", prettyNum(as.numeric(colnames(qvit.melt.df)) * 100)
          , "pctl", sep = "")
qvit.melt.df <-
    data.frame(qvit.melt.df
               ,age = sapply(strsplit(rownames(qvit.melt.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(qvit.melt.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )

vital.init.est <-
    melt(BKFem.Recon.MCMC$fixed.params$mean.fert.rate[nzfr,])
vital.init.est <-
    rename.vars(vital.init.est, from = c("X1", "X2", "value")
                ,to = c("age", "year", "fert.rate.50pctl"))
vital.init.est.melt.df <-
    data.frame(vital.init.est, fert.rate.97.5pctl = NA
               ,fert.rate.2.5pctl = NA
               ,legend = "init. est."
               )

plot.df <- rbind(vital.init.est.melt.df, qvit.melt.df)
plot.df$age <- as.numeric(plot.df$age)
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")


##
## Plot quantiles
##
print(
      ggplot(data = plot.df, aes(x = age, y = fert.rate.50pctl, color = legend)) +
      facet_wrap(~ year) +
      geom_line(size = 1) +
      geom_point() +
      geom_ribbon(aes(ymin = fert.rate.2.5pctl
                      ,ymax = fert.rate.97.5pctl, fill = legend), alpha = 0.15
                  ,color = NA
                  ) +
      ylab("fert. rate")
      )



###################################################
### code chunk number 9: PLOT-age-specific-surv-prop
###################################################
getOption("SweaveHooks")[["fig"]]()

######################################################################
###
### Calculate posterior quantiles for age-specific survival
### proportion and plot
###
######################################################################

require(ggplot2)
require(gdata)
##
## Posterior quantiles
##
vital.chain <- BKFem.Recon.MCMC$surv.prop.mcmc
q.to.plot = c(0.025, 0.5, 0.975)
q.vital <- apply(vital.chain, 2, function(z) quantile(z, probs = q.to.plot))
dimnames(q.vital) <- list(as.character(q.to.plot), colnames(vital.chain))
##
## Age, year labels
##
colspl <- strsplit(colnames(vital.chain), ".", fixed = TRUE)
years <- unique(sapply(colspl, FUN = function(z) z[1]))
surv.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
surv.ages.numeric <- as.numeric(gsub("[^0-9]", "", surv.ages))
##
## Reshape data frame
##
qvit.melt <- melt(q.vital)
qvit.melt.col <- cbind(qvit.melt
                   ,expand.grid(quant = q.to.plot, ages = surv.ages.numeric
                               ,years = years)
                   )
##
## Initial estimates
##
vital.init.est <-
    BKFem.Recon.MCMC$fixed.params$mean.surv.prop
vital.init.est.melt.col <-
    cbind(value = melt(vital.init.est)$value
          ,expand.grid(ages = surv.ages.numeric
                       ,years = years, quant = 5) # use quant=5 for init.est
        )
##
## Prepare data sets
##
alpha <- BKFem.Recon.MCMC$fixed.params$alpha.surv.prop
beta <- BKFem.Recon.MCMC$fixed.params$beta.surv.prop

qvit.melt.df <- t(q.vital)
colnames(qvit.melt.df) <-
    paste("surv.prop.", prettyNum(as.numeric(colnames(qvit.melt.df)) * 100)
          , "pctl", sep = "")
qvit.melt.df <-
    data.frame(qvit.melt.df
               ,age = sapply(strsplit(rownames(qvit.melt.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(qvit.melt.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )

vital.init.est <-
    melt(BKFem.Recon.MCMC$fixed.params$mean.surv.prop)
vital.init.est <-
    rename.vars(vital.init.est, from = c("X1", "X2", "value")
                ,to = c("age", "year", "surv.prop.50pctl"))
vital.init.est.melt.df <-
    data.frame(vital.init.est, surv.prop.97.5pctl = NA
               ,surv.prop.2.5pctl = NA
               ,legend = "init. est."
               )

plot.df <- rbind(vital.init.est.melt.df, qvit.melt.df)
plot.df$age <- as.numeric(plot.df$age)
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")
##
## Plot quantiles
##
print(
      ggplot(data = plot.df, aes(x = age, y = surv.prop.50pctl, color = legend)) +
      facet_wrap(~ year) +
      geom_line(size = 1) +
      geom_point() +
      geom_ribbon(aes(ymin = surv.prop.2.5pctl
                      ,ymax = surv.prop.97.5pctl, fill = legend), alpha = 0.15
                  ,color = NA) +
      ylab("surv. prop")
      )



###################################################
### code chunk number 10: PLOT-age-specific-mig-prop
###################################################
getOption("SweaveHooks")[["fig"]]()

######################################################################
###
### Calculate posterior quantiles for age-specific migration
### proportion and plot
###
######################################################################

require(ggplot2)
require(gdata)
##
## Posterior quantiles
##
vital.chain <- BKFem.Recon.MCMC$mig.prop.mcmc
q.to.plot = c(0.025, 0.5, 0.975)
q.vital <- apply(vital.chain, 2, function(z) quantile(z, probs = q.to.plot))
dimnames(q.vital) <- list(as.character(q.to.plot), colnames(vital.chain))
##
## Age, year labels
##
colspl <- strsplit(colnames(vital.chain), ".", fixed = TRUE)
years <- unique(sapply(colspl, FUN = function(z) z[1]))
mig.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
mig.ages.numeric <- as.numeric(gsub("[^0-9]", "", mig.ages))
##
## Reshape data frame
##
qvit.melt <- melt(q.vital)
qvit.melt.col <- cbind(qvit.melt
                   ,expand.grid(quant = q.to.plot, ages = mig.ages.numeric
                               ,years = years)
                   )
##
## Initial estimates
##
vital.init.est <-
    BKFem.Recon.MCMC$fixed.params$mean.mig.prop
vital.init.est.melt.col <-
    cbind(value = melt(vital.init.est)$value
          ,expand.grid(ages = mig.ages.numeric
                       ,years = years, quant = 5) # use quant=5 for init.est
        )
##
## Prepare data sets
##
alpha <- BKFem.Recon.MCMC$fixed.params$alpha.mig.prop
beta <- BKFem.Recon.MCMC$fixed.params$beta.mig.prop

qvit.melt.df <- t(q.vital)
colnames(qvit.melt.df) <-
    paste("mig.prop.", prettyNum(as.numeric(colnames(qvit.melt.df)) * 100)
          , "pctl", sep = "")
qvit.melt.df <-
    data.frame(qvit.melt.df
               ,age = sapply(strsplit(rownames(qvit.melt.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(qvit.melt.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )

vital.init.est <-
    melt(BKFem.Recon.MCMC$fixed.params$mean.mig.prop)
vital.init.est <-
    rename.vars(vital.init.est, from = c("X1", "X2", "value")
                ,to = c("age", "year", "mig.prop.50pctl"))
vital.init.est.melt.df <-
    data.frame(vital.init.est, mig.prop.97.5pctl = NA
               ,mig.prop.2.5pctl = NA
               ,legend = "init. est."
               )

plot.df <- rbind(vital.init.est.melt.df, qvit.melt.df)
plot.df$age <- as.numeric(plot.df$age)
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")
##
## Plot quantiles
##
print(
      ggplot(data = plot.df, aes(x = age, y = mig.prop.50pctl, color = legend)) +
      facet_wrap(~ year) +
      geom_line(size = 1) +
      geom_point() +
      geom_ribbon(aes(ymin = mig.prop.2.5pctl
                      ,ymax = mig.prop.97.5pctl, fill = legend), alpha = 0.15
                  ,color = NA) +
      ylab("mig. prop")
      )



###################################################
### code chunk number 11: PLOT-age-specific-baseline-count
###################################################
getOption("SweaveHooks")[["fig"]]()

######################################################################
###
### Calculate posterior quantiles for age-specific baseline count
### and plot
###
######################################################################

require(ggplot2)
require(gdata)
##
## Posterior quantiles
##
vital.chain <- BKFem.Recon.MCMC$baseline.count.mcmc
q.to.plot = c(0.025, 0.5, 0.975)
q.vital <- apply(vital.chain, 2, function(z) quantile(z, probs = q.to.plot))
dimnames(q.vital) <- list(as.character(q.to.plot), colnames(vital.chain))
##
## Age, year labels
##
colspl <- strsplit(colnames(vital.chain), ".", fixed = TRUE)
years <- unique(sapply(colspl, FUN = function(z) z[1]))
baseline.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
baseline.ages.numeric <- as.numeric(gsub("[^0-9]", "", baseline.ages))
##
## Reshape data frame
##
qvit.melt <- melt(q.vital)
qvit.melt.col <- cbind(qvit.melt
                   ,expand.grid(quant = q.to.plot, ages = baseline.ages.numeric
                               ,years = years)
                   )
##
## Initial estimates
##
vital.init.est <-
    BKFem.Recon.MCMC$fixed.params$mean.baseline.count
vital.init.est.melt.col <-
    cbind(value = melt(vital.init.est)$value
          ,expand.grid(ages = baseline.ages.numeric
                       ,years = years, quant = 5) # use quant=5 for init.est
        )
##
## Prepare data sets
##
alpha <- BKFem.Recon.MCMC$fixed.params$alpha.population.count
beta <- BKFem.Recon.MCMC$fixed.params$beta.population.count

qvit.melt.df <- t(q.vital)
colnames(qvit.melt.df) <-
    paste("baseline.count.", prettyNum(as.numeric(colnames(qvit.melt.df)) * 100)
          , "pctl", sep = "")
qvit.melt.df <-
    data.frame(qvit.melt.df
               ,age = sapply(strsplit(rownames(qvit.melt.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(qvit.melt.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )

vital.init.est <-
    melt(BKFem.Recon.MCMC$fixed.params$mean.baseline.count)
vital.init.est <-
    rename.vars(vital.init.est, from = c("X1", "X2", "value")
                ,to = c("age", "year", "baseline.count.50pctl"))
vital.init.est.melt.df <-
    data.frame(vital.init.est, baseline.count.97.5pctl = NA
               ,baseline.count.2.5pctl = NA
               ,legend = "init. est."
               )

plot.df <- rbind(vital.init.est.melt.df, qvit.melt.df)
plot.df$age <- as.numeric(plot.df$age)
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")
##
## Plot quantiles
##
print(
      ggplot(data = plot.df, aes(x = age, y = baseline.count.50pctl, color = legend)) +
      facet_wrap(~ year) +
      geom_line(size = 1) +
      geom_point() +
      geom_ribbon(aes(ymin = baseline.count.2.5pctl
                      ,ymax = baseline.count.97.5pctl, fill = legend), alpha = 0.15
                  ,color = NA) +
      ylab("baseline. count")
      )



###################################################
### code chunk number 12: calculate-tfr-posterior
###################################################
getOption("SweaveHooks")[["fig"]]()

######################################################################
###
### Calculate posterior quantiles for TFR and plot
###
######################################################################

require(ggplot2)
q.to.plot = c(0.025, 0.5, 0.975)


###
### Posterior
###
dn <- list(NULL,
            unique(sapply(strsplit(colnames(BKFem.Recon.MCMC$fert.rate.mcmc)
                                   ,"\\."), FUN = function(z) z[[1]])
                   )
            )
BKFem.Recon.tfr <-
  matrix(0, nrow = nrow(BKFem.Recon.MCMC$fert.rate.mcmc)
         ,ncol = length(dn[[2]])
         ,dimnames = dn
         )

fert.rate.mcmc.colYrs <-
  sapply(strsplit(colnames(BKFem.Recon.MCMC$fert.rate.mcmc)
                         ,"\\."), FUN = function(z) z[[1]])
##
## calculate tfr
##
for(i in 1:ncol(BKFem.Recon.tfr)) {
  colYrs.index <- fert.rate.mcmc.colYrs == colnames(BKFem.Recon.tfr)[i]
  BKFem.Recon.tfr[,i] <-
    apply(BKFem.Recon.MCMC$fert.rate.mcmc[,colYrs.index]
        ,1
        ,FUN = function(z) 5 * sum(z)
        )
}
##
## tfr quantiles
##
BKFem.Recon.tfrQuant <- apply(BKFem.Recon.tfr, 2, FUN = function(z)
                        {
                          quantile(z, probs = q.to.plot)
                        })
BKFem.Recon.tfrQuant.df <-
    as.data.frame(t(BKFem.Recon.tfrQuant))
colnames(BKFem.Recon.tfrQuant.df) <-
    paste("tfr.", strsplit(colnames(BKFem.Recon.tfrQuant.df), split = "%")
          ,"pctl", sep = "")
BKFem.Recon.tfrQuant.df$legend = "posterior"
BKFem.Recon.tfrQuant.df$year = as.numeric(rownames(BKFem.Recon.tfrQuant.df))

###
### Initial estimates
###
BKFem.Recon.tfr.init.est <-
    data.frame(year = colnames(BKFem.Recon.MCMC$fixed.params$mean.fert.rate)
               ,tfr.50pctl =
               melt(5 * colSums(BKFem.Recon.MCMC$fixed.params$mean.fert.rate))[,1]
               )
BKFem.Recon.tfr.init.est$legend <- "init. est."
BKFem.Recon.tfr.init.est$tfr.2.5pctl <-
    BKFem.Recon.tfr.init.est$tfr.97.5pctl <- NA


###
### Plot
###
plot.df <-
    rbind(BKFem.Recon.tfrQuant.df, BKFem.Recon.tfr.init.est)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")
plot.df$year <- as.numeric(plot.df$year)
print(ggplot(data = plot.df, aes(x = year, y = tfr.50pctl, color = legend)) +
      geom_line(size = 1) +
      geom_point() +
      geom_point() + geom_ribbon(aes(ymin = tfr.2.5pctl
                                     ,ymax = tfr.97.5pctl, fill = legend)
                                 ,alpha = 0.15, color = NA) +
      ylab("total fertility rate")
      )



###################################################
### code chunk number 13: calculate-e0-posterior
###################################################
getOption("SweaveHooks")[["fig"]]()

######################################################################
###
### Calculate posterior quantiles for life expectancy at birth and
### plot
###
######################################################################

require(ggplot2)
q.to.plot = c(0.025, 0.5, 0.975)

surv.prop.years <-
    sapply(strsplit(colnames(BKFem.Recon.MCMC$surv.prop.mcmc), "\\."), "[[", 1)

message("Calculating life expectancy at birth ...")
BKFem.leb.stationary.df <-
    apply(BKFem.Recon.MCMC$surv.prop.mcmc[,], 1, function(z) {
        tapply(z, INDEX = surv.prop.years, FUN = "life.expectancy.stationary")
      })
message("... done")

BKFem.leb.stationary.Quantiles <-
    apply(BKFem.leb.stationary.df, 1, "quantile", probs = q.to.plot)

BKFem.leb.stationary.Quantiles.df <-
    as.data.frame(t(BKFem.leb.stationary.Quantiles))

colnames(BKFem.leb.stationary.Quantiles.df) <-
    paste("leb."
          , strsplit(colnames(BKFem.leb.stationary.Quantiles.df)
                                            , split = "%")
          ,"pctl", sep = "")

BKFem.leb.stationary.Quantiles.df$legend <- "posterior"
BKFem.leb.stationary.Quantiles.df$year <-
    as.numeric(rownames(BKFem.leb.stationary.Quantiles.df))


###
### Prior by converting posterior quantiles of survival and assuming
### stationary population relation holds
###
lebp.yrs <- as.numeric(colnames(BKFem.Recon.MCMC$fixed.params$mean.surv.prop))
BKFem.lebPrior.stationary.df <-
    data.frame(year = lebp.yrs
               ,leb.50pctl = apply(BKFem.Recon.MCMC$fixed.params$mean.surv.prop
                ,2
                ,FUN = "life.expectancy.stationary"
                ))
BKFem.lebPrior.stationary.df$leb.2.5pctl <-
    BKFem.lebPrior.stationary.df$leb.97.5pctl <- NA
BKFem.lebPrior.stationary.df$legend <- "init. est."

###
### Plot
###
plot.df <-
    rbind(BKFem.lebPrior.stationary.df, BKFem.leb.stationary.Quantiles.df)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")
print(ggplot(data = plot.df, aes(x = year, y = leb.50pctl, color = legend)) +
      geom_line(alpha = 0.65) +
      geom_point() +
      geom_ribbon(aes(ymin = leb.2.5pctl
                                     ,ymax = leb.97.5pctl, fill = legend)
                                 ,alpha = 0.15, color = NA) +
      ylab("life expectancy at birth (years)")
      )



###################################################
### code chunk number 14: calculate-tot-avg-ann-net-mig
###################################################
getOption("SweaveHooks")[["fig"]]()

#######################################################################
###
### Calculate posterior quantiles for average annual total net
### number of migrants
###
######################################################################

require(ggplot2)
q.to.plot = c(0.025, 0.5, 0.975)

## NB: Can't simply sum migration proportions because they are based
##     on different population totals. Need to get net number of migrants
##     and convert back into proportions. Use Leslie matrix formula from
##     article draft.
##

###
### Posterior distribution
###

##
## Prepare output matrix
##
BKFem.Recon.netMig <-
  matrix(0, nrow = nrow(BKFem.Recon.MCMC$mig.prop.mcmc)
         ,ncol = ncol(BKFem.Recon.MCMC$fixed.params$mean.mig.prop)
         ,dimnames = list(NULL,
            unique(sapply(strsplit(colnames(BKFem.Recon.MCMC$mig.prop.mcmc)
                                   ,"\\."), FUN = function(z) z[[1]])
                   )
            )
         )

##
## The 5-year sub-intervals to be used as an index into the columns of
## BKFem.Recon.netMig
##
mig.prop.mcmc.colYrs <-
  sapply(strsplit(colnames(BKFem.Recon.MCMC$mig.prop.mcmc)
                         ,"\\."), FUN = function(z) z[[1]])
mig.prop.mcmc.colYrsUniq <- unique(mig.prop.mcmc.colYrs)

##
## Years used in survival proportions
##
surv.prop.mcmc.colYrs <-
    sapply(strsplit(colnames(BKFem.Recon.MCMC$surv.prop.mcmc)
                    ,"\\."), FUN = function(z) z[[1]])

##
## Concatenate baseline and lx to get a single matrix with population
## counts
##

pop.mat <- cbind(BKFem.Recon.MCMC$baseline.count.mcmc
                ,BKFem.Recon.MCMC$lx.mcmc)

##
## Index for population years
##

pop.mat.colYrs <- sapply(strsplit(colnames(pop.mat)
                         ,"\\."), FUN = function(z) z[[1]])
pop.mat.colYrsUniq <- unique(pop.mat.colYrs)

message("Calculating net number of migrants ...")
for(k in 1:nrow(BKFem.Recon.MCMC$mig.prop.mcmc)) {
    if(k %% 1000 == 0)
        message(paste("row ", k, " of "
                      ,nrow(BKFem.Recon.MCMC$mig.prop.mcmc), sep = "")
                )
    ##
    ## cycle through years
    for(i in 1:ncol(BKFem.Recon.netMig)) {
        ##
        ## 5-year sub-intervals for indexing columns
        mig.colYrs.index <-
            colnames(BKFem.Recon.netMig) == mig.prop.mcmc.colYrsUniq[i]
        surv.colYrs.index <-
            surv.prop.mcmc.colYrs == mig.prop.mcmc.colYrsUniq[i]
        fert.colYrs.index <-
            fert.rate.mcmc.colYrs == mig.prop.mcmc.colYrsUniq[i]
        pop.colYrs.index1 <-
            pop.mat.colYrs == mig.prop.mcmc.colYrsUniq[i]
        pop.colYrs.index2 <-
            pop.mat.colYrs == as.numeric(mig.prop.mcmc.colYrsUniq[i]) + 5
        ##
        ## get vital rates and make leslie matrix
        sk <- BKFem.Recon.MCMC$surv.prop.mcmc[k,surv.colYrs.index]
        fk <- rep(0, nrow(BKFem.Recon.MCMC$fixed.params$mean.fert.rate))
        fk[BKFem.Recon.MCMC$alg.params$non.zero.fert.rows] <-
            BKFem.Recon.MCMC$fert.rate.mcmc[k,fert.colYrs.index]
        popk1 <- pop.mat[k,pop.colYrs.index1]
        popk2 <- pop.mat[k,pop.colYrs.index2]
        Lk <- make.leslie.matrix(pop = popk1, surv = sk, fert = fk, srb = 1.05
                           ,age.int = 5)
        ##
        ## calculate net number of migrants
        netMigk <- net.number.migrants(n1 = popk1, n2 = popk2, L = Lk)
        ##
        ## store
        BKFem.Recon.netMig[k, mig.colYrs.index] <- sum(netMigk)
    }
}
message("... done")

##
## Posterior quantiles
##
BKFem.nmig.post.quant <-
    apply(BKFem.Recon.netMig, 2, FUN = function(z)
                        {
                          quantile(z, probs = q.to.plot)
                        })
BKFem.nmig.post.quant.df <-
    as.data.frame(t(BKFem.nmig.post.quant))

colnames(BKFem.nmig.post.quant.df) <-
    paste("total.mig.count.", strsplit(colnames(BKFem.nmig.post.quant.df)
                                            , split = "%")
          ,"pctl", sep = "")

BKFem.nmig.post.quant.df$legend <- "posterior"
BKFem.nmig.post.quant.df$year <-
    as.numeric(rownames(BKFem.nmig.post.quant.df))

###
### Initial estimates
###

##
## Prepare output matrix
##
BKFem.nmig.input <- rep(0, ncol(BKFem.Recon.MCMC$fixed.params$mean.mig.prop))
names(BKFem.nmig.input) <-
    colnames(BKFem.Recon.MCMC$fixed.params$mean.mig.prop)
##
## Input population counts
##
pop.input.mat <-
    popRecon.ccmp.female(pop=BKFem.Recon.MCMC$fixed.params$mean.baseline.count
                      ,surv=BKFem.Recon.MCMC$fixed.params$mean.surv.prop
                      ,fert=BKFem.Recon.MCMC$fixed.params$mean.fert.rate
                      ,mig=BKFem.Recon.MCMC$fixed.params$mean.mig.prop
                      )
##
## Calculate input net migration
##
for(k in 1:(ncol(pop.input.mat)-1)) {
    Lk <- make.leslie.matrix(pop = pop.input.mat[,k]
                       ,surv = BKFem.Recon.MCMC$fixed.params$mean.surv.prop[,k]
                       ,fert = BKFem.Recon.MCMC$fixed.params$mean.fert.rate[,k]
                       ,srb = 1.05
                       ,age.int = 5)
        netMigk <- net.number.migrants(n1 = pop.input.mat[,k]
                                ,n2 = pop.input.mat[,k+1]
                                ,L = Lk)
    BKFem.nmig.input[k] <- sum(netMigk)
}

BKFem.nmig.input.df <-
    data.frame(year = as.numeric(names(BKFem.nmig.input))
               ,total.mig.count.50pctl = BKFem.nmig.input
               )
BKFem.nmig.input.df$total.mig.count.2.5pctl <- NA
BKFem.nmig.input.df$total.mig.count.97.5pctl <- NA
BKFem.nmig.input.df$legend <- "init. est."

###
### Plot
###
plot.df <- rbind(BKFem.nmig.input.df, BKFem.nmig.post.quant.df)
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")
print(ggplot(data = plot.df, aes(x = year
             , y = total.mig.count.50pctl/1E3, color = legend)) +
      geom_line(alpha = 0.65) +
      geom_point() +
      geom_point() + geom_ribbon(aes(ymin = total.mig.count.2.5pctl/1E3
                                     ,ymax = total.mig.count.97.5pctl/1E3, fill = legend)
                                 ,alpha = 0.15, color = NA) +
      ylab("net number of migrants (000s)")
      )



