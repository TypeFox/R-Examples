# $Id: interval_test.R 169 2006-01-05 05:02:44Z andrewr $

rm(list=ls())

plots <- read.csv("../data/fia_plots.csv")

plots$forest <- factor(plots$forest)

reps <- 10

ls()

## The following code extracts the plot-level estimates from the FIA database.

#load("/home/andrewr/0.rosetta/history/1.Research/Equivalence/Data/fia.wo.VSEP16-03.RData")

#dim(fia.wo)

#names(fia.wo)

#fia.wo[1:10,]

#plot(fia.wo$DBH[1:100], fia.wo$BA[1:100])
#plot(fia.wo$DBH[1:100], fia.wo$BA[1:100]*fia.wo$TF[1:100])

#plots <- aggregate(x=list(ba=fia.wo$sba),
#                   by=list(forest=fia.wo$FORCODE,
#                     plot=fia.wo$PLOT),
#                   FUN=mean)

#table(plots$forest)

#save(plots, file="../data/plots.RData")

means <- tapply(plots$ba, plots$forest, mean)

require(boot)
require(survey)

boot.mean <- function(x, index)
  c(mean(x[index]), (length(x)-1)*var(x[index])/length(x)^2)

within <- function(conf, x)
  (x >= conf[length(conf)-1] & x < conf[length(conf)])
int.above <- function(conf, x)
  (x >= conf[length(conf)])
int.below <- function(conf, x)
  (x < conf[length(conf)-1])
int.length <- function(conf)
  conf[length(conf)] - conf[length(conf)-1]

## Goal: variance stabilization.  Square root might work.

opar <- par(mfrow=c(9,2), mar=c(5,5,0,0))
for (i in 1:length(levels(plots$forest))) {
  boot.object <- boot(plots$ba[plots$forest==levels(plots$forest)[i]],
                      boot.mean, R=1999)
  scatter.smooth(boot.object$t[,1], sqrt(boot.object$t[,2]), main="",
       ylab="Standard Error",
       xlab=paste("Bootstrap Mean, Forest:", levels(plots$forest)[i]))
  boot.object <- boot(sqrt(plots$ba[plots$forest==levels(plots$forest)[i]]),
                      boot.mean, R=1999)
  scatter.smooth(boot.object$t[,1], sqrt(boot.object$t[,2]), main="",
       ylab="Standard Error",
       xlab=paste("Bootstrap Mean, Forest:", levels(plots$forest)[i]))
}
par(opar)

## Sample with replacement 

n.s <- c(16)
store <- 0

total.reps <- reps * length(n.s) * length(levels(plots$forest))

interval.lengths <- list(rep(NA, total.reps))
interval.coverage <- list(rep(NA, total.reps))
interval.above <- list(rep(NA, total.reps))
interval.below <- list(rep(NA, total.reps))

for (j in 1:length(unique(n.s))) {            # j <- 1
  n <- n.s[j]
  for (i in 1:length(levels(plots$forest))) { # i <- 1
    for (b in 1:reps) {                       # b <- 1
      store <- store + 1
      the.sample <-
        sample(plots$ba[plots$forest==levels(plots$forest)[i]],
               size=n,
               replace=TRUE)
      boot.object <- boot(the.sample, boot.mean, R=1999)
      ci <- boot.ci(boot.object)
      log.ci <- boot.ci(boot.object, type=c("norm","basic","stud"),
                        h=log, hdot=function(x) 1/x, hinv=exp)
      sqrt.ci <- boot.ci(boot.object, type=c("norm","basic","stud"),
                        h=sqrt, hdot=function(x) 2*x, hinv=function(x) x^2)
      classical <- mean(the.sample) + c(-1, 1) *
        sd(the.sample) / sqrt(n) * qt(p=0.975, df=n-1)
      interval.coverage[[store]] <-
        c(within(classical, means[i]),
          as.numeric(lapply(ci[4:8], within, x=means[i])),
          as.numeric(lapply(log.ci[4:6], within, x=means[i])),
          as.numeric(lapply(sqrt.ci[4:6], within, x=means[i])))
      interval.lengths[[store]] <- 
        c(int.length(classical),
          as.numeric(lapply(ci[4:8], int.length)),
          as.numeric(lapply(log.ci[4:6], int.length)),
          as.numeric(lapply(sqrt.ci[4:6], int.length)))
      interval.above[[store]] <- 
        c(int.above(classical, means[i]),
          as.numeric(lapply(ci[4:8], int.above, x=means[i])),
          as.numeric(lapply(log.ci[4:6], int.above, x=means[i])),
          as.numeric(lapply(sqrt.ci[4:6], int.above, x=means[i])))
      interval.below[[store]] <- 
        c(int.below(classical, means[i]),
          as.numeric(lapply(ci[4:8], int.below, x=means[i])),
          as.numeric(lapply(log.ci[4:6], int.below, x=means[i])),
          as.numeric(lapply(sqrt.ci[4:6], int.below, x=means[i])))
    }
    cat("\nFinished simulation for Forest", levels(plots$forest)[i],
        "for n =", n, "(", reps, "reps )")
  }
  cat("\n")
}

output <- 
  expand.grid(interval=c("classical",
                "normal", "basic", "studentized", "percentile", "bca",
                "log.normal", "log.basic", "log.studentized",
                "sqrt.normal", "sqrt.basic", "sqrt.studentized"),
              replicate=1:reps,
              forest = levels(plots$forest),
              n = n.s)

output$length <- unlist(interval.lengths)
output$above <- unlist(interval.above) 
output$coverage <- unlist(interval.coverage) 
output$below <- unlist(interval.below) 

results <-
  aggregate(x=list(length=output$length,
              low=output$above,
              contained=output$coverage,
              high=output$below),
            by=list(interval=output$interval,
              forest=output$forest,
              n=output$n),
            FUN=mean)

results$scaled <- results$length /
  qt(results$contained, df=as.numeric(as.character(results$n))-1) *
  qt(0.95, df=as.numeric(as.character(results$n))-1) 

results$overCover <- format((qt(0.95,
                                df=as.numeric(as.character(results$n))-1) /
                             qt(results$contained, 
                                df=as.numeric(as.character(results$n))-1))^2,
                            dig=3)

results$length <- format(results$length, dig=3)
results$low <- format(results$low, dig=3)
results$contained <- format(results$contained, dig=3)
results$high <- format(results$high, dig=3)
results$scaled <- format(results$scaled, dig=3)

results[1:12,]

overall.results <-
  aggregate(x=list(length=output$length,
              low=output$above,
              contained=output$coverage,
              high=output$below),
            by=list(interval=output$interval,
              n=output$n),
            FUN=mean)

overall.results$scaled <- overall.results$length /
  qt(overall.results$contained,
     df=as.numeric(as.character(overall.results$n))-1) *
  qt(0.95, df=as.numeric(as.character(overall.results$n))-1) 

overall.results$overCover <-
  format((qt(0.95,
             df=as.numeric(as.character(overall.results$n))-1) /
          qt(overall.results$contained, 
             df=as.numeric(as.character(overall.results$n))-1))^2,
         dig=3)

overall.results$length <- format(overall.results$length, dig=3)
overall.results$low <- format(overall.results$low, dig=3)
overall.results$contained <- format(overall.results$contained, dig=3)
overall.results$high <- format(overall.results$high, dig=3)
overall.results$scaled <- format(overall.results$scaled, dig=3)

overall.results

save.image("booted_srs.RData")

shell("./closure")
