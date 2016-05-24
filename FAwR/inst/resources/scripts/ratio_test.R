# $Id: ratio_test.R 172 2006-01-10 05:16:07Z andrewr $

rm(list=ls())

plots <- read.csv("../data/fia_plots.csv")

plots$forest <- factor(plots$forest)

reps <- 2000
reps <- 2

n.s <- c(16, 32, 64, 96, 128)
n.s <- c(16, 32)

ls()

## The following code extracts the plot-level estimates from the FIA database.

#load("/home/andrewr/1.research/FVS\ Refit/fia.wo.VJAN09-06.RData")

#dim(fia.wo)

#names(fia.wo)

#fia.wo[1:10,]

#plot(fia.wo$DBH[1:100], fia.wo$BA[1:100])
#plot(fia.wo$DBH[1:100], fia.wo$BA[1:100]*fia.wo$TF[1:100])

plots <- aggregate(x=list(ba=fia.wo$SBA,
                    ht=fia.wo$HT),
                   by=list(forest=fia.wo$FORCODE,
                     plot=fia.wo$PLOT),
                   FUN=mean, na.rm=T)

#table(plots$forest)

#xyplot(ba ~ ht | forest, data=plots)

#write.csv(plots, file="../data/fia_plots.csv")

require(boot)
require(survey)

boot.ratio <- function(data, index) {
  y <- data$ba
  x <- data$ht
  R <- sum(y[index]) / sum(x[index])
  c(R,
    (var(y[index]) +
     R^2*var(x[index]) -
     2 * R * cov(y[index], x[index])) / length(x))
}

within <- function(conf, x)
  (x >= conf[length(conf)-1] & x < conf[length(conf)])
int.above <- function(conf, x)
  (x >= conf[length(conf)])
int.below <- function(conf, x)
  (x < conf[length(conf)-1])
int.length <- function(conf)
  conf[length(conf)] - conf[length(conf)-1]

## Sample with replacement 

store <- 0

ba.means <- tapply(plots$ba, plots$forest, mean)
ht.means <- tapply(plots$ht, plots$forest, mean)
  
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
      indices <- sample(1:table(plots$forest)[i], size=n, replace=TRUE)
      the.sample <-
        plots[plots$forest==levels(plots$forest)[i],][indices,]
      boot.object <- boot(the.sample, boot.ratio, R=1999)
      ci <- boot.ci(boot.object)
      log.ci <- boot.ci(boot.object, type=c("norm","basic","stud"),
                        h=log, hdot=function(x) 1/x, hinv=exp)
      sqrt.ci <- boot.ci(boot.object, type=c("norm","basic","stud"),
                        h=sqrt, hdot=function(x) 2*x, hinv=function(x) x^2)

      R <- sum(the.sample$ba)/sum(the.sample$ht)
      classical <- R * ht.means[i] + c(-1,1) * qt(p=0.975, df=n-1) * 
                    sqrt((var(the.sample$ba) + R^2*var(the.sample$ht) -
                          2 * R * cov(the.sample$ba, the.sample$ht)) / n)
      
      ratio.lin <- svydesign(id=~1, data=the.sample, weights=1)
      ratio.lin.hat <- svyratio(~ba, ~ht, design=ratio.lin)
      lin.int <-
        matrix(unlist(predict(ratio.lin.hat, total=ht.means[i], se=TRUE)),
               nrow=1, byrow=TRUE) %*%
                 matrix(c(1,1,c(-1,1) * qt(p=0.975, df=n-1)),
                        nrow=2, byrow=TRUE)
          # Yes, I really thought that this would be compact!
      ratio.jack <- as.svrepdesign(ratio.lin)
      ratio.jack.hat <- svyratio(~ba, ~ht, design=ratio.jack)
      jack.int <-
        matrix(unlist(predict(ratio.jack.hat, total=ht.means[i], se=TRUE)),
               nrow=1, byrow=TRUE) %*%
                 matrix(c(1,1,c(-1,1) * qt(p=0.975, df=n-1)),
                        nrow=2, byrow=TRUE)
      
      interval.coverage[[store]] <-
        c(within(classical, ba.means[i]),
          within(lin.int, ba.means[i]),
          within(jack.int, ba.means[i]),
          as.numeric(lapply(ci[4:8], within, x=ba.means[i]/ht.means[i])),
          as.numeric(lapply(log.ci[4:6], within, x=ba.means[i]/ht.means[i])),
          as.numeric(lapply(sqrt.ci[4:6], within, x=ba.means[i]/ht.means[i])))
      interval.lengths[[store]] <- 
        c(int.length(classical),
          int.length(lin.int),
          int.length(jack.int),
          as.numeric(lapply(ci[4:8], int.length)) * ht.means[i],
          as.numeric(lapply(log.ci[4:6], int.length)) * ht.means[i],
          as.numeric(lapply(sqrt.ci[4:6], int.length)) * ht.means[i])
      interval.above[[store]] <- 
        c(int.above(classical, ba.means[i]),
          int.above(lin.int, ba.means[i]),
          int.above(jack.int, ba.means[i]),
          as.numeric(lapply(ci[4:8], int.above, x=ba.means[i]/ht.means[i])),
          as.numeric(lapply(log.ci[4:6], int.above,
                            x=ba.means[i]/ht.means[i])),
          as.numeric(lapply(sqrt.ci[4:6], int.above,
                            x=ba.means[i]/ht.means[i])))
      interval.below[[store]] <- 
        c(int.below(classical, ba.means[i]),
          int.below(lin.int, ba.means[i]),
          int.below(jack.int, ba.means[i]),
          as.numeric(lapply(ci[4:8], int.below, x=ba.means[i]/ht.means[i])),
          as.numeric(lapply(log.ci[4:6], int.below,
                            x=ba.means[i]/ht.means[i])),
          as.numeric(lapply(sqrt.ci[4:6], int.below,
                            x=ba.means[i]/ht.means[i])))
    }
    cat("\nFinished simulation for Forest", levels(plots$forest)[i],
        "for n =", n, "(", reps, "reps )")
  }
  cat("\n")
}

output <- 
  expand.grid(interval=c("classical", "linearized","jackknife",
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

results$overCover <-
  (qt(0.95, df=as.numeric(as.character(results$n))-1) /
   qt(results$contained, df=as.numeric(as.character(results$n))-1))^2

results.tab <- results

results.tab$length <- format(results$length, dig=3)
results.tab$low <- format(results$low, dig=3)
results.tab$contained <- format(results$contained, dig=3)
results.tab$high <- format(results$high, dig=3)
results.tab$scaled <- format(results$scaled, dig=3)
results.tab$overCover <- format(results$overCover, dig=3)

results.tab[1:12,]

overall.results <-
  aggregate(x=list(length=results$length,
              low=results$low,
              contained=results$contained,
              high=results$high),
            by=list(interval=results$interval,
              n=results$n),
            FUN=mean)

overall.results$scaled <- overall.results$length /
  qt(overall.results$contained,
     df=as.numeric(as.character(overall.results$n))-1) *
  qt(0.95, df=as.numeric(as.character(overall.results$n))-1) 

overall.results$overCover <-
  (qt(0.95, df=as.numeric(as.character(overall.results$n))-1) /
   qt(overall.results$contained,
      df=as.numeric(as.character(overall.results$n))-1))^2

overall.results.tab <- overall.results

overall.results.tab$length <- format(overall.results$length, dig=3)
overall.results.tab$low <- format(overall.results$low, dig=3)
overall.results.tab$contained <- format(overall.results$contained, dig=3)
overall.results.tab$high <- format(overall.results$high, dig=3)
overall.results.tab$scaled <- format(overall.results$scaled, dig=3)
overall.results.tab$overCover <- format(overall.results$overCover, dig=3)

overall.results.tab

save(results, overall.results, plots, "booted_ratio.RData")

system("./closure")
