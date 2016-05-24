#' @demo Exploratory and inferential analysis of the monks3 benchmark
#'   experiment.
#'
#' @details
#'   Analysis used in 'Exploratory and Inferential Analysis of
#'   Benchmark Experiments' by Manuel J. A. Eugster, Torsten Hothorn
#'   and Friedrich Leisch.

library("benchmark")

data("monks3raw", package = "benchmark")

w <- as.warehouse.array4dim(monks3raw)



### Common analysis of benchmark experiments: ########################

apm <- w$viewAlgorithmPerformance(performances = "Misclassification")


### Performance point estimations:
ddply(apm, "algorithms",
      function(x)
      c(mean = mean(x$value), sd = sd(x$value),
        median = median(x$value), max = max(x$value)))


### Bootstrap percentile intervals:
pci <- PercintTest$new(apm)$pairwiseTest()$percint(0.05)

benchmark:::plot.percint(pci)


### Average-case scenario ranking:
rp1 <- as.relation(paircomp(apm, family = GenericPointPaircomp,
                            type = "<", estimator = "mean"))
rp1


### Worst-case scenario ranking:
rp2 <- as.relation(paircomp(apm, family = GenericPointPaircomp,
                            type = "<", estimator = "max"))
rp2



### Exploratory analysis: ############################################

### Basic plots:
stripchart(apm)
boxplot(apm)
densityplot(apm)


### Benchmark experiment plot:
set.seed(38)
beplot0(apm)

set.seed(38)
beplot0(apm, lines.show = TRUE, lines.col = rep(gray(0.7), 6))



### Inference: #######################################################

### Non-parametric -- Friedman test based:
p1 <- paircomp(apm, family = FriedmanTestPaircomp, type = "<",
               significance = 0.05)
p1


## Decision base:
p1$base$globaltest
statistic(p1$base$globaltest$test, type = "linear")

p1$base$pairwisetest
statistic(p1$base$pairwisetest$test, type = "linear")
pvalue(p1$base$pairwisetest$test, method = "single-step")


### Parametric -- Mixed effects model based:
p2 <- paircomp(apm, family = LmerTestPaircomp, type = "<",
               significance = 0.05)
p2


## Decision base:
p2$base$model
fixef(p2$base$model)
VarCorr(p2$base$model)

p2$base$globaltest

p2$base$pairwisetest
plot(p2$base$pairwisetest$test)


## Model diagnostic:
m <- p2$base$model

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
boxplot(fitted(m) ~ m@frame$algorithms, ylab = "Fitted values")
boxplot(m@frame$value ~ m@frame$algorithms, ylab = "True values")
par(op)

boxplot(resid(m) ~ m@frame$algorithms, ylab = "Residuals")


## Significance versus relevance:
p2b <- paircomp(apm, family = LmerTestPaircomp, type = "<",
                significance = 0.05, relevance = 0.01)
p2b



### Preference relations (with package "relations"): #################

### Point estimation based preferences see above:
rp1
rp2


### Statistical test based preferences:
rt1 <- as.relation(p1)
rt1
plot(rt1)

rt2 <- as.relation(p2)
rt2
plot(rt2)

rt2b <- as.relation(p2b)
rt2b
plot(rt2b)



### Preference combination (with package "relations"): ###############

### Computation time:
apm <- w$viewAlgorithmPerformance(performances = "Time")
stripchart(apm)

rp3 <- as.relation(paircomp(apm, family = GenericPointPaircomp,
                            estimator = "mean"))
rp3



### Consensus:
R <- relation_ensemble(Rm = rt2, Rw = rp2, Rc = rp3)

as.ranking(relation_consensus(R, method = "SD/O"))

as.ranking(relation_consensus(R, method = "SD/P"),
           weights = c(1, 1.5, 0.2))

