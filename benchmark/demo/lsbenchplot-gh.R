#' @demo Exploratory and inferential analysis of the Grasshopper
#'   domain benchmark experiment.
#'
#' @details
#'   Analysis used in 'Analysis of Domain-based Benchmark Experiments'
#'   by Manuel J. A. Eugster, Torsten Hothorn and Friedrich Leisch.
#'
#'   Note that domain-based benchmark experiments are not fully
#'   integrated in the workflow provided by benchmark
#'   v0.3-2. Therefore, the one or the other HACK might be available.

library("benchmark")

data("ghraw", package = "benchmark")
data("ghrel", package = "benchmark")

alg_cols <- c(lda = "#984EA3", rf = "#FF7F00",
              knn = "#FFFF33", rpart = "#E41A1C",
              svm = "#377EB8", nb = "#A65628")

alg_light_cols <-c(lda = "#DECBE4", rf = "#FED9A6",
                   knn = "#FFFFCC", rpart = "#FBB4AE",
                   svm = "#B3CDE3", nb = "#E5D8BD")

alg_light_cols2 <- c(none = "white", alg_light_cols)



### Consensus: #######################################################

r <- relation_consensus(ghrel, method = "SD/L")
plot(r)


### Visualizations: ##################################################

### Trellis plot:

ggplot(subset(ghraw, perf == "vmiscl"), aes(alg, value, fill = alg, colour = alg)) +
  geom_boxplot(outlier.size = 1) + facet_wrap(~ ds) +
  ylab("Misclassification") + xlab("Algorithm") +
  scale_fill_manual(value = alg_light_cols) +
  scale_colour_manual(value = alg_cols)



### Benchmark summary graph:

library("Rgraphviz")
library("vcd")

med <- ddply(subset(ghraw, perf == "vmiscl"), c("ds"),
             function(x) sapply(split(x$value, x$alg), median))
medalg <- apply(as.matrix(med[, -1]), 1, which.min)


d <- relation_dissimilarity(ghrel)

pdf("lsbenchplot-bsgraph.pdf", width = 5, height = 4, pointsize = 8)
n <- 6
bsgraph0(d, ndists.show = n,
         edge.col = terrain_hcl(n, c=c(65,0), l=c(45,90), power=c(1/2,1.5)),
         edge.lwd = seq(4, 1, length.out=n),
         node.fill = alg_light_cols[colnames(med)[-1][medalg]])
dev.off()



### Linear mixed-effects model: ######################################

library("lme4")
library("multcomp")

ghraw$dssamp <- ghraw$ds:ghraw$samp


### Model:

mod1 <- lmer(value ~ alg - 1 + (1 | ds) + (alg - 1 | ds) + (1 | dssamp),
             data = subset(ghraw, perf == "vmiscl"))

print(mod1, corr = FALSE)
ranef(mod1)$ds



### Pairwise tests:

pt <- glht(mod1, linfct = mcp(alg = "Tukey"))
plot(pt)



### Pairwise tests as relation (HACK):

as.PaircompDecision <- function(test, model, algorithms,
                                significance = 0.05,
                                relevance = 0) {

  result <- emptyLeDecision(algorithms)

  pt <- LmerPairwiseTestResult$new(test)
  ci <- pt$getConfint(1 - significance)

  desc <- !(ci[, 'lwr'] < 0 & ci[, 'upr'] > 0)
  desc <- desc & !(ci[, 'lwr'] > -relevance & ci[, 'upr'] < relevance)

  sigdirs <- sign(ci[desc, 'Estimate'])
  sigpairs <- strsplit(rownames(ci)[desc], ' - ')
  sigpairs[sigdirs == 1] <- lapply(sigpairs[sigdirs == 1], rev)

  for ( p in sigpairs )
    result[p[1], p[2]] <- 1

  PaircompDecision(result, "<",
                   list(model = model, globaltest = NULL,
                        pairwisetest = pt, confint = ci))
}


pd <- as.PaircompDecision(pt, mod1, names(alg_cols)[-6])
r <- relation(incidence = pd$decision)

plot(r)
