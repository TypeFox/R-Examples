#' @demo Exploratory and inferential analysis of the UCI domain
#' benchmark experiment.
#'
#' @details
#'   Analysis used in 'Analysis of Domain-based Benchmark Experiments'
#'   by Manuel J. A. Eugster, Torsten Hothorn and Friedrich Leisch.
#'
#'   Note that domain-based benchmark experiments are not fully
#'   integrated in the workflow provided by benchmark
#'   v0.3-2. Therefore, the one or the other HACK might be available.

library("benchmark")

data("uci621raw", package = "benchmark")
data("uci621rel", package = "benchmark")


alg_cols <- c(lda = "#984EA3", rf = "#FF7F00",
              knn = "#FFFF33", rpart = "#E41A1C",
              svm = "#377EB8", nnet = "#4DAF4A")

alg_light_cols <-c(lda = "#DECBE4", rf = "#FED9A6",
                   knn = "#FFFFCC", rpart = "#FBB4AE",
                   svm = "#B3CDE3", nnet = "#CCEBC5")

alg_light_cols2 <- c(none = "white", alg_light_cols)



### Visualizing the performance distributions ########################

dat1 <- melt(uci621raw[, , "Misclassification", ])
dat1 <- na.omit(dat1)
dat1$samp <- factor(dat1$samp)


### Trellis plot:

ggplot(dat1, aes(alg, value, fill = alg, colour = alg)) +
    facet_wrap(~ ds) + geom_boxplot(outlier.size = 1) +
    scale_fill_manual(value = alg_light_cols) +
    scale_colour_manual(value = alg_cols)



### Dendrogram:

d <- relation_dissimilarity(uci621rel)

hc <- hclust(d, method = "complete")

plot(hc, main = NULL, sub = "", xlab = "")




### Stacked bar plot:

m <- apply(uci621raw[, , "Misclassification", ,drop=FALSE], 4,
           function(x) {
             apply(x, 2, mean, na.rm = TRUE)
           })

barplot(m[, hc$order], horiz = TRUE, border = 1, col = alg_light_cols)



### Summary plot:

bsplot0(uci621rel, col = alg_light_cols, sig.lwd = 4,
        stat = m, stat.col = alg_cols,
        ds.order = names(uci621rel)[hc$order],
        alg.order = dimnames(uci621raw)[[2]][c(5, 1:4, 6)],
        ylab.las = 2, ylab = "")


### Summary graph:

library("Rgraphviz")
library("vcd")

rm <- benchmark:::bsranking(uci621rel)
uw <- apply(rm, 2,
            function(x) {
              w <- which(x == 1)
              ifelse(length(w) == 1,
                     names(w), "none")
            })


n <- 10

bsgraph0(d, ndists.show = n,
         edge.col = terrain_hcl(n, c=c(65,0), l=c(45,90), power=c(1/2,1.5)),
         edge.lwd = seq(4, 1, length.out=n),
         node.fill = alg_light_cols[uw])

symbols(180, 1060, circles = 185, inches = FALSE, add = TRUE,
        lwd = 1, lty = 2)
text(60, 850, "A")

abline(a = 500, b = 0.5, lwd = 1, lty = 2)
text(60, 600, "C")
text(60, 470, "B")



### Mixed effects model: #############################################

library("lme4")
library("multcomp")


dat1 <- uci621raw[, , "Misclassification", ]
na <- unique(which(is.na(dat1), arr.ind = TRUE)[,1])
dat1 <- melt(dat1[-na, , ])
dat1$samp <- factor(dat1$samp)
dat1$alg <- relevel(dat1$alg, "lda")
dat1$dssamp <- dat1$ds:dat1$samp


### Model:

mod1 <- lmer(value ~ alg - 1 + (1 | ds) + (alg - 1 | ds) + (1 | dssamp), data = dat1)
print(mod1, corr = FALSE)


## Random effects:
ranef(mod1)$ds
round(ranef(mod1)$ds[which.min(ranef(mod1)$ds[, "(Intercept)"]), ], 4)
round(ranef(mod1)$ds[which.max(ranef(mod1)$ds[, "(Intercept)"]), ], 4)


## Fixed effects:
round(fixef(mod1), 4)



### Pairwise tests:

pt <- glht(mod1, linfct = mcp(alg = "Tukey"))
plot(pt, main = "")



### Pairwise tests as relation (HACK):

as.PaircompDecision <- function(test, model, algorithms,
                                significance = 0.05,
                                relevance = 0.10) {

  result <- benchmark:::emptyLeDecision(algorithms)

  pt <- benchmark:::LmerPairwiseTestResult$new(test)
  ci <- pt$getConfint(1 - significance)

  desc <- !(ci[, 'lwr'] < 0 & ci[, 'upr'] > 0)
  desc <- desc & !(ci[, 'lwr'] > -relevance & ci[, 'upr'] < relevance)

  sigdirs <- sign(ci[desc, 'Estimate'])
  sigpairs <- strsplit(rownames(ci)[desc], ' - ')
  sigpairs[sigdirs == 1] <- lapply(sigpairs[sigdirs == 1], rev)

  for ( p in sigpairs )
    result[p[1], p[2]] <- 1

  benchmark:::PaircompDecision(result, "<",
                               list(model = model, globaltest = NULL,
                                    pairwisetest = pt, confint = ci))
}

pd <- as.PaircompDecision(pt, mod1, names(alg_cols))

r <- relation(incidence = pd$decision)

plot(r, main = "", attrs = attrs)
plot(transitive_reduction((r & dual(r))))




### Consensus of local relations: ####################################

lapply(uci621rel, as.ranking)

plot(relation_consensus(uci621rel, method = "SD/L", all = TRUE))
plot(relation_consensus(uci621rel, method = "SD/O", all = TRUE))

