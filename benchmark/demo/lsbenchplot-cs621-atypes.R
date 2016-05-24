#' @demo Archetypal analysisof the UCI domain benchmark experiment.
#'
#' @details
#'   Note that domain-based benchmark experiments are not fully
#'   integrated in the workflow provided by benchmark
#'   v0.3-2. Therefore, the one or the other HACK might be available.

library("benchmark")
library("archetypes")

data("uci621raw", package = "benchmark")


alg_cols <- c(lda = "#984EA3", rf = "#FF7F00",
              knn = "#FFFF33", rpart = "#E41A1C",
              svm = "#377EB8", nnet = "#4DAF4A")

alg_light_cols <-c(lda = "#DECBE4", rf = "#FED9A6",
                   knn = "#FFFFCC", rpart = "#FBB4AE",
                   svm = "#B3CDE3", nnet = "#CCEBC5")



### Aggregation of the data: #########################################

dat <- do.call(rbind, lapply(1:6, function(i) uci621raw[, i, 1, ]))
alg <- rep(dimnames(uci621raw)$alg, each = 250)

dat <- na.omit(dat)
alg <- alg[-attr(dat, "na.action")]


### Data set order based on hierarchical clustering:

# data("uci621rel", package = "benchmark")
# d <- relation_dissimilarity(uci621rel)
# hc <- hclust(d, method = "complete")
# dat <- dat[, hc$order]



### Archetypes: ######################################################

set.seed(1234)
as <- stepArchetypes(data = dat, k = 2:10)

screeplot(as)


### Go with k = 4:

k <- 4
a <- bestModel(as[[k-1]])


## Archetypes:
barplot(a, dat)


## Data versus archetypes, parallel coordinates plot:
pcplot(a, dat,
       atypes.col = atype_cols,
       data.col = alg_light_cols[alg])


## Coefficient matrix alpha, parallel coordinates plot:

alpha <- coef(a, "alphas")
colnames(alpha) <- sprintf("A%s", 1:k)

pcplot(alpha, col = alg_cols[alg])

