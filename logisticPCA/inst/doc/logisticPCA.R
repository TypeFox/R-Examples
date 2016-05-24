## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=5)

## ----setup---------------------------------------------------------------
library(logisticPCA)
library(ggplot2)
data("house_votes84")

## ----table, echo=FALSE---------------------------------------------------
df = data.frame(
  Formulation = c("Exponential Family PCA", "Logistic PCA", "Convex Logistic PCA"),
  Function = c("logisticSVD", "logisticPCA", "convexLogisticPCA"),
  Class = c("`lsvd`", "`lpca`", "`clpca`"),
  Returns = c("`mu`, `A`, `B`", "`mu`, `U`, (`PCs`)", "`mu`, `H`, (`U`, `PCs`)"),
  Specify_M = c("No", "Yes", "Yes")
)
knitr::kable(df, col.names = c("Formulation", "Function", "Class", "Returns", "Specify m?"))

## ----lsvd----------------------------------------------------------------
logsvd_model = logisticSVD(house_votes84, k = 2)

## ----printlsvd-----------------------------------------------------------
logsvd_model

## ----cvlpca--------------------------------------------------------------
logpca_cv = cv.lpca(house_votes84, ks = 2, ms = 1:10)
plot(logpca_cv)

## ----lpca----------------------------------------------------------------
logpca_model = logisticPCA(house_votes84, k = 2, m = which.min(logpca_cv))
clogpca_model = convexLogisticPCA(house_votes84, k = 2, m = which.min(logpca_cv))

## ----clpca_trace---------------------------------------------------------
plot(clogpca_model, type = "trace")

## ----lsvd_trace----------------------------------------------------------
plot(logsvd_model, type = "trace")

## ----plot, warning=FALSE-------------------------------------------------
party = rownames(house_votes84)
plot(logsvd_model, type = "scores") + geom_point(aes(colour = party)) + 
  ggtitle("Exponential Family PCA") + scale_colour_manual(values = c("blue", "red"))
plot(logpca_model, type = "scores") + geom_point(aes(colour = party)) + 
  ggtitle("Logistic PCA") + scale_colour_manual(values = c("blue", "red"))
plot(clogpca_model, type = "scores") + geom_point(aes(colour = party)) + 
  ggtitle("Convex Logistic PCA") + scale_colour_manual(values = c("blue", "red"))

## ----fitted--------------------------------------------------------------
head(fitted(logpca_model, type = "response"))

## ----fake----------------------------------------------------------------
d = ncol(house_votes84)
votes_fake = matrix(sample(c(0, 1), 5 * d, replace = TRUE), 5, d,
                    dimnames = list(NULL, colnames(house_votes84)))


## ----predict-------------------------------------------------------------
predict(logpca_model, votes_fake, type = "PCs")

