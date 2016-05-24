## ----, eval=2------------------------------------------------------------
install.packages("rnmf_0.5.tar.gz", repos = NULL, type = "source")
library(rNMF)

## ------------------------------------------------------------------------
data(Symbols_c)

## ----, fig.height=3.2, fig.width=3.2-------------------------------------
see(Symbols_c, title = "Corrupted data set")

## ----, echo=TRUE---------------------------------------------------------
res <- rnmf(Symbols_c, k = 4, showprogress = FALSE, my.seed = 100)

## ----, fig.height=3.2, fig.width=3.2-------------------------------------
see(res$fit, title = "Regular NMF reconstruction with k = 4")

## ----, fig.height=1, fig.width=3, warning=FALSE--------------------------
see(res$W, title = "Regular NMF basis", layout = c(1,4))

## ----, echo=TRUE---------------------------------------------------------
res2 <- rnmf(Symbols_c, k = 4, gamma = 0.03, showprogress = FALSE, 
             my.seed = 100, tol = 0.0001, maxit = 50)

## ----, fig.height=3.2, fig.width=3.2-------------------------------------
see(res2$fit, title = "rNMF reconstruction with k = 4")

## ----, fig.height=1, fig.width=3, warning=FALSE--------------------------
see(res2$W, title = "rNMF basis vectors", layout = c(1,4))

## ----, fig.height=3.5, fig.width=3.5-------------------------------------
outliers <- matrix(0, nrow = nrow(Symbols_c), ncol = ncol(Symbols_c))
outliers[res2$trimmed[[res2$niter]]] <- 1
see(outliers, title = "Outliers extracted by rNMF")

## ------------------------------------------------------------------------
data(face)

## ----, fig.height=5/1.1, fig.width=4/1.1---------------------------------
see(face, title = "Corrupted face image", col = "grey", input = "single")

## ----, echo=TRUE---------------------------------------------------------
res <- rnmf(face, k = 10, showprogress = FALSE, my.seed = 100)

## ----, fig.height=5/1.1, fig.width=4/1.1---------------------------------
see(res$fit, title = "NMF compression", col = "grey", input = "single")

## ----, echo=TRUE---------------------------------------------------------
res2 <- rnmf(face, k = 10, gamma = 0.025, showprogress = FALSE, my.seed = 100)

## ----, fig.height=5/1.1, fig.width=4/1.1---------------------------------
see(res2$fit, title = "rNMF compression", col = "grey", input = "single")

