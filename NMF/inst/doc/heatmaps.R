
## ----pkgmaker_preamble, echo=FALSE, results='asis'-----------------------
pkgmaker::latex_preamble()


## ----bibliofile, echo=FALSE, results='asis'------------------------------
pkgmaker::latex_bibliography('NMF')	


## ----options, include=FALSE, verbose=TRUE--------------------------------
#options(prompt=' ')
#options(continue=' ')
set.seed(123456)
library(NMF)


## ----data----------------------------------------------------------------
# random data that follow an 3-rank NMF model (with quite some noise: sd=2)
X <- syntheticNMF(100, 3, 20, noise=2)

# row annotations and covariates
n <- nrow(X)
d <- rnorm(n)
e <- unlist(mapply(rep, c('X', 'Y', 'Z'), 10))
e <- c(e, rep(NA, n-length(e)))
rdata <- data.frame(Var=d, Type=e)

# column annotations and covariates
p <- ncol(X)
a <- sample(c('alpha', 'beta', 'gamma'), p, replace=TRUE)
# define covariates: true groups and some numeric variable
c <- rnorm(p)
# gather them in a data.frame
covariates <- data.frame(a, X$pData, c)


## ----figoptions, include=FALSE-------------------------------------------
library(knitr)
opts_chunk$set(fig.width=14, fig.height=7)


## ----heatmap_data--------------------------------------------------------
par(mfrow=c(1,2))
aheatmap(X, annCol=covariates, annRow=X$fData)
aheatmap(X)


## ----model, cache=TRUE---------------------------------------------------
res <- nmf(X, 3, nrun=10)
res


## ----coefmap_res, fig.keep='last'----------------------------------------
opar <- par(mfrow=c(1,2))
# coefmap from multiple run fit: includes a consensus track
coefmap(res)
# coefmap of a single run fit: no consensus track
coefmap(minfit(res))
par(opar)


## ----coefmap_default, eval=FALSE-----------------------------------------
## Rowv = NA
## Colv = TRUE
## scale = 'c1'
## color = 'YlOrRd:50'
## annCol = predict(object) + predict(object, 'consensus')


## ----coefmap_custom, fig.keep='last', tidy=FALSE-------------------------
opar <- par(mfrow=c(1,2))
# removing all automatic annotation tracks
coefmap(res, tracks=NA)
# customized plot
coefmap(res, Colv = 'euclidean'
	, main = "Metagene contributions in each sample", labCol = NULL
	, annRow = list(Metagene=':basis'), annCol = list(':basis', Class=a, Index=c)
	, annColors = list(Metagene='Set2')
	, info = TRUE)
par(opar)


## ----basismap_res, fig.keep='last'---------------------------------------
opar <- par(mfrow=c(1,2))
# default plot
basismap(res)
# customized plot: only use row special annotation track.  
basismap(res, main="Metagenes", annRow=list(d, e), tracks=c(Metagene=':basis'))
par(opar)


## ----basismap_default, eval=FALSE----------------------------------------
## Colv = NA
## scale = 'r1'
## color = 'YlOrRd:50'
## annRow = predict(object, 'features')


## ----consensusmap_res, fig.keep='last'-----------------------------------
opar <- par(mfrow=c(1,2))
# default plot
consensusmap(res)
# customized plot
consensusmap(res, annCol=covariates, annColors=list(c='blue')
		, labCol='sample ', main='Cluster stability'
		, sub='Consensus matrix and all covariates')
par(opar)


## ----cmap_default, eval=FALSE--------------------------------------------
## distfun = function(x) as.dist(1-x) # x being the consensus matrix
## hclustfun = 'average'
## Rowv = TRUE
## Colv = "Rowv"
## color = '-RdYlBu'


## ----estimate, cache=TRUE------------------------------------------------
res2_7 <- nmf(X, 2:7, nrun=10, .options='v')
class(res2_7)


## ----consensusmap_estimate, fig.keep='last'------------------------------
consensusmap(res2_7)


## ----fit_methods, cache=TRUE---------------------------------------------
res_methods <- nmf(X, 3, list('lee', 'brunet', 'nsNMF'), nrun=10)
class(res_methods)


## ----consensusmap_methods, fig.width=10, fig.height=7, fig.keep='last'----
consensusmap(res_methods)	


## ----demo_hm, eval=FALSE-------------------------------------------------
## demo('aheatmap')
## # or
## demo('heatmaps')


## ----sessionInfo, echo=FALSE, results='asis'-----------------------------
toLatex(sessionInfo())


