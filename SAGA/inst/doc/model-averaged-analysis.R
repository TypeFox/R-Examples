## ---- echo=FALSE, results='asis'-----------------------------------------
foo <- read.csv(system.file("cmatrix.xy.csv", package ="SAGA"))[1:15,1:13]
knitr::kable(foo, row.names=T, output=T)

## ------------------------------------------------------------------------
NEW.cmat <- SAGA::DisplayCmatrix(table = "XY")

## ------------------------------------------------------------------------
colnames(NEW.cmat)

## ------------------------------------------------------------------------
NEW.cmat <- NEW.cmat[, c(1:3,6:12,20:21)]

## ---- echo=FALSE, results='asis'-----------------------------------------
data(per.inf, package="SAGA")
colnames(per.inf) <- c("Cohort ID", "Mean", "SE") 
knitr::kable(per.inf, row.names=F, output=T)

## ---- echo=FALSE---------------------------------------------------------
cohorts <- SAGA::cohortID()
prows <- nrow(cohorts)/2
cohorts <- cbind(cohorts[1:prows,],rep(" ",prows),rep(" ",prows),cohorts[(prows+1):(prows*2),])
colnames(cohorts) <- c("ID","cohorts","","","ID","cohorts")
knitr::kable(cohorts, output=T)

## ----fig.cap='Model averaged estimate of genetic architecture.', echo=TRUE, warning=FALSE,results='markup', fig.width=6.5----
# we will need the plotrix package for plotting
library(plotrix)
results <- SAGA::AnalyzeCrossesMM(per.inf, graph=T, cex.names=.8)

## ---- fig.cap='Unconditional estimate of genetic architecture.', echo=TRUE, warning=FALSE, results='markup', fig.width=7----
library(plotrix)
#Sperm receptacle length in Drosophila mojavensis
data(SR, package="SAGA")
#Because we are using cohorts where we know the distribution of sexes we set sexed=T.
results2 <- SAGA::AnalyzeCrossesMM(SR, even.sex=T, graph=T, cex.names=.7)

## ------------------------------------------------------------------------
names(results)

## ---- eval=F-------------------------------------------------------------
#  results[[2]]

## ---- echo=F-------------------------------------------------------------
results[[2]][1:nrow(results[[2]]),1:ncol(results[[2]])] <- round(as.numeric(results[[2]]), 3)
results[[2]]

## ------------------------------------------------------------------------
results[[4]]

## ----fig.cap='Subset of model averaged estimate of genetic architecture.', echo=TRUE, warning=FALSE, results='hide', fig.width=5, fig.height=3----
    # here we extract the 4 largest composite effects found in the first analysis
    estimates <- as.numeric(results[[2]][1, c(3, 7, 8, 9)])
    names(estimates) <- colnames(results[[2]])[c(3, 7, 8, 9)]
    barplot(estimates, main = "Estimate for composite effects",
            names.arg = names(estimates))

## ----fig.cap='The standard full plot returned from the analysis', fig.width = 6, fig.height = 2.5, echo=F----
# first the standard full plot returned from the analysis
SAGA::plot.genarch(results, min.vi=0)

## ----fig.cap='Now we have reduced the plot to include just those CGEs with a *vi* of at least .25', fig.width = 6, fig.height = 2.5, echo=F----
SAGA::plot.genarch(results, min.vi=.25)

## ---- fig.cap='Now have switched to the viridis color pallete to make it color blind friendly, adjusted the Y axis, and changed the main title.', fig.width = 5.5, fig.height = 2.5, echo=F----
SAGA::plot.genarch(results, maxval=85, min.vi=.253, main="A nicer plot", viridis=T)

## ------------------------------------------------------------------------
# first lets find the best two models
order(results[[3]])[1:2]

## ---- warning=FALSE, results='hide', fig.width = 3.5, fig.height = 2.5, fig.cap='Estimates conditional on model 166'----
SAGA::EvaluateModel(results, 166, cex.names=.7, cex.main=.7)

## ---- warning=FALSE, results='hide', fig.width = 3.5, fig.height = 2.5, fig.cap='Estimates conditional on model 110'----
SAGA::EvaluateModel(results, 110, cex.names=.7, cex.main=.7)

## ---- fig.cap='Observed line means.', fig.width = 4, fig.height = 4, echo=F----
    SAGA::plotObserved(SR, results)

## ---- fig.cap='Distribution of akaike weights across model space for _Tribolium_ dataset.', echo=FALSE, warning=FALSE, results='hide', fig.width=4.5, fig.height=4----
SAGA::VisModelSpace(results, cex.u=1.6)

## ---- fig.cap='Distribution of akaike weights across model space for _Drosophila_ sperm receptacle length.', echo=FALSE, warning=FALSE, results='hide', fig.width=4.5, fig.height=4----
SAGA::VisModelSpace(results2, cex.u=.4)

