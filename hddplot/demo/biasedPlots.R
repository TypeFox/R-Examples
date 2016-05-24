## High-dimensional data, classification, and plots

##                    What groups are of interest?
data(Golub)
data(golubInfo)      # 7129 rows by 72 columns
with(golubInfo, table(cancer, tissue.mf))
attach(golubInfo)
## Identify allB samples for that are BM:f or BM:m or PB:m
subsetB <- cancer=="allB" & tissue.mf%in%c("BM:f","BM:m","PB:m")
## Form vector that identifies these as BM:f or BM:m or PB:m
tissue.mfB <- tissue.mf[subsetB, drop=TRUE]
## Separate off the relevant columns of the matrix Golub
GolubB <- Golub[, subsetB]
detach(golubInfo)

## Select the 15 features that are the "best discriminators",
## as judged by the aov F-statistic
## Select relevant rows & transpose to observations by features
ord15 <- orderFeatures(GolubB, cl=tissue.mfB)[1:15]
dfB.ord <- data.frame(t(GolubB[ord15, ]))
## Calculations for the left panel
dfB15 <- data.frame(t(GolubB[ord15, ]))
library(MASS)
dfB15.lda <-  lda(dfB15, grouping=tissue.mfB)
scores <- predict(dfB15.lda, dimen=2)$x
## Scores for the single PB:f observation
attach(golubInfo)
df.PBf <- data.frame(t(Golub[ord15, tissue.mf=="PB:f"
                                    & cancer=="allB", drop=FALSE]))
scores.PBf <- predict(dfB15.lda, newdata=df.PBf, dimen=2)$x
detach(golubInfo)
## Now create two plots
oldpar <- par(mfrow=c(1,2), pty="s")
## Plot that is based on Golub data
scoreplot(scorelist=list(scores=scores, cl=tissue.mfB,
                      other=scores.PBf, cl.other=factor("PB:f"),
                      nfeatures=15),
                      params=list(other=list(cex=1.2,col=4,pch=4)))
## Plot that is based on random data
simscores <- simulateScores(nrows=7129, cl=rep(1:3, c(19,10,2)),
                             cl.other=4, nfeatures=15, seed=41)
  # Returns list elements: scores, cl, scores.other & cl.other
scoreplot(simscores, params=list(other=list(cex=0.9,col=4,pch=4)),
          prefix.title="Random data:")

## Alternatively, set:
## dimen <- dim(GolubB); rsetB <- array(rnorm(prod(dimen)), dim=dimen)
## Then replace GolubB with rsetB, and repeat the lines above.
par(oldpar)
