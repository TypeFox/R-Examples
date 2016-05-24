## ----setup, include=FALSE---------------------------------
library("knitr")
### Set global chunk options
opts_chunk$set(eval=TRUE,
   ## text results
   echo=TRUE,
   results=c('markup', 'asis', 'hold', 'hide')[1],
   collapse=FALSE,
   warning=TRUE, message=TRUE, error=TRUE,
   split=FALSE, include=TRUE, strip.white=TRUE,
   ## code decoration
   tidy=FALSE, prompt=FALSE, comment='##',
   highlight=TRUE, size='normalsize',
   background=c('#F7F7F7', colors()[479], c(0.1, 0.2, 0.3))[1],
   ## cache
   cache=FALSE,
   ## plots
   fig.path=c('figure', 'figure/minimal-')[1],
   fig.keep=c('high', 'none', 'all', 'first', 'last')[1],
   fig.align=c('center', 'left', 'right', 'default')[1],
   fig.show=c('hold', 'asis', 'animate', 'hide')[1],
   dev=c('pdf', 'png', 'tikz')[1],
   fig.width=7, fig.height=7, #inches
   fig.env=c('figure', 'marginfigure')[1],
   fig.pos=c('', 'h', 't', 'b', 'p', 'H')[1])
### Set R options
options(formatR.arrow=TRUE, width=60)

## ----prob-------------------------------------------------
p <- seq(5)
h <- c(1, 3, 1.5, 3, 1)
plot(p, h, type="b", 
     col="blue", axes=FALSE,
     xlab="probability",
     ylab="size of statistic")
axis(1, at=p, labels=c("<0.1", "0.1-0.3", "0.3-0.7", "0.7-0.9", ">0.9"))
axis(2, at=c(1, 1.5, 2, 3), labels=c("low", "low-med", "med", "high"))
dChisq <- c(2, 1.5, 2)
points(c(2, 3, 4), dChisq, type="b", col="green")
lines(x=c(1, 2), y=c(3, 2), col="green")
lines(x=c(1, 2), y=c(1, 2),col="green")
lines(x=c(4, 5), y=c(2, 1), col="green")
lines(x=c(4, 5), y=c(2, 3),col="green")
dBhat <- c(1, 3, 2, 3, 1)
points(p, dBhat, type="b", col="red")
legend(2.5, y=1.4, legend=c("h", "dChisq", "dBhat"), 
       fill=c("blue", "green", "red"))
mtext("Probability vs. size of statistic")

## ----plots------------------------------------------------
library("LogisticDx")
## H&L 2nd ed. Table 4.9. Figures 5.5-5.8. Pages 177-180.
data(uis)
uis <- within(uis, {
    NDRGFP1 <- 10 / (NDRGTX + 1)
    NDRGFP2 <- NDRGFP1 * log((NDRGFP1 + 1) / 10)
})
summary(g1 <- glm(DFREE ~ AGE + NDRGFP1 + NDRGFP2 + IVHX +
                  RACE + TREAT + SITE +
                  AGE:NDRGFP1 + RACE:SITE,
                  family=binomial, data=uis))
plot(g1, devNew=FALSE)

