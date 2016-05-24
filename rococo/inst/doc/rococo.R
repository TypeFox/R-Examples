## ----Init,echo=FALSE,message=FALSE,results='hide'---------------------
options(width=72)
knitr::opts_knit$set(width=72)
set.seed(0)
library(rococo, quietly=TRUE)
rococoVersion <- packageDescription("rococo")$Version
rococoDateRaw <- packageDescription("rococo")$Date
rococoDateYear <- as.numeric(substr(rococoDateRaw, 1, 4))
rococoDateMonth <- as.numeric(substr(rococoDateRaw, 6, 7))
rococoDateDay <- as.numeric(substr(rococoDateRaw, 9, 10))
rococoDate <- paste(month.name[rococoDateMonth], " ",
                     rococoDateDay, ", ",
                     rococoDateYear, sep="")

## ----InstallRoCoCo,eval=FALSE-----------------------------------------
#  install.packages("rococo")

## ----LoadRoCoCo,eval=FALSE--------------------------------------------
#  library(rococo)

## ----OpenVignette,eval=FALSE------------------------------------------
#  vignette("rococo")

## ----ShowHelp,eval=FALSE----------------------------------------------
#  help(rococo)

## ----ToyData,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
x1 <- rnorm(15)
y1 <- 2 * x1 + rnorm(length(x1), sd=0.25)
plot(x1, y1)

## ----SimpleRococo-----------------------------------------------------
rococo(x1, y1)

## ----SimpleRococoTest-------------------------------------------------
rococo.test(x1, y1, alternative="two.sided")

## ----RoCoCoTestFormula,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
data(iris)
plot(~ Sepal.Length + Petal.Length, iris)
rococo.test(~ Sepal.Length + Petal.Length, iris, alternative="two.sided")

## ----FlatNoisyData,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
x2 <- rnorm(15)
f2 <- function(x) ifelse(x > 0.8, x - 0.8, ifelse(x < -0.8, x + 0.8, 0))
y2 <- f2(x2) + rnorm(length(x2), sd=0.1)
plot(x2, y2)

## ----FlatNoisyDefault-------------------------------------------------
rococo.test(x2, y2, alternative="greater")

## ----FlatNoisyOther---------------------------------------------------
rococo.test(x2, y2, similarity="gauss",
            alternative="greater")
rococo.test(x2, y2, similarity=c("classical", "gauss"),
            alternative="greater")

## ----FlatNoisySameR---------------------------------------------------
rococo.test(x2, y2, similarity="gauss", r=0.1, alternative="greater")

## ----FlatNoisyDifferentR----------------------------------------------
rococo.test(x2, y2, similarity=c("linear", "gauss"), r=c(0.05, 0.1),
            alternative="greater")

## ----CheckIQR---------------------------------------------------------
rococo.test(x2, y2, similarity=c("linear", "gauss"), r=0,
            alternative="greater")
IQR(x2) * 0.1
IQR(y2) * 0.1

## ----CheckIQR2--------------------------------------------------------
rococo.test(x2, y2, similarity=c("linear", "gauss"), tnorm="prod",
            alternative="greater")

## ----YagertNorm-------------------------------------------------------
DrastictNorm <- function(x, y)
{
    if (x == 1) y
    else if (y == 1) x
    else 0
}
YagertNorm <- function(lambda)
{
    fun <- function(x, y)
    {
        if (lambda == 0)
            DrastictNorm(x, y)
        else if (is.infinite(lambda))
            min(x, y)
        else
            max(0, 1 - ((1 - x)^lambda + (1 - y)^lambda)^(1 / lambda))
    }

    attr(fun, "name") <- paste("Yager t-norm with lambda =", lambda)

    fun
}
rococo(x2, y2, tnorm=YagertNorm(0.5))
rococo.test(x2, y2, tnorm=YagertNorm(0.2))

## ----LargeNoShuffles,echo=TRUE----------------------------------------
res <- rococo.test(x2, y2, numtests=100000, storeValues=TRUE)
res

## ----ExactpValueExample-----------------------------------------------
rococo.test(x2[1:8], y2[1:8], exact=TRUE)

## ----LargeNoShufflesPic,fig.width=8,fig.height=5.5,out.width='0.7\\textwidth'----
hist(res@perm.gamma, breaks=100, probability=TRUE, xlab="gamma",
     main="Distribution of gamma for random shuffles")
plot(function(x) dnorm(x, mean=res@H0gamma.mu, sd=res@H0gamma.sd),
     min(res@perm.gamma), max(res@perm.gamma), col="red", lwd=2, add=TRUE)

## ----LargeNoShufflespVal,echo=TRUE------------------------------------
res@p.value.approx

## ----SimpleGaussCorTest-----------------------------------------------
gauss.cor(x1, y1)
gauss.cor.test(x1, y1, alternative="two.sided")
gauss.cor.test(~ Sepal.Length + Petal.Length, iris, alternative="two.sided")

## ----GetBibTeX,eval=FALSE---------------------------------------------
#  toBibtex(citation("rococo"))

