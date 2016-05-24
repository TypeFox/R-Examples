## ----setup, include=FALSE------------------------------------------------
library(knitr)
library(tabplot)

## ----message=FALSE,  results='hide', warning=FALSE-----------------------
require(ggplot2)
data(diamonds)
## add some NA's
is.na(diamonds$price) <- diamonds$cut == "Ideal"
is.na(diamonds$cut) <- (runif(nrow(diamonds)) > 0.8)

## ----message=FALSE, results='hide', fig.height=6.5, fig.width=10, warning=FALSE----
tableplot(diamonds)

## ----results='hide', fig.height=6.5, fig.width=10, warning=FALSE---------
tableplot(diamonds, select = c(carat, price, cut, color, clarity), sortCol = price)

## ----results='hide', fig.height=6.5, fig.width=10, tidy=TRUE, warning=FALSE----
tableplot(diamonds, select = c(carat, price, cut, color, clarity), sortCol = price, from = 0, to = 5)

## ----results='hide', fig.height=6.5, fig.width=10, warning=FALSE---------
tableplot(diamonds, subset = price < 5000 & cut == "Premium")

## ----fig.height=3.5, fig.width=10, warning=FALSE-------------------------
tablePalettes()

## ----results='hide', fig.height=6.5, fig.width=10, warning=FALSE---------
tableplot(diamonds, pals = list(cut="Set1(6)", color="Set5", clarity=rainbow(8)))

## ----message=FALSE,  results='hide', warning=FALSE-----------------------
diamonds$carat_class <- num2fac(diamonds$carat, n=20)
diamonds$price_class <- num2fac(diamonds$price, n=100)

## ----results='hide', fig.height=6.5, fig.width=10, warning=FALSE---------
tableplot(diamonds, select=c(carat, price, carat_class, price_class))

## ----warning=FALSE-------------------------------------------------------
# create large dataset
large_diamonds <- diamonds[rep(seq.int(nrow(diamonds)), 10),]

system.time({
	p <- tablePrepare(large_diamonds)
})

system.time({
	tableplot(p, plot=FALSE)
})

system.time({
	tableplot(p, sortCol=price, nBins=200, plot=FALSE)
})

## ----warning=FALSE-------------------------------------------------------
system.time({
	tableplot(large_diamonds, plot=FALSE)
})

system.time({
	tableplot(large_diamonds, sortCol=price, nBins=200, plot=FALSE)
})

## ----fig.height=6.5, fig.width=10, warning=FALSE-------------------------
system.time({
	tableplot(p, sample=TRUE)
})

## ----fig.height=6.5, fig.width=10, warning=FALSE-------------------------
# calculate normalized carats to be used as sample probabilities
carat.norm <- with(diamonds, carat / max(diamonds$carat))

# draw samples
exp.diamonds <- diamonds[sample(1:nrow(diamonds), size=10000, prob=carat.norm, replace=TRUE),]
chp.diamonds <- diamonds[sample(1:nrow(diamonds), size=10000, prob=1-carat.norm, replace=TRUE),]

tp1 <- tableplot(exp.diamonds, plot=FALSE)
tp2 <- tableplot(chp.diamonds, plot=FALSE)

plot(tp2 - tp1)

## ----results='hide', warning=FALSE---------------------------------------
tab <- tableplot(diamonds, plot = FALSE)

## ----size="scriptsize", fig.keep='none', warning=FALSE-------------------
summary(tab)
plot(tab)

## ----results='hide', fig.height=6.5, fig.width=10, warning=FALSE---------
tableplot(diamonds, select = 1:7, fontsize = 14, legend.lines = 8, title = "Shine on you crazy Diamond", fontsize.title = 18)

## ----fig.height=6.5, fig.width=10, warning=FALSE-------------------------
tab2 <- tableChange(tab, select_string = c("carat", "price", "cut", "color", "clarity"), pals = list(cut="Set1(2)"))
plot(tab2)

## ----warning=FALSE-------------------------------------------------------
tableSave(tab, filename = "diamonds.png", width = 5, height = 3, fontsize = 6, legend.lines = 6)

