## ----include=FALSE,message=FALSE-----------------------------------------
require(fastR)
require(mosaicData)
require(knitr)
opts_chunk$set(
			   fig.show="hold",
			   fig.width=7, fig.height=4,
			   out.width=".47\\textwidth"
			   )

## ------------------------------------------------------------------------
require(mosaicData)

## ----eval=FALSE----------------------------------------------------------
#  goal( formula, data = mydata, ...)

## ----eval=FALSE----------------------------------------------------------
#  goal( ~ x,  data = mydata)       # for single variable summaries
#  goal( y ~ x, data = mydata)      # for two-variable summaries and linear models
#  goal( y ~ x | z, data = mydata)  # for multi-variable summaries and faceting in plots

## ------------------------------------------------------------------------
require(fastR)
trellis.par.set(theme = col.mosaic())   # change default colors, etc.
table( iris $ Species )
tally( ~ Species, data = iris )

## ------------------------------------------------------------------------
tally( ~ Species, data = iris, margins = FALSE )

## ------------------------------------------------------------------------
tally( ~ Species, data = iris, format = "count")
tally( ~ Species, data = iris, format = "percent")
tally( ~ Species, data = iris, format = "proportion")

## ------------------------------------------------------------------------
mean( ~ Sepal.Length, data = iris)
median( ~ Sepal.Length, data = iris)
sd( ~ Sepal.Length, data = iris)
iqr( ~ Sepal.Length, data = iris)
favstats( ~ Sepal.Length, data = iris)

## ------------------------------------------------------------------------
mean( Sepal.Length ~ Species, data = iris )
favstats( Sepal.Length ~ Species, data = iris )

## ----eval=FALSE, tidy=FALSE----------------------------------------------
#  ?mean

## ------------------------------------------------------------------------
qdata( ~ Sepal.Length, p = 0.5, data = iris )
median( ~ Sepal.Length, data = iris )
pdata( ~ Sepal.Length, 5,  data = iris )
tally( ~ (Sepal.Length <= 5), data = iris, format = "proportion")

## ------------------------------------------------------------------------
bargraph( ~ substance, data = HELPrct)
bargraph( ~ substance, data = HELPrct, groups = sex )

## ------------------------------------------------------------------------
histogram( ~ week1, data = fumbles, width = 1 )

## ------------------------------------------------------------------------
histogram( ~ Sepal.Length, data = iris, groups = Sepal.Length > 5, h = c(.1,.2) )
histogram( ~ Sepal.Length | Species, data = iris, fit = "normal", v = 6 )

## ----eval=FALSE----------------------------------------------------------
#  mPlot(iris)
#  mPlot(HELPrct, "density")

## ------------------------------------------------------------------------
rflip(10)
do(3) * rflip(10)
Flips <- do(1000) * rflip(10)
tally( ~ heads, data = Flips)
histogram( ~ heads, data = Flips, width = 1)

## ------------------------------------------------------------------------
plotDist("binom", params = list(size = 10,prob = .5))
plotDist("binom", params = list(size = 10,prob = .5), kind = 'cdf')

## ------------------------------------------------------------------------
plotDist("binom", params = list(size = 10,prob = .5), kind = 'hist')
plotDist("binom", params = list(size = 10,prob = .5), kind = 'qq')

## ------------------------------------------------------------------------
plotDist("chisq", params = list(df = 4)) 
plotDist("chisq", params = list(df = 4), kind = 'cdf')

## ------------------------------------------------------------------------
binom.test( ~ sex, data = HELPrct )

## ------------------------------------------------------------------------
pval( binom.test( ~ sex, data = HELPrct ) )
confint( binom.test( ~ sex, data = HELPrct ) )

## ------------------------------------------------------------------------
f <- makeFun( x^2 ~ x )
f(3)
g <- makeFun( A*x^2 + B*x + C ~ x, A = 1, B = 2, C = 3 )
g(2)
g(2, A = 3, B = 2, C = 1)

## ------------------------------------------------------------------------
fprime <- D(f(x) ~ x)
fprime(2)
fprime
gprime <- D(g(x) ~ x)
gprime(3)
gprime(3, A = 3, B = 2, C = 1)
h <- makeFun( sin(x^2) ~ x )
hprime <- D( h(x) ~ x )
plotFun(hprime(x) ~ x, col = "red", x.lim = c(0,pi))
plotFun( h(x) ~ x, x.lim = c(0,pi), add = TRUE )

## ------------------------------------------------------------------------
plotFun(f(x) ~ x, type = "h")
F <- antiD( f(x) ~ x )
F(1) - F(0)

## ------------------------------------------------------------------------
t.test( ~ age, data = HELPrct )

## ------------------------------------------------------------------------
snippet("mom-beta01")   # to define beta.mom
results <- do(1000) * beta.mom(rbeta(50,2,5))
head(results, 2)
histogram( ~shape1, data = results, type = 'density', v = 2 )
histogram( ~shape2, data = results, type = 'density', v = 5 )

## ------------------------------------------------------------------------
ball.model <- lm( time ~ sqrt(height), data = balldrop)
time <- makeFun(ball.model)
time( height = 0.8 )
time( height = 0.8, interval = "confidence" )

## ------------------------------------------------------------------------
xyplot( time ~ height, data = balldrop )
plotFun( time(height) ~ height, add = TRUE )

## ------------------------------------------------------------------------
# TukeyHSD() can take a model created by lm()
model <- lm( pollution ~ location, data = airpollution)
TukeyHSD(model)
# we can even let TukeyHSD build the model for us
TukeyHSD( pollution ~ location, data = airpollution)

