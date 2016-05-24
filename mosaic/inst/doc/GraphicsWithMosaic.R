## ----include=FALSE-------------------------------------------------------
# Don't delete this chunk if you are using the mosaic package
# This loads the mosaic and dplyr packages
require(mosaic)
require(mosaicData)
require(NHANES)

## ----include=FALSE-------------------------------------------------------
# Some customization.  You can alter or delete as desired (if you know what you are doing).

# This changes the default colors in lattice plots.
trellis.par.set(theme=theme.mosaic())  

# knitr settings to control how R chunks work.
require(knitr)
set.seed(1)
opts_chunk$set(
  tidy=FALSE,     # display code as typed
  size="small",   # slightly smaller font for code
  fig.width=5,
  fig.height=3
)

## ------------------------------------------------------------------------
histogram(~ rbinom( 500, 20, .3), width=1, fit="normal", v=c(6,10), h=0.1 )

## ---- fig.keep="last"----------------------------------------------------
xyplot( rnorm(100) ~ rnorm(100) )
ladd( grid.text("Here is some text", x=0, y=0, default.units="native") )
ladd( panel.abline( a=0, b=1, col="red", lwd=3, alpha=.4 ) )
ladd( panel.rect(x=-1, y=-1, width=1, height=1, col="gray80", fill="lightsalmon"))
ladd( panel.rect(x=0, y=0, width=2, height=2, col="gray80", fill="lightskyblue"), 
      under=TRUE)

## ---- fig.height=8-------------------------------------------------------
require(gridExtra)
mod <- lm( width ~ length * sex, data=KidsFeet )
mplot(mod, which=1:7, multiplot = TRUE, ncol=2)

## ---- fig.height=8-------------------------------------------------------
mplot(mod, which=1:7, system="ggplot", ncol=2)

## ------------------------------------------------------------------------
mplot(mod, which=7)
mplot(mod, which=7, rows=-1)
mplot(mod, which=7, rows=c("sexG", "length", "length:sexG"), 
      title="Custom titles are supported")

## ------------------------------------------------------------------------
mod <- lm(age ~ substance, data=HELPrct)
mplot(TukeyHSD(mod))
mplot(TukeyHSD(mod), system="ggplot")

## ---- fig.keep="last"----------------------------------------------------
mod <- lm(width ~ length* sex, data = KidsFeet)
L <- makeFun(mod)
L( length=15, sex="B")
L( length=15, sex="G")
xyplot(width ~ length, groups = sex, data = KidsFeet, auto.key=TRUE)
plotFun( L(length, sex="B") ~ length, add=TRUE, col=1 )
plotFun( L(length, sex="G") ~ length, add=TRUE, col=2 )

## ---- label="logistic", fig.keep="last", fig.height=5--------------------
mod <- glm( SmokeNow =="Yes" ~ Age + Race3, data=NHANES, family=binomial())
SmokerProb <- makeFun(mod)
xyplot( SmokeNow=="Yes" ~ Age, groups=Race3, data=NHANES, alpha=.01, xlim=c(20,90) )
plotFun(SmokerProb(Age, Race3="Black") ~ Age, col="black", add=TRUE)
plotFun(SmokerProb(Age, Race3="White") ~ Age, col="red", add=TRUE) 
ladd(grid.text("Black", x=25, y=SmokerProb(25, Race="Black"),hjust = 0, vjust=-0.2,
               gp=gpar(col="black"),
               default.units="native"))
ladd(grid.text("White", x=25, y=SmokerProb(25, Race="White"),hjust = 0, vjust=-0.2,
               gp=gpar(col="red"),
               default.units="native"))

## ------------------------------------------------------------------------
f <- makeFun(sin(x) ~ x)
plotFun( f(x) ~ x, xlim = c( -2 * pi, 2 * pi) )

## ------------------------------------------------------------------------
plotFun( x * sin(1/x) ~ x, xlim=c(-1,1) )
plotFun( x * sin(1/x) ~ x, xlim=c(-1,1), npts=10000 )

## ------------------------------------------------------------------------
plotDist("chisq", df=3)
plotDist("chisq", df=3, kind="cdf")

## ------------------------------------------------------------------------
xpnorm(80, mean=100, sd=15)
xpnorm(c(80,120), mean=100, sd=15)

## ------------------------------------------------------------------------
pdist("chisq", 4, df=3)
pdist("f", 3, df1=5, df2=20)
qdist("t", c(.025, .975) , df=5)

## ------------------------------------------------------------------------
histogram( ~ rbinom(1000, 20, .4), width=1, v=20 * .4 )
SD <- sqrt(20 * .4 * .6)
plotDist("norm", mean=.4*20, sd=SD, add=TRUE, alpha=.7)

## ---- fig.keep="last"----------------------------------------------------
plotDist("norm", col="blue", mean=2, xlim=c(-4,8))
plotDist("norm", mean=5, col="green", kind='histogram', add=TRUE)  # add, overtop
plotDist("norm", mean=0, col="red", kind='histogram', under=TRUE)  # add, but underneath!

## ------------------------------------------------------------------------
# we need to get state names into the data frame and then fix two of them with
# wrong state abbreviations.  Then we are ready to make maps
sAnscombe <- car::Anscombe %>%
group_by(state = rownames(car::Anscombe)) %>%
summarise(income = sum(income)) %>%
mutate(state = standardName(state, c(IO = "IA", KA = "KS"), quiet=TRUE))

mUSMap(sAnscombe, key="state", fill="income")
mUSMap(sAnscombe, key="state", fill="income", style="real") 

# A sillier example
if (require(mapproj)) {
Countries %>% mutate( nletters = nchar(gapminder) ) %>%
  mWorldMap( key="gapminder", fill="nletters") + coord_map()
} else {
Countries %>% mutate( nletters = nchar(gapminder) ) %>%
  mWorldMap( key="gapminder", fill="nletters") 
}

