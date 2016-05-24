## ----setupRnw, include=FALSE----------------------------------------
# source('../setup.R')
## @knitr setup
#setCacheDir("cache")
require(grDevices) 
require(datasets) 
require(stats) 
require(lattice)
require(grid) 
require(mosaic) 
require(mosaicData) 
trellis.par.set(theme=col.mosaic(bw=FALSE))
trellis.par.set(fontsize=list(text=9))
options(keep.blank.line=FALSE) 
options(width=70)
require(vcd)
require(knitr)

opts_chunk$set( 
  fig.align="center",
  fig.show="hold"
)

## ----eval=FALSE, echo=TRUE------------------------------------------
#  install.packages("mosaic")

## ----moreSetup------------------------------------------------------
require(mosaic)
require(mosaicData) 
options(digits=3)

## -------------------------------------------------------------------
Mustangs <- read.file("http://www.mosaic-web.org/go/datasets/MustangPrice.csv")

## ----dot, fig.width=4, fig.height=2---------------------------------
histogram(~Price, data=Mustangs)

## ----meanmust-------------------------------------------------------
mean(~Price, data=Mustangs)

## ----eval=FALSE-----------------------------------------------------
#  mean(Mustangs$Price)

## ----eval=FALSE-----------------------------------------------------
#  mean(~Price, data=Mustangs)

## -------------------------------------------------------------------
simple=c(1, 2, 3, 4, 5)
resample(simple)
resample(simple)
resample(simple)

## ----results="hide"-------------------------------------------------
resample(Mustangs)

## ----echo=FALSE-----------------------------------------------------
set.seed(104)
head(resample(Mustangs))
cat("... and so on")

## -------------------------------------------------------------------
mean(~Price, data=resample(Mustangs))

## -------------------------------------------------------------------
mean(~Price, data=resample(Mustangs))

## -------------------------------------------------------------------
do(5) * mean(~Price, data=resample(Mustangs))

## -------------------------------------------------------------------
Mustangs.Price.boot <- do(1000) * mean(~Price, data=resample(Mustangs))

## ----hist, fig.width=4, fig.height=2--------------------------------
histogram(~ mean, data=Mustangs.Price.boot, xlab="Mean Mustang Price (in thousand dollars)")

## ----confint--------------------------------------------------------
confint(Mustangs.Price.boot, level=0.90, method="quantile")
confint(Mustangs.Price.boot, level=0.90, method="stderr")

## ----qdata----------------------------------------------------------
qdata(~mean, c(.05, .95), data=Mustangs.Price.boot)
# alternative
cdata(~mean, .90, data=Mustangs.Price.boot)


## ----tstar----------------------------------------------------------
tstar <- qt(.95, df=24)
zstar <- qnorm(.95)

## ----margin---------------------------------------------------------
tstar * sd(~mean, data=Mustangs.Price.boot)
zstar * sd(~mean, data=Mustangs.Price.boot)

## ----proptable------------------------------------------------------
prop(~ rbinom(1000, prob=0.5, size=428) >= 240)

## ----proptable2-----------------------------------------------------
prop(~ rbinom(1000, prob=0.5, size=428) >= 240)

## ----pbinom---------------------------------------------------------
pbinom(239, prob=0.5, size=428)

## -------------------------------------------------------------------
binom.test(240, 248)

## ----coinflip-------------------------------------------------------
do(1) * rflip(428)

## ----flips1, fig.width=4, fig.height=3------------------------------
NFL.null <- do(1000) * rflip(428)
prop(~ heads >= 240, data=NFL.null)

## ----flips, fig.width=4, fig.height=2-------------------------------
histogram(~ heads, groups=(heads >= 240), data=NFL.null)

## ----fetchsleep-----------------------------------------------------
Sleep <- read.file("http://www.mosaic-web.org/go/datasets/SleepCaffeine.csv")

## ----obsmean--------------------------------------------------------
mean(Words ~ Group, data=Sleep)
obs <- diff(mean(Words ~ Group, data=Sleep))
obs

## ----fig.width=4, fig.height=2--------------------------------------
bwplot(Group ~ Words, data=Sleep)

## -------------------------------------------------------------------
diff(mean(Words ~ shuffle(Group), data=Sleep))

## -------------------------------------------------------------------
diff(mean(Words ~ shuffle(Group), data=Sleep))

## ----setseed134, echo=FALSE-----------------------------------------
set.seed(134) # make sure the result below is "typical"

## ----sleep, tidy=FALSE, fig.width=4, fig.height=2-------------------
Sleep.null <- do(1000) * diff(mean(Words ~ shuffle(Group), data=Sleep))
histogram(~ Sleep, groups=(Sleep >= obs),  data=Sleep.null, width=0.4,
  xlab="Distribution of difference in means\nunder the null hypothesis")

## -------------------------------------------------------------------
cor(Price, Miles, data=Mustangs)

## -------------------------------------------------------------------
Mustangs.cor.boot <- do(1000) * cor(Price, Miles, data=resample(Mustangs))
quantiles <- qdata(~cor, c(.025, .975), data=Mustangs.cor.boot); quantiles

## ----cor, fig.width=4, fig.height=2, tidy=FALSE---------------------
histogram(~ cor, data=Mustangs.cor.boot,
  groups=cut(cor, c(-Inf, quantiles$quantile, Inf)),
  n=30) 
confint(Mustangs.cor.boot)

## ----price-mileage-graph, fig.width=4, fig.height=3-----------------
xyplot( Price ~ Miles, data=Mustangs )

## ----mustangregression----------------------------------------------
lm(Price ~ Miles, data=Mustangs)

## -------------------------------------------------------------------
mean( ~Price, data=Mustangs )

## -------------------------------------------------------------------
lm( Price ~ 1, data=Mustangs)

## -------------------------------------------------------------------
KidsFeet %>% select(-name, -birthmonth) %>% rescale() -> KidsFeet2
KidsFeet %>% select(-name, -birthmonth) %>% rescale() -> KidsFeet2
mean( Price ~ 1, data=Mustangs)

## -------------------------------------------------------------------
mean( Words ~ 1, data=Sleep )

## -------------------------------------------------------------------
mean( Words ~ Group, data=Sleep )

## -------------------------------------------------------------------
lm( Words ~ Group, data=Sleep )

## -------------------------------------------------------------------
diffmean( Words ~ Group, data=Sleep )

## -------------------------------------------------------------------
prop( homeless ~ 1, data=HELPrct)

## -------------------------------------------------------------------
prop( homeless ~ sex, data=HELPrct )


## -------------------------------------------------------------------
diffprop( homeless ~ sex, data=HELPrct )

## -------------------------------------------------------------------
lm( homeless=="homeless" ~ 1, data=HELPrct )

## -------------------------------------------------------------------
lm(homeless=="homeless" ~ sex, data=HELPrct)

## -------------------------------------------------------------------
Mustangs.lm.boot <- do(1000) * lm(Price ~ Miles, data=resample(Mustangs))
confint(Mustangs.lm.boot)

## -------------------------------------------------------------------
HELPrct.null <- do(1000) * lm(homeless=="homeless" ~ shuffle(sex), data=HELPrct)
prop(~ abs(sex.male) > 0.1146, data=HELPrct.null)

## -------------------------------------------------------------------
Mustangs.boot1 <- do(1000) * lm( Price ~ Age, data=resample(Mustangs))
Mustangs.boot2 <- do(1000) * lm( Price ~ Miles, data=resample(Mustangs))
Mustangs.boot3 <- do(1000) * lm( Price ~ Miles + Age, data=resample(Mustangs))

## -------------------------------------------------------------------
confint(Mustangs.boot1)

## -------------------------------------------------------------------
confint(Mustangs.boot2)

## -------------------------------------------------------------------
confint(Mustangs.boot3)

## -------------------------------------------------------------------
anova( lm( Price ~ Miles + Age, data=Mustangs))

## -------------------------------------------------------------------
anova( lm( Price ~ Age + Miles, data=Mustangs))

## -------------------------------------------------------------------
do(1) * lm( Price ~ Miles, data=Mustangs)
do(1) * lm( Price ~ Miles + Age, data=Mustangs)

## -------------------------------------------------------------------
Mustangs.Age.boot <- do(1000) * lm( Price ~ Miles + shuffle(Age), data=Mustangs )
confint(Mustangs.Age.boot)

## -------------------------------------------------------------------
Mustangs.Miles.boot <- do(1000) * lm( Price ~ shuffle(Miles) + Age, data=Mustangs )
confint(Mustangs.Miles.boot)

## -------------------------------------------------------------------
chisq.test( tally( ~ homeless + sex, 
                   data=HELPrct, margins=FALSE))

## -------------------------------------------------------------------
pval( chisq.test( tally( ~ homeless + sex, 
                   data=HELPrct, margins=FALSE)) )

## -------------------------------------------------------------------
pval( chisq.test( tally( ~ shuffle(homeless) + sex, 
                         data=HELPrct, margins=FALSE)))

## -------------------------------------------------------------------
Chisq.null <- do(1000)* pval( chisq.test( tally( ~ shuffle(homeless) + sex, 
                         data=HELPrct, margins=FALSE)))

## ----fig.width=4, fig.height=2--------------------------------------
prop( ~(p.value < 0.05), data=Chisq.null)
histogram( ~p.value, data=Chisq.null, width=0.10 )
qqmath( ~p.value, data=Chisq.null, dist=qunif)

## -------------------------------------------------------------------
HELP.logistic.boot <- do(1000) * 
   glm( homeless=="homeless" ~ age + sex, 
     data=resample(HELPrct), family="binomial")
confint(HELP.logistic.boot)

