## ----setup,echo=FALSE,message=FALSE,include=FALSE------------------------
#source('setup.R')
require(mosaic)
require(mosaicData)
require(parallel)
trellis.par.set(theme=col.mosaic())
set.seed(123)
#knit_hooks$set(inline = function(x) {
#	if (is.numeric(x)) return(knitr:::format_sci(x, 'latex'))
#	x = as.character(x)
#	h = knitr:::hilight_source(x, 'latex', list(prompt=FALSE, size='normalsize'))
#	h = gsub("([_#$%&])", "\\\\\\1", h)
#	h = gsub('(["\'])', '\\1{}', h)
#	gsub('^\\\\begin\\{alltt\\}\\s*|\\\\end\\{alltt\\}\\s*$', '', h)
#})
opts_chunk$set(
	dev="pdf",
	eval=FALSE,
	tidy=FALSE,
	fig.align='center',
	fig.show='hold',
	message=FALSE
	)

## ------------------------------------------------------------------------
#  apropos()
#  ?
#  ??
#  example()

## ------------------------------------------------------------------------
#  # basic ops: + - * / ^ ( )
#  log(); exp(); sqrt()

## ----highlight=FALSE-----------------------------------------------------
#  log10(); abs(); choose()

## ------------------------------------------------------------------------
#  goal(y ~ x | z, data=...,
#                   groups=...)

## ------------------------------------------------------------------------
#  favstats()   # mosaic
#  tally()      # mosaic
#  mean()       # mosaic augmented
#  median()     # mosaic augmented
#  sd()         # mosaic augmented
#  var()        # mosaic augmented
#  diffmean()   # mosaic

## ----highlight=FALSE-----------------------------------------------------
#  quantile()   # mosaic augmented
#  prop()       # mosaic
#  perc()       # mosaic
#  rank()
#  IQR()        # mosaic augmented
#  min(); max() # mosaic augmented

## ------------------------------------------------------------------------
#  bwplot()
#  xyplot()
#  histogram() # mosaic augmented
#  densityplot()
#  freqpolygon() # mosaic
#  qqmath()
#  makeFun()   # mosaic
#  plotFun()   # mosaic

## ----highlight=FALSE-----------------------------------------------------
#  ladd()      # mosaic
#  dotPlot()   # mosaic
#  bargraph()  # mosaic
#  xqqmath()   # mosaic

## ----eval=FALSE----------------------------------------------------------
#  mplot(data=HELPrct, 'scatter')
#  mplot(data=HELPrct, 'boxplot')
#  mplot(data=HELPrct, 'histogram')

## ------------------------------------------------------------------------
#  rflip()     # mosaic
#  do()        # mosaic
#  sample()    # mosaic augmented
#  resample()  # with replacement
#  shuffle()   # mosaic

## ----highlight=FALSE-----------------------------------------------------
#  rbinom()
#  rnorm()     # etc, if needed

## ------------------------------------------------------------------------
#  pbinom(); pnorm();
#  xpnorm()    # mosaic augmented
#  pchisq(); pt()
#  qbinom(); qnorm();
#  qchisq(); qt()
#  plotDist()  # mosaic

## ------------------------------------------------------------------------
#  t.test()      # mosaic augmented
#  binom.test()  # mosaic augmented
#  prop.test()   # mosaic augmented
#  xchisq.test() # mosaic
#  fisher.test()
#  pval()        # mosaic
#  model <- lm() # linear models
#  summary(model)
#  coef(model)
#  confint(model) # mosaic augmented
#  anova(model)
#  makeFun(model) # mosaic
#  resid(model); fitted(model)
#  mplot(model)   # mosaic

## ----highlight=FALSE-----------------------------------------------------
#  mplot(TukeyHSD(model))
#  model <- glm() # logistic reg.

## ------------------------------------------------------------------------
#  require(mosaicData)  # load package
#  read.file()    # mosaic
#  nrow(); ncol(); dim()
#  summary()
#  str()
#  names()
#  head(); tail()
#  with()
#  factor()

## ----highlight=FALSE-----------------------------------------------------
#  ntiles()      # mosaic
#  cut()
#  c()
#  cbind(); rbind()
#  colnames()
#  rownames()
#  relevel()
#  reorder()

## ----highlight=FALSE-----------------------------------------------------
#  rep()
#  seq()
#  sort()
#  rank()

## ----highlight=FALSE-----------------------------------------------------
#  select()      # dplyr
#  mutate()      # dplyr
#  filter()      # dplyr
#  arrange()     # dplyr
#  summarise()   # dplyr
#  group_by()    # dplyr
#  left_join()   # dplyr
#  inner_join()  # dplyr
#  merge()

## ----more-hooks,eval=TRUE,echo=FALSE-------------------------------------
opts_chunk$set(
	eval=TRUE, 
  size='small',
	fig.width=4,
	fig.height=1.9,
	fig.align="center",
	out.width=".25\\textwidth",
	out.height=".125\\textwidth",
	tidy=TRUE,
	comment=NA
)

## ----echo=FALSE-----------------------
options(width=40)
options(show.signif.stars=FALSE)

## ----coins,fig.keep="last"------------
rflip(6)
do(2) * rflip(6)
coins <- do(1000)* rflip(6)
tally(~ heads, data=coins)

## -------------------------------------
tally(~ heads, data=coins, format="perc")
tally(~ (heads>=5 | heads<=1) , data=coins)

## ----coins-hist,fig.keep="last"-------
histogram(~ heads, data=coins, width=1,
            groups = (heads>=5 | heads<=1))

## ----tally----------------------------
require(mosaicData)   # this package contains data sets used below
tally(sex ~ substance, data=HELPrct)
mean(age ~ sex, data=HELPrct)
diffmean(age ~ sex, data=HELPrct)
favstats(age ~ sex, data=HELPrct)

## ----densityplot,fig.height=2.4-------
densityplot(~ age | sex, groups=substance, 
               data=HELPrct, auto.key=TRUE)

## ----bwplot---------------------------
bwplot(age ~ substance | sex, data=HELPrct)

## ----message=FALSE--------------------
pval(binom.test(~ sex, data=HELPrct))
confint(t.test(~ age, data=HELPrct))

## ----tidy=FALSE-----------------------
model <- lm(age ~ sex + substance, 
			data=HELPrct) 
anova(model)

## ----tidy=FALSE-----------------------
xyplot(Sepal.Length ~ Sepal.Width, 
        groups=Species, data=iris) 

## ----fig.keep="last", tidy=FALSE, fig.height=2.3----
model <- lm(length ~ width + sex, 
            data=KidsFeet)
ln <- makeFun(model)
ln( width=8.25, sex="B")
xyplot(length ~ width, groups=sex, 
       data=KidsFeet)
plotFun(ln(w, sex="B") ~ w, 
        add=TRUE, col="skyblue")
plotFun(ln(w,sex="G") ~ w, 
        add=TRUE, col="navy")

## ----fig.height=1.75------------------
plotDist("chisq", df=4)

## ----include=FALSE--------------------
tally(homeless ~ sex, data=HELPrct)

## ----include=FALSE--------------------
chisq.test(tally(homeless ~ sex, data=HELPrct))
prop.test(homeless ~ sex, data=HELPrct)

