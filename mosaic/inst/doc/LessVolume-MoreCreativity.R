## ----include=FALSE----------------------------------------------------------------------
require(lubridate)
require(dplyr)
require(mosaic)
require(mosaicData)
trellis.par.set(theme=col.mosaic())
require(knitr)
opts_chunk$set(
  size='tiny', 
  tidy=FALSE,
  fig.width=6, 
  fig.height=3,
  fig.align="center",
  out.width="70%"
)
options(width=90)

## ---------------------------------------------------------------------------------------
require(mosaic)
require(mosaicData)

## ----eval=FALSE-------------------------------------------------------------------------
#  # simpler version
#  goal( ~ x, data = mydata )
#  # fancier version
#  goal( y ~ x | z , data = mydata )
#  # unified version
#  goal( formula , data = mydata )

## ----echo=FALSE, out.width="60%", out.height="35%"--------------------------------------
xyplot( births ~ date, data=Births78) 

## ----echo=TRUE, out.width="60%"---------------------------------------------------------
xyplot( births ~ date, data=Births78) 

## ----echo=FALSE-------------------------------------------------------------------------
bwplot( age ~ substance, data=HELPrct, xlab="substance" )

## ----echo=TRUE--------------------------------------------------------------------------
bwplot( age ~ substance, data=HELPrct)

## ----echo=FALSE-------------------------------------------------------------------------
bwplot( substance ~ age, data=HELPrct)

## ----echo=TRUE--------------------------------------------------------------------------
bwplot( substance ~ age, data=HELPrct )

## ---------------------------------------------------------------------------------------
histogram( ~ age, data=HELPrct) 

## ----eval=FALSE, tidy=FALSE-------------------------------------------------------------
#    histogram( ~age, data=HELPrct )
#  densityplot( ~age, data=HELPrct )
#       bwplot( ~age, data=HELPrct )
#       qqmath( ~age, data=HELPrct )
#  freqpolygon( ~age, data=HELPrct )
#     bargraph( ~sex, data=HELPrct )

## ----eval=FALSE, tidy=FALSE-------------------------------------------------------------
#  xyplot(  i1 ~ age,       data=HELPrct )
#  bwplot( age ~ substance, data=HELPrct )
#  bwplot( substance ~ age, data=HELPrct )

## ----eval=FALSE-------------------------------------------------------------------------
#  names(KidsFeet)    # 4th graders' feet
#  ?KidsFeet

## ----eval=FALSE-------------------------------------------------------------------------
#  names(Utilities)   # utility bill data
#  ?Utilities

## ----eval=FALSE-------------------------------------------------------------------------
#  require(NHANES)    # load package
#  names(NHANES)      # body shape, etc.
#  ?NHANES

## ---- tidy=FALSE------------------------------------------------------------------------
densityplot( ~ age | sex, data=HELPrct,  
               groups=substance,  
               auto.key=TRUE)   

## ----out.width="85%", tidy=FALSE--------------------------------------------------------
require(lubridate)
xyplot( births ~ date, data=Births78,  
  groups=wday(date, label=TRUE, abbr=TRUE), 
  type='l',
  auto.key=list(columns=4, lines=TRUE, points=FALSE),
  par.settings=list(
    superpose.line=list( lty=1 ) ))

## ----eval=FALSE, include=FALSE----------------------------------------------------------
#  xyplot( births ~ date, data=Births78,
#          groups=wday(date, label=TRUE, abbr=TRUE), type='l',
#          auto.key=list(columns=4),
#          superpose.symbol=list(
#              pch=16, cex=1.2, alpha=.8)))

## ----fig.show='hold'--------------------------------------------------------------------
histogram( ~ age, data=HELPrct )  # width=5 (or 10) might be good here
     mean( ~ age, data=HELPrct )

## ---------------------------------------------------------------------------------------
favstats( ~ age, data=HELPrct )

## ---------------------------------------------------------------------------------------
tally( ~ sex, data=HELPrct)
tally( ~ substance, data=HELPrct)

## ----eval = FALSE-----------------------------------------------------------------------
#  sd(   age ~ substance, data=HELPrct )
#  sd( ~ age | substance, data=HELPrct )
#  sd( ~ age, groups=substance, data=HELPrct )

## ---- echo=FALSE------------------------------------------------------------------------
sd( ~ age, groups=substance, data=HELPrct )

## ---------------------------------------------------------------------------------------
tally( sex ~ substance, data=HELPrct )
tally( ~ sex + substance, data=HELPrct )

## ---------------------------------------------------------------------------------------
tally( sex ~ substance,   data=HELPrct, format="proportion" )
tally( substance ~ sex,   data=HELPrct, format="proportion", margins=TRUE )
tally( ~ sex + substance, data=HELPrct, format="proportion", margins=TRUE )
tally( sex ~ substance,   data=HELPrct, format="percent" )

## ----echo=FALSE-------------------------------------------------------------------------
HELPrct <- transform(HELPrct, sex=factor(sex, labels=c('F','M')),
                     substance = factor(substance, labels=c('A', 'C', 'H')))

## ----size='small'-----------------------------------------------------------------------
mean( age ~ substance | sex, data=HELPrct )
mean( age ~ substance | sex, data=HELPrct, .format="table" )

## ----echo=FALSE-------------------------------------------------------------------------
rm(HELPrct)
data(HELPrct)

## ----eval=FALSE-------------------------------------------------------------------------
#    mean( age ~ sex, data=HELPrct )
#  bwplot( age ~ sex, data=HELPrct )
#      lm( age ~ sex, data=HELPrct )

## ----echo=FALSE-------------------------------------------------------------------------
  mean( age ~ sex, data=HELPrct )
    coef(lm( age ~ sex, data=HELPrct ))

## ---------------------------------------------------------------------------------------
xpnorm( 700, mean=500, sd=100)

## ---------------------------------------------------------------------------------------
xpnorm( c(300, 700), mean=500, sd=100)

## ---- echo=FALSE------------------------------------------------------------------------
phs <- cbind(c(104,189),c(10933,10845))
colnames(phs) <- c("heart attack","no heart attack")
rownames(phs) <- c("aspirin","placebo")

## ---------------------------------------------------------------------------------------
xchisq.test(phs)

## ---------------------------------------------------------------------------------------
model <- lm(width ~ length * sex, 
            data=KidsFeet)
Width <- makeFun(model)
Width( length=25, sex="B")
Width( length=25, sex="G")

## ---- include=FALSE---------------------------------------------------------------------
trellis.par.set(
  superpose.symbol=list(col=c('navy','red'), pch=16), 
  superpose.line=list(lty=1, col=c('navy','red'))
)

## ---- fig.keep='last'-------------------------------------------------------------------
xyplot( width ~ length, data=KidsFeet, 
        groups=sex, auto.key=TRUE )
plotFun( Width(length, sex="B") ~ length, 
         col=1, add=TRUE)
plotFun( Width(length, sex="G") ~ length, 
         col=2, add=TRUE)

## ---- include=FALSE---------------------------------------------------------------------
trellis.par.set(theme=col.mosaic())

## ----echo=FALSE-------------------------------------------------------------------------
require(mosaic)
trellis.par.set(theme=col.mosaic())
require(knitr)
opts_chunk$set(size='small', cache=TRUE)
options(width=90)
set.seed(12345)

