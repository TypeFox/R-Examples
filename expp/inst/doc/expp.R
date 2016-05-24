## ----eval=FALSE----------------------------------------------------------
#  install.packages("expp")

## ----echo=FALSE, message=FALSE-------------------------------------------
require(rgeos); require(sp); require(spdep); require(spatstat); require(deldir)
par(mar = c(0,0,1,0) )

## ----,message=FALSE------------------------------------------------------
require(expp)

## ----eval=FALSE----------------------------------------------------------
#  help(bluetit_breeding)
#  help(bluetit_epp)
#  help(bluetit_boundary)

## ----, eval = FALSE------------------------------------------------------
#  data(bluetit_breeding)
#  head(bluetit_breeding[bluetit_breeding$year_ == 2011, ])

## ----, echo=FALSE, results='asis'----------------------------------------
data(bluetit_breeding)
knitr::kable(head(bluetit_breeding[bluetit_breeding == 2011, ]))
data(bluetit_breeding)

## ----, eval = FALSE------------------------------------------------------
#  data(bluetit_epp)
#  head(bluetit_epp[bluetit_epp$year_ == 2011, ])

## ----, echo=FALSE, results='asis'----------------------------------------
data(bluetit_epp)
knitr::kable(head(bluetit_epp[bluetit_epp == 2011, ]))

## ------------------------------------------------------------------------
data(bluetit_boundary)
summary(bluetit_boundary)

## ------------------------------------------------------------------------
b = split(bluetit_breeding, bluetit_breeding$year_)
e = split(bluetit_epp, bluetit_epp$year_) 

# sample sizes by year

# number of breeding pairs
sapply(b, nrow)

# number of extra-pair events
sapply(e, nrow)

# For the sake of conciseness only two years are used in the folowing analyses
b = b[c("2009", "2010")]
e = e[c("2009", "2010")]
p = bluetit_boundary[bluetit_boundary$year_ %in% c("2009", "2010"), ]


## ----tidy=FALSE----------------------------------------------------------
breedingDat = lapply(b, SpatialPointsBreeding, coords= ~x+y, id='id', breeding= ~male + female, 
  proj4string = CRS(proj4string(p)))

eppDat = lapply(e, eppMatrix, pairs = ~ male + female)



## ------------------------------------------------------------------------
polygonsDat = mapply(DirichletPolygons, x = breedingDat, boundary = split(p, p$year_)) 

## ------------------------------------------------------------------------
maxlag = 10
eppOut = mapply(FUN = epp, breedingDat, polygonsDat, eppDat, maxlag)

## ----results='hide', dpi=100, fig.width=7, fig.height=10, fig.align='left', warning=FALSE----
op = par(mar = c(0,0,2,0))

for(year in c("2009", "2010") ) { 
  plot(eppOut[[year]], cex = 0.7, lwd = .5, border = "navy" )
  title(main = year)
  }



## ----, warning=FALSE, dpi=100, fig.width=6, fig.height=8, fig.align='left'----
year = '2010'
box = 110
eppOut10 = eppOut[[year]]
plot(eppOut10 , zoom = box, maxlag = 2,cex = .7,  border = 'white', col = 'grey70', zoom.col = "bisque")

par(op)

## ----results='hide',fig.width=8, fig.height=6----------------------------
op = par(mfrow = c(1,2))
    
barplot(eppOut[[1]],relativeValues = TRUE, main = 2009) 
legend(x="topright", legend = c('Observed', 'Potential'), lty = c(1, 2),bty='n')
barplot(eppOut[[2]], relativeValues = TRUE, main = 2010)

par(op)

## ------------------------------------------------------------------------
dat = lapply(eppOut, as.data.frame) # a list of data.frame(s)
dat = do.call(rbind, dat)
dat$year_ = dat$year__MALE; dat$year__FEMALE = NULL


## ------------------------------------------------------------------------
dat$rank = dat$rank - min(dat$rank)
table(dat$rank)

## ------------------------------------------------------------------------
center = function(x) { return(x - mean(x, na.rm = TRUE)) }
scale2 = function(x) { return(x/(2*sd(x, na.rm = TRUE))) }

# Compute asynchrony
dat$asynchrony = abs(dat$layingDate_MALE - dat$layingDate_FEMALE)

#a Compute relative within-rank asynchrony
MALE_splitBy = paste(dat$year_, dat$id_MALE, dat$male, dat$rank, sep = "_")
dat$relative_asynchrony_MALE = unsplit(lapply(split(dat$asynchrony, MALE_splitBy), center), MALE_splitBy)
dat$relative_asynchrony_MALE = scale2(dat$relative_asynchrony_MALE)

FEMALE_splitBy = paste(dat$year_, dat$id_FEMALE, dat$female, dat$rank, sep = "_")
dat$relative_asynchrony_FEMALE = unsplit(lapply(split(dat$asynchrony, FEMALE_splitBy), center), FEMALE_splitBy)
dat$relative_asynchrony_FEMALE = scale2(dat$relative_asynchrony_FEMALE)

## ----eval=FALSE----------------------------------------------------------
#  table(dat$epp, dat$year_) #extra-pair frequency by year.

## ----echo = FALSE, results='asis'----------------------------------------
knitr::kable(table(dat$epp, dat$year_))

## ----eval=FALSE----------------------------------------------------------
#  require(lme4)
#  fm = glmer(epp ~ rank + male_age_MALE + relative_asynchrony_MALE + relative_asynchrony_FEMALE +
#               (1|male) + (1|female) + (1|year_), data = dat, family = binomial)
#  summary(fm)

