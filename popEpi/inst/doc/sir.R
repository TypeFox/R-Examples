## ---- echo=TRUE, warning=FALSE, message=FALSE----------------------------
library(popEpi)
library(Epi)
library(splines)

## ------------------------------------------------------------------------
data(sire)
data(popmort)
c <- lexpand( sire, status = status, birth = bi_date, exit = ex_date, entry = dg_date,
              breaks = list(per = 1950:2013, age = 1:100, fot = c(0,10,20,Inf)), 
              aggre = list(fot, agegroup = age, year = per, sex) )

se <- sir( coh.data = c, coh.obs = 'from0to2', coh.pyrs = 'pyrs',
           ref.data = popmort, ref.rate = 'haz', 
           adjust = c('agegroup','year','sex'), print ='fot')
se

## ------------------------------------------------------------------------
c <- lexpand( sire, status = status %in% 1:2, birth = bi_date, exit = ex_date, entry = dg_date,
              breaks = list(per = 1950:2013, age = 1:100, fot = c(0,10,20,Inf)), 
              aggre = list(fot, agegroup = age, year = per, sex) )

se <- sir( coh.data = c, coh.obs = 'from0to1', coh.pyrs = 'pyrs',
           ref.data = popmort, ref.rate = 'haz', 
           adjust = c('agegroup','year','sex'), print ='fot')
se

## ---- fig.height=3, fig.width=6------------------------------------------
plot(se, col = 2:3)
title('SMR for follow-up categories')

## ---- fig.height=5, fig.width=6------------------------------------------
c <- lexpand( sire, status = status %in% 1:2, birth = bi_date, exit = ex_date, entry = dg_date,
              breaks = list(per = 1950:2013, age = 1:100, fot = 0:50), 
              aggre = list(fot, agegroup = age, year = per, sex) )

sf <- sirspline( coh.data = c, coh.obs = 'from0to1', coh.pyrs = 'pyrs', 
                 ref.data = popmort, ref.rate = 'haz', 
                 adjust = c('agegroup','year','sex'),
                 spline = c('agegroup','fot'), dependent.splines=FALSE)

st <- sirspline( coh.data = c, coh.obs = 'from0to1', coh.pyrs = 'pyrs', 
                 ref.data = popmort, ref.rate = 'haz', 
                 adjust = c('agegroup','year','sex'),
                 spline = c('agegroup','fot'), dependent.splines = TRUE)

plot(sf, col=2, log=TRUE)
title('Splines fitted in different models')

plot(st, col=4, log=TRUE)
title('Splines are dependent')

## ---- results='hide', fig.height=5, fig.width=6--------------------------
c$year.cat <- ifelse(c$year < 2002, 1, 2)
sy <- sirspline( coh.data = c, coh.obs = 'from0to1', coh.pyrs = 'pyrs', 
                 ref.data = popmort, ref.rate = 'haz', 
                 adjust = c('agegroup','year','sex'),
                 spline = c('agegroup'), print = 'year.cat')
plot(sy, log=TRUE)
legend('topright', c('before 2002','after 2002'), lty=1, col=c(1,2))

## ------------------------------------------------------------------------
print(sy)

