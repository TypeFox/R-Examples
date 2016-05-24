require(multilevelPSA)
require(party)
data(pisana)

prop.table(table(pisana$CNT, pisana$PUBPRIV, useNA='ifany'), 1) * 100

cnt = 'USA' #Can change this to USA, MEX, or CAN
pisana2 = pisana[pisana$CNT == cnt,]

prop.table(table(pisana2$PUBPRIV, useNA='ifany')) * 100

pisana2$treat <- as.integer(pisana2$PUBPRIV) %% 2

lr.results <- glm(treat ~ ., data=pisana2[,c('treat',pisa.psa.cols)], family='binomial')
st = data.frame(ps=fitted(lr.results), 
				math=apply(pisana2[,paste('PV', 1:5, 'MATH', sep='')], 1, mean), 
				pubpriv=pisana2$treat)
st$treat = as.logical(st$pubpriv)

loess.plot(st$ps, response=st$math, treatment=st$treat, 
		   percentPoints.control = 0.4, percentPoints.treat=0.4,
		   plot.strata=10, method='loess')

