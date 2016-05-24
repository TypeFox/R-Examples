require(PSAboot)
data(pisa.psa.cols)

##### United States
data(pisausa)
nrow(pisausa)
table(pisausa$PUBPRIV, useNA='ifany')
prop.table(table(pisausa$PUBPRIV, useNA='ifany')) * 100

bm.usa <- PSAboot(Tr=as.integer(pisausa$PUBPRIV) - 1,
			  Y=pisausa$Math,
			  X=pisausa[,pisa.psa.cols],
			  control.ratio=5, M=100, seed=2112)

(bootsum <- summary(bm.usa))
as.data.frame(bootsum)

plot(bm.usa)
plot(bm.usa, sort='none')
plot(bm.usa, sort='Stratification')
plot(bm.usa, sort='Matching')
plot(bm.usa, sort='MatchIt')

matrixplot(bm.usa)
boxplot(bm.usa)
hist(bm.usa)

(bm.usa.bal <- balance(bm.usa))
plot(bm.usa.bal) + geom_vline(xintercept=.1, linetype=2)
boxplot(bm.usa.bal)

##### Luxembourg
data(pisalux)
levels(pisalux$PUBPRIV)

t.test(Math ~ PUBPRIV, data=pisalux)
table(as.integer(pisalux$PUBPRIV) - 1)

bm.lux <- PSAboot(Tr=as.integer(pisalux$PUBPRIV) - 1,
				  Y=pisalux$Math,
				  X=pisalux[,pisa.psa.cols],
				  control.ratio=4, M=100, seed=2112)
summary(bm.lux)
plot(bm.lux)
matrixplot(bm.lux)
boxplot(bm.lux)
hist(bm.lux)

(bm.lux.bal <- balance(bm.lux))
plot(bm.lux.bal)
boxplot(bm.lux.bal)
