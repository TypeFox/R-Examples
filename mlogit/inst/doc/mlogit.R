### R code from vignette source 'mlogit.rnw'

###################################################
### code chunk number 1: mlogit.rnw:175-176
###################################################
options(prompt= "R> ", useFancyQuotes = FALSE)


###################################################
### code chunk number 2: mlogit.rnw:181-182
###################################################
library("mlogit")


###################################################
### code chunk number 3: mlogit.rnw:224-226
###################################################
data("Fishing", package = "mlogit")
head(Fishing, 3)


###################################################
### code chunk number 4: mlogit.rnw:239-241
###################################################
data("TravelMode", package="AER")
head(TravelMode)


###################################################
### code chunk number 5: mlogit.rnw:260-261
###################################################
Fish <- mlogit.data(Fishing, shape="wide", varying=2:9, choice="mode")


###################################################
### code chunk number 6: mlogit.rnw:274-275
###################################################
head(Fish, 5)


###################################################
### code chunk number 7: mlogit.rnw:287-288
###################################################
head(index(Fish))


###################################################
### code chunk number 8: mlogit.rnw:313-316
###################################################
TM <- mlogit.data(TravelMode,
choice = "choice", shape = "long", alt.levels = c("air", "train",
"bus", "car"))


###################################################
### code chunk number 9: mlogit.rnw:323-325
###################################################
TM <- mlogit.data(TravelMode ,choice = "choice", shape = "long",
                  alt.var = "mode")


###################################################
### code chunk number 10: mlogit.rnw:332-335
###################################################
TM <- mlogit.data(TravelMode, choice = "choice",
                  shape = "long", chid.var = "individual",
                  alt.levels = c("air", "train", "bus", "car"))


###################################################
### code chunk number 11: mlogit.rnw:340-342
###################################################
TM <- mlogit.data(TravelMode, choice = "choice", shape = "long",
                 chid.var = "individual", alt.var = "mode")


###################################################
### code chunk number 12: mlogit.rnw:348-351
###################################################
TM <- mlogit.data(TravelMode, choice = "choice", shape = "long",
                 chid.var = "individual", alt.var = "mode", drop.index = TRUE)
head(TM)


###################################################
### code chunk number 13: mlogit.rnw:357-359
###################################################
data("Train", package="mlogit")
head(Train, 3)


###################################################
### code chunk number 14: mlogit.rnw:368-372
###################################################
Tr <- mlogit.data(Train, shape = 'wide', choice="choice", 
                  varying=4:11, sep="", alt.levels=c(1, 2), id = "id")
head(Tr, 3)
head(index(Tr), 3)


###################################################
### code chunk number 15: mlogit.rnw:456-457
###################################################
f <- mFormula(choice ~ vcost | income + size | travel)


###################################################
### code chunk number 16: mlogit.rnw:465-467
###################################################
f2 <- mFormula(choice ~ vcost + travel | income + size)
f2 <- mFormula(choice ~ vcost + travel | income + size | 0)


###################################################
### code chunk number 17: mlogit.rnw:470-472
###################################################
f3 <- mFormula(choice ~ 0 | income | 0)
f3 <- mFormula(choice ~ 0 | income)


###################################################
### code chunk number 18: mlogit.rnw:475-478
###################################################
f4 <- mFormula(choice ~ vcost + travel)
f4 <- mFormula(choice ~ vcost + travel | 1)
f4 <- mFormula(choice ~ vcost + travel | 1 | 0)


###################################################
### code chunk number 19: mlogit.rnw:484-488
###################################################
f5 <- mFormula(choice ~ vcost | 0 | travel)
f6 <- mFormula(choice ~ vcost | income + 0 | travel)
f6 <- mFormula(choice ~ vcost | income -1 | travel)
f7 <- mFormula(choice ~ 0 | income -1 | travel)


###################################################
### code chunk number 20: mlogit.rnw:495-497
###################################################
f <- mFormula(choice ~ vcost | income  | travel)
head(model.matrix(f, TM))


###################################################
### code chunk number 21: mlogit.rnw:1289-1293
###################################################
data("Train", package="mlogit")
Tr <- mlogit.data(Train, shape = 'wide', choice="choice", 
                  varying=4:11, sep="", alt.levels=c(1, 2), id = "id")



###################################################
### code chunk number 22: mlogit.rnw:1299-1301
###################################################
Tr$price <- Tr$price / 100 * 2.20371
Tr$time <- Tr$time / 60


###################################################
### code chunk number 23: mlogit.rnw:1308-1310
###################################################
ml.Train <- mlogit(choice~price+time+change+comfort | -1, Tr)
summary(ml.Train)


###################################################
### code chunk number 24: mlogit.rnw:1319-1320
###################################################
coef(ml.Train)[-1]/coef(ml.Train)[1]


###################################################
### code chunk number 25: mlogit.rnw:1331-1333
###################################################
Fish <- mlogit.data(Fishing, shape="wide", varying=2:9, choice="mode")
ml.Fish <- mlogit(mode~price | income | catch, Fish)


###################################################
### code chunk number 26: mlogit.rnw:1340-1342
###################################################
ml.Fish <- mlogit(mode~price | income | catch, Fishing, shape = "wide", varying = 2:9)
summary(ml.Fish)


###################################################
### code chunk number 27: mlogit.rnw:1349-1351
###################################################
head(fitted(ml.Fish))
head(fitted(ml.Fish, outcome=FALSE))


###################################################
### code chunk number 28: mlogit.rnw:1368-1370
###################################################
mlogit(mode~price | income | catch, Fish, reflevel='charter', 
       alt.subset=c('beach', 'pier', 'charter'))


###################################################
### code chunk number 29: mlogit.rnw:1417-1423
###################################################
data("Game", package = "mlogit")
data("Game2", package = "mlogit")
head(Game,2)
head(Game2, 7)
nrow(Game)
nrow(Game2)


###################################################
### code chunk number 30: mlogit.rnw:1433-1437
###################################################
G <- mlogit.data(Game2, shape="long", choice="ch", alt.var='platform', ranked=TRUE)
G <- mlogit.data(Game, shape="wide", choice="ch", varying=1:12, ranked=TRUE)
head(G)
nrow(G)


###################################################
### code chunk number 31: mlogit.rnw:1447-1448
###################################################
summary(mlogit(ch~own|hours+age, G, reflevel="PC"))


###################################################
### code chunk number 32: mlogit.rnw:1554-1555
###################################################
data("ModeCanada", package = "mlogit")


###################################################
### code chunk number 33: mlogit.rnw:1562-1569
###################################################
busUsers <- with(ModeCanada, case[choice == 1 & alt == 'bus'])
Bhat <- subset(ModeCanada, !case %in% busUsers & alt != 'bus' & nchoice == 4)
Bhat$alt <- Bhat$alt[drop = TRUE]
Bhat <- mlogit.data(Bhat, shape='long', chid.var = 'case',
                    alt.var = 'alt', choice='choice',
                    drop.index=TRUE)



###################################################
### code chunk number 34: mlogit.rnw:1574-1577
###################################################
ml.MC <- mlogit(choice ~ freq + cost + ivt + ovt | urban + income, Bhat, reflevel = 'car')
hl.MC <- mlogit(choice ~ freq + cost + ivt + ovt | urban + income, Bhat, reflevel = 'car', heterosc = TRUE)
summary(hl.MC)


###################################################
### code chunk number 35: mlogit.rnw:1589-1598
###################################################
data("TravelMode",package="AER")
TravelMode <- mlogit.data(TravelMode,choice="choice",shape="long",
                          alt.var="mode",chid.var="individual")
TravelMode$avinc <- with(TravelMode,(mode=='air')*income)
ml.TM <- mlogit(choice ~ wait + gcost + avinc, TravelMode, 
                reflevel = "car")
hl.TM <- mlogit(choice ~ wait + gcost + avinc, TravelMode, 
             reflevel = "car", heterosc = TRUE)
summary(hl.TM)


###################################################
### code chunk number 36: mlogit.rnw:1608-1609
###################################################
c(coef(hl.TM)[7:9], sp.car = 1)*pi/sqrt(6)


###################################################
### code chunk number 37: mlogit.rnw:1851-1859
###################################################
data("TravelMode",package="AER")
TravelMode <- mlogit.data(TravelMode,choice="choice",shape="long",
                          alt.var="mode",chid.var="individual")
TravelMode$avinc <- with(TravelMode,(mode=='air')*income)
nl.TM <- mlogit(choice ~ wait + gcost + avinc, TravelMode, reflevel = "car", 
                nests = list(fly = "air", ground = c("train", "bus", "car")), 
                unscaled=TRUE)
summary(nl.TM)


###################################################
### code chunk number 38: mlogit.rnw:1870-1873
###################################################
data("HC", package = "mlogit")
HC <- mlogit.data(HC, varying = c(2:8, 10:16), choice = "depvar", shape = "wide")
head(HC)


###################################################
### code chunk number 39: mlogit.rnw:1881-1883
###################################################
cooling.modes <- HC$alt %in% c("gcc", "ecc", "erc","hpc")
HC$icca[!cooling.modes] <- HC$occa[!cooling.modes] <- 0


###################################################
### code chunk number 40: mlogit.rnw:1888-1893
###################################################
ml.HC <- mlogit(depvar~occa+icca+och+ich, HC)
nl.HC <- mlogit(depvar~occa+icca+och+ich, HC, 
                nests = list(cooling = c('ecc', 'erc', 'gcc', 'hpc'), 
                  noncool = c('ec', 'gc', 'er')))
summary(nl.HC)


###################################################
### code chunk number 41: mlogit.rnw:1902-1903
###################################################
nl.HC.u <- update(nl.HC, un.nest.el = TRUE)


###################################################
### code chunk number 42: mlogit.rnw:2152-2155
###################################################
pcl <- mlogit(choice~freq+cost+ivt+ovt, Bhat, reflevel='car',
              nests='pcl', constPar=c('iv.train.air'))
summary(pcl)


###################################################
### code chunk number 43: mlogit.rnw:2482-2487
###################################################
data("Train", package = "mlogit")
Tr <- mlogit.data(Train, shape = "wide", varying = 4:11, 
                  choice = "choice", sep = "", 
                  opposite = c("price", "time", "change", "comfort"),
                  alt.levels=c("choice1", "choice2"), id="id")


###################################################
### code chunk number 44: mlogit.rnw:2489-2494 (eval = FALSE)
###################################################
## Train.ml <- mlogit(choice ~ price + time + change + comfort, Tr)
## Train.mxlc <- mlogit(choice ~ price + time + change + comfort, Tr,
##                panel = TRUE, rpar = c(time = "cn", change = "n", comfort = "ln"),
##                correlation = TRUE, R = 100, halton = NA)
## Train.mxlu <- update(Train.mxlc, correlation = FALSE)


###################################################
### code chunk number 45: mlogit.rnw:2497-2498
###################################################
data("Examples", package = "mlogit")


###################################################
### code chunk number 46: mlogit.rnw:2512-2517
###################################################
time.value <- rpar(Train.mxlc, "time", norm = "price")
summary(time.value)
med(time.value)
mean(time.value)
stdev(time.value)


###################################################
### code chunk number 47: mlogit.rnw:2521-2524
###################################################
cor.mlogit(Train.mxlc)
cov.mlogit(Train.mxlc)
stdev(Train.mxlc)


###################################################
### code chunk number 48: mlogit.rnw:2868-2871 (eval = FALSE)
###################################################
## data("Fishing", package = "mlogit")
## Fish <- mlogit.data(Fishing, shape="wide", varying=2:9, choice="mode")
## Fish.mprobit <- mlogit(mode~price | income | catch, Fish, probit = TRUE, alt.subset=c('beach', 'boat','pier'))


###################################################
### code chunk number 49: mlogit.rnw:2874-2875
###################################################
summary(Fish.mprobit)


###################################################
### code chunk number 50: mlogit.rnw:2985-2986
###################################################
ml.Fish <- mlogit(mode~price | income | catch, Fishing, shape = "wide", varying = 2:9)


###################################################
### code chunk number 51: mlogit.rnw:2994-2995
###################################################
ml.Fish.c <- update(ml.Fish, . ~ . | . - income | .)


###################################################
### code chunk number 52: mlogit.rnw:3001-3004
###################################################
waldtest(ml.Fish, ml.Fish.c)
lrtest(ml.Fish, ml.Fish.c)
scoretest(ml.Fish.c, ml.Fish)


###################################################
### code chunk number 53: mlogit.rnw:3009-3019 (eval = FALSE)
###################################################
## lrtest(ml.Fish, . ~ . | . - income | .)
## lrtest(ml.Fish, mode ~ price | 1 | catch)
## lrtest(ml.Fish.c, . ~ . | . + income | .)
## lrtest(ml.Fish.c, mode ~ price | income | catch)
## waldtest(ml.Fish, . ~ . | . - income | .)
## waldtest(ml.Fish, mode ~ price | 1 | catch)
## waldtest(ml.Fish.c, . ~ . | . + income | .)
## waldtest(ml.Fish.c, mode ~ price | income | catch)
## scoretest(ml.Fish.c, . ~ . | . + income | .)
## scoretest(ml.Fish.c, mode ~ price | income | catch)


###################################################
### code chunk number 54: mlogit.rnw:3034-3036
###################################################
lrtest(hl.MC, ml.MC)
waldtest(hl.MC, heterosc = FALSE)


###################################################
### code chunk number 55: mlogit.rnw:3039-3041
###################################################
lrtest(hl.MC)
waldtest(hl.MC)


###################################################
### code chunk number 56: mlogit.rnw:3045-3047
###################################################
library("car")
linearHypothesis(hl.MC, c('sp.air=1', 'sp.train=1'))


###################################################
### code chunk number 57: mlogit.rnw:3054-3055
###################################################
scoretest(ml.MC, heterosc = TRUE)


###################################################
### code chunk number 58: mlogit.rnw:3066-3069
###################################################
c(wald = waldtest(hl.TM)$statistic, 
  lr = lrtest(hl.TM)$Chisq[2],
  score = scoretest(ml.TM, heterosc = TRUE)$statistic)


###################################################
### code chunk number 59: mlogit.rnw:3092-3096
###################################################
lrtest(nl.HC)
waldtest(nl.HC)
scoretest(ml.HC, nests = list(cooling = c('ecc', 'erc', 'gcc', 'hpc'), 
                   noncool = c('ec', 'gc', 'er')))


###################################################
### code chunk number 60: mlogit.rnw:3101-3102
###################################################
linearHypothesis(nl.HC, c("iv.cooling=1", "iv.noncool=1"))


###################################################
### code chunk number 61: mlogit.rnw:3111-3115
###################################################
lrtest(nl.HC,nl.HC.u)
waldtest(nl.HC, un.nest.el = TRUE)
scoretest(nl.HC.u, un.nest.el = FALSE)
linearHypothesis(nl.HC, "iv.cooling=iv.noncool")


###################################################
### code chunk number 62: mlogit.rnw:3136-3140
###################################################
lrtest(Train.mxlu, Train.ml)
waldtest(Train.mxlu)
scoretest(Train.ml, rpar = c(time = "n", change = "n", comfort = "n"), R = 100,     
          correlation = FALSE, halton = NA, panel = TRUE)


###################################################
### code chunk number 63: mlogit.rnw:3145-3149
###################################################
lrtest(Train.mxlc, Train.ml)
waldtest(Train.mxlc)
scoretest(Train.ml, rpar = c(time = "n", change = "n", comfort = "n"), R = 100,     
          correlation = TRUE, halton = NA, panel = TRUE)


###################################################
### code chunk number 64: mlogit.rnw:3155-3158
###################################################
lrtest(Train.mxlc, Train.mxlu)
waldtest(Train.mxlc, correlation = FALSE)
scoretest(Train.mxlu, correlation = TRUE)


