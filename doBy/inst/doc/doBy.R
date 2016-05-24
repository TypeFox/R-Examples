### R code from vignette source 'doBy.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: doBy.Rnw:22-25
###################################################
require( doBy )
prettyVersion <- packageDescription("doBy")$Version
prettyDate <- format(Sys.Date())


###################################################
### code chunk number 2: doBy.Rnw:77-81
###################################################
dir.create("figures")
oopt <- options()
options("digits"=4, "width"=80, "prompt"=" ", "continue"="  ")
options(useFancyQuotes="UTF-8")


###################################################
### code chunk number 3: doBy.Rnw:104-105
###################################################
library(doBy)


###################################################
### code chunk number 4: doBy.Rnw:123-128
###################################################
data(CO2)
CO2 <- transform(CO2, Treat=Treatment, Treatment=NULL)
levels(CO2$Treat) <- c("nchil","chil")
levels(CO2$Type)  <- c("Que","Mis")
CO2 <- subset(CO2, Plant %in% c("Qn1", "Qc1", "Mn1", "Mc1"))


###################################################
### code chunk number 5: doBy.Rnw:140-141
###################################################
airquality <- subset(airquality, Month %in% c(5,6))


###################################################
### code chunk number 6: doBy.Rnw:169-171
###################################################
myfun1 <- function(x){c(m=mean(x), v=var(x))}
summaryBy( conc + uptake ~ Plant, data=CO2, FUN=myfun1)


###################################################
### code chunk number 7: doBy.Rnw:179-180
###################################################
summaryBy( list(c("conc","uptake"), "Plant"), data=CO2, FUN=myfun1)


###################################################
### code chunk number 8: doBy.Rnw:186-188
###################################################
myfun2 <- function(x){c(mean(x), var(x))}
summaryBy( conc + uptake ~ Plant, data=CO2, FUN=myfun2)


###################################################
### code chunk number 9: doBy.Rnw:196-197
###################################################
summaryBy( conc + uptake ~ Plant, data=CO2, FUN=list( mean, var ) )


###################################################
### code chunk number 10: doBy.Rnw:204-205
###################################################
summaryBy(uptake~Plant, data=CO2, FUN=list( mean, var, myfun1 ))


###################################################
### code chunk number 11: doBy.Rnw:210-212
###################################################
summaryBy(uptake~Plant, data=CO2, FUN=list( mean, var, myfun1 ),
          fun.names=c("mean","var","mm","vv"))


###################################################
### code chunk number 12: doBy.Rnw:222-224
###################################################
summaryBy(log(uptake) + I(conc+uptake) + conc+uptake ~ Plant, data=CO2,
          FUN=myfun1)


###################################################
### code chunk number 13: doBy.Rnw:230-232
###################################################
summaryBy(log(uptake) + I(conc+uptake) + conc + uptake ~ Plant, data=CO2,
          FUN=myfun1, var.names=c("log.upt", "conc+upt", "conc", "upt"))


###################################################
### code chunk number 14: doBy.Rnw:242-244
###################################################
summaryBy(log(uptake)+I(conc+uptake)~Plant, data=CO2, p2d=TRUE,
FUN=myfun1)


###################################################
### code chunk number 15: doBy.Rnw:258-260
###################################################
summaryBy(conc+uptake~Plant, data=CO2, FUN=myfun1, id=~Type+Treat)
summaryBy(conc+uptake~Plant, data=CO2, FUN=myfun1, id=c("Type","Treat"))


###################################################
### code chunk number 16: doBy.Rnw:279-280
###################################################
summaryBy(log(uptake)+I(conc+uptake)+. ~Plant, data=CO2, FUN=myfun1)


###################################################
### code chunk number 17: doBy.Rnw:291-292
###################################################
summaryBy(log(uptake) ~Plant+., data=CO2, FUN=myfun1)


###################################################
### code chunk number 18: doBy.Rnw:301-302
###################################################
summaryBy(log(uptake) ~ 1, data=CO2, FUN=myfun1)


###################################################
### code chunk number 19: doBy.Rnw:314-316
###################################################
summaryBy(conc+uptake+log(uptake)~Plant,
data=CO2, FUN=mean, id=~Type+Treat, keep.names=TRUE)


###################################################
### code chunk number 20: doBy.Rnw:329-330
###################################################
x<-orderBy(~Temp+Month, data=airquality)


###################################################
### code chunk number 21: doBy.Rnw:334-335
###################################################
head(x)


###################################################
### code chunk number 22: doBy.Rnw:341-343
###################################################
x<-orderBy(~-Temp+Month, data=airquality)
head(x)


###################################################
### code chunk number 23: doBy.Rnw:354-356
###################################################
x<-splitBy(~Month, data=airquality)
x


###################################################
### code chunk number 24: doBy.Rnw:362-363
###################################################
x[['5']]


###################################################
### code chunk number 25: doBy.Rnw:368-369
###################################################
attr(x,"groupid")


###################################################
### code chunk number 26: doBy.Rnw:379-380
###################################################
sampleBy(~1, frac=0.5, data=airquality)


###################################################
### code chunk number 27: doBy.Rnw:386-387
###################################################
sampleBy(~Month, frac=0.2, data=airquality,systematic=T)


###################################################
### code chunk number 28: doBy.Rnw:398-399
###################################################
subsetBy(~Month, subset=Wind>mean(Wind), data=airquality)


###################################################
### code chunk number 29: doBy.Rnw:412-414
###################################################
transformBy(~Month, data=airquality, minW=min(Wind), maxW=max(Wind),
    chg=sum(range(Wind)*c(-1,1)))


###################################################
### code chunk number 30: doBy.Rnw:428-433
###################################################
data(dietox)
dietox <- orderBy(~Pig+Time, data=dietox)
FEfun  <- function(d){c(NA, diff(d$Weight)/diff(d$Feed))}
v      <- lapplyBy(~Pig, data=dietox, FEfun)
dietox$FE <- unlist(v)


###################################################
### code chunk number 31: doBy.Rnw:438-442
###################################################
dietox <- orderBy(~Pig+Time, data=dietox)
wdata  <- splitBy(~Pig, data=dietox)
v      <- lapply(wdata, FEfun)
dietox$FE <- unlist(v)


###################################################
### code chunk number 32: doBy.Rnw:450-454
###################################################
x<-scaleBy( list(c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
                 "Species"),     data=iris)
head(x)
head(iris)


###################################################
### code chunk number 33: doBy.Rnw:463-467
###################################################
mydata <- data.frame(y=rnorm(32), x=rnorm(32),
g1=factor(rep(c(1,2),each=16)), g2=factor(rep(c(1,2), each=8)),
g3=factor(rep(c(1,2),each=4)))
head(mydata)


###################################################
### code chunk number 34: doBy.Rnw:472-480
###################################################
## Based on the formula interface to t.test
t.testBy1 <- function(formula, group, data, ...){
  formulaFunBy(formula, group, data, FUN=t.test, class="t.testBy1", ...)
}
## Based on the default interface to t.test
t.testBy2 <- function(formula, group, data, ...){
  xyFunBy(formula, group, data, FUN=t.test, class="t.testBy1", ...)
}


###################################################
### code chunk number 35: doBy.Rnw:487-489
###################################################
t.testBy1(y~g1, ~g2, data=mydata)
t.testBy2(y~x,  ~g2, data=mydata)


###################################################
### code chunk number 36: doBy.Rnw:508-514
###################################################
ff  <- function(a,b=2,c=4){a+b+c}
ff1 <- specialize(ff, arglist=list(a=1, b=7, yy=123))
ff1
gg  <- rnorm
gg1 <- specialize(gg, list(n=10))
gg1


###################################################
### code chunk number 37: doBy.Rnw:519-522
###################################################
f  <- function(a) {a <- a + 1; a}
f1 <- specialize(f, list(a = 10))
f1


###################################################
### code chunk number 38: doBy.Rnw:532-535
###################################################
x <- c(1,1,1,2,2,2,1,1,1,3)
firstobs(x)
lastobs(x)


###################################################
### code chunk number 39: doBy.Rnw:540-542
###################################################
firstobs(~Plant, data=CO2)
lastobs(~Plant, data=CO2)


###################################################
### code chunk number 40: doBy.Rnw:551-554
###################################################
x <- c(1:4,0:5,11,NA,NA)
which.maxn(x,3)
which.minn(x,5)


###################################################
### code chunk number 41: doBy.Rnw:563-568
###################################################
x <- c(1,1,2,2,2,1,1,3,3,3,3,1,1,1)
subSeq(x)
subSeq(x, item=1)
subSeq(letters[x])
subSeq(letters[x],item="a")


###################################################
### code chunk number 42: doBy.Rnw:576-580
###################################################
x <- c("dec","jan","feb","mar","apr","may")
src1 <- list(c("dec","jan","feb"), c("mar","apr","may"))
tgt1 <- list("winter","spring")
recodeVar(x,src=src1,tgt=tgt1)


###################################################
### code chunk number 43: doBy.Rnw:587-589
###################################################
head(renameCol(CO2, 1:2, c("kk","ll")))
head(renameCol(CO2, c("Plant","Type"), c("kk","ll")))


###################################################
### code chunk number 44: doBy.Rnw:597-598
###################################################
yvar <- c(0,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0)


###################################################
### code chunk number 45: doBy.Rnw:605-606
###################################################
tvar <- seq_along(yvar) + c(0.1,0.2)


###################################################
### code chunk number 46: doBy.Rnw:611-612
###################################################
tse<- timeSinceEvent(yvar,tvar)


###################################################
### code chunk number 47: doBy.Rnw:627-632
###################################################
plot(sign.tse~tvar, data=tse, type="b")
grid()
rug(tse$tvar[tse$yvar==1], col='blue',lwd=4)
points(scale(tse$run), col=tse$run, lwd=2)
lines(abs.tse+.2~tvar, data=tse, type="b",col=3)


###################################################
### code chunk number 48: doBy.Rnw:636-641
###################################################
plot(tae~tvar, data=tse, ylim=c(-6,6),type="b")
grid()
lines(tbe~tvar, data=tse, type="b", col='red')
rug(tse$tvar[tse$yvar==1], col='blue',lwd=4)
lines(run~tvar, data=tse, col='cyan',lwd=2)


###################################################
### code chunk number 49: doBy.Rnw:645-649
###################################################
plot(ewin~tvar, data=tse,ylim=c(1,4))
rug(tse$tvar[tse$yvar==1], col='blue',lwd=4)
grid()
lines(run~tvar, data=tse,col='red')


###################################################
### code chunk number 50: doBy.Rnw:655-656
###################################################
tse$tvar[tse$abs<=1]


###################################################
### code chunk number 51: doBy.Rnw:663-666
###################################################
lynx <- as.numeric(lynx)
tvar <- 1821:1934
plot(tvar,lynx,type='l')


###################################################
### code chunk number 52: doBy.Rnw:672-676
###################################################
yyy <- lynx>mean(lynx)
head(yyy)
sss <- subSeq(yyy,TRUE)
sss


###################################################
### code chunk number 53: doBy.Rnw:680-682
###################################################
plot(tvar,lynx,type='l')
rug(tvar[sss$midpoint],col='blue',lwd=4)


###################################################
### code chunk number 54: doBy.Rnw:687-690
###################################################
yvar <- rep(0,length(lynx))
yvar[sss$midpoint] <- 1
str(yvar)


###################################################
### code chunk number 55: doBy.Rnw:694-696
###################################################
tse <- timeSinceEvent(yvar,tvar)
head(tse,20)


###################################################
### code chunk number 56: doBy.Rnw:702-705
###################################################
len1 <- tapply(tse$ewin, tse$ewin, length)
len2 <- tapply(tse$run, tse$run, length)
c(median(len1),median(len2),mean(len1),mean(len2))


###################################################
### code chunk number 57: doBy.Rnw:710-713
###################################################
tse$lynx <- lynx
tse2 <- na.omit(tse)
plot(lynx~tae, data=tse2)


###################################################
### code chunk number 58: doBy.Rnw:717-720
###################################################
plot(tvar,lynx,type='l',lty=2)
mm <- lm(lynx~tae+I(tae^2)+I(tae^3), data=tse2)
lines(fitted(mm)~tvar, data=tse2, col='red')


###################################################
### code chunk number 59: doBy.Rnw:733-734
###################################################
options(oopt)


###################################################
### code chunk number 60: doBy.Rnw:749-750
###################################################
CO2


###################################################
### code chunk number 61: doBy.Rnw:755-756
###################################################
head(airquality, n=20)


