library(ISwR)
.make.epsf <- Sys.getenv("EPSF")=="y"
ps.options(height=3.5, width=4.4, pointsize=8, horiz=F)
if (.make.epsf) X11(height=3.5,width=4.4,pointsize=8) else postscript()
dev.copy2eps <- function(...) invisible(grDevices::dev.copy2eps(...))
par(mar=c(4,4,3,2)+.1)
options(width=66, useFancyQuotes="TeX")
suppressWarnings(RNGversion("1.5.1")) #Yes, Kinderman-Ramage was buggy...
set.seed(310367)
#Rprof(interval=.001)
plot(rnorm(500))
2 + 2
exp(-2)
rnorm(15)
x <- 2
x
x + x
weight <- c(60, 72, 57, 90, 95, 72)
weight
height <- c(1.75, 1.80, 1.65, 1.90, 1.74, 1.91)
bmi <- weight/height^2
bmi
sum(weight)
sum(weight)/length(weight)
xbar <- sum(weight)/length(weight)
weight - xbar
(weight - xbar)^2
sum((weight - xbar)^2)
sqrt(sum((weight - xbar)^2)/(length(weight) - 1))
mean(weight)
sd(weight)
t.test(bmi, mu=22.5)
plot(height,weight)
if (.make.epsf) dev.copy2eps(file="h-w.ps")
plot(height, weight, pch=2)
if (.make.epsf) dev.copy2eps(file="h-w-triangle.ps")
hh <- c(1.65, 1.70, 1.75, 1.80, 1.85, 1.90)
lines(hh, 22.5 * hh^2)
if (.make.epsf) dev.copy2eps(file="h-w-line.ps")
args(plot.default)
c("Huey","Dewey","Louie")
c('Huey','Dewey','Louie')
c(T,T,F,T)
bmi > 25
cat(c("Huey","Dewey","Louie"))
cat("Huey","Dewey","Louie", "\n")
cat("What is \"R\"?\n")
c(42,57,12,39,1,3,4)
x <- c(1, 2, 3)
y <- c(10, 20)
c(x, y, 5)
x <- c(red="Huey", blue="Dewey", green="Louie")
x
names(x)
c(FALSE, 3)
c(pi, "abc")
c(FALSE, "abc")
seq(4,9)
seq(4,10,2)
4:9
oops <- c(7,9,13)
rep(oops,3)
rep(oops,1:3)
rep(1:2,c(10,15))
x <- 1:12
dim(x) <- c(3,4)
x
matrix(1:12,nrow=3,byrow=T)
x <- matrix(1:12,nrow=3,byrow=T)
rownames(x) <- LETTERS[1:3]
x
t(x)
cbind(A=1:4,B=5:8,C=9:12)
rbind(A=1:4,B=5:8,C=9:12)
pain <- c(0,3,2,2,1)
fpain <- factor(pain,levels=0:3)
levels(fpain) <- c("none","mild","medium","severe")
fpain
as.numeric(fpain)
levels(fpain)
intake.pre <- c(5260,5470,5640,6180,6390,
6515,6805,7515,7515,8230,8770)
intake.post <- c(3910,4220,3885,5160,5645,
4680,5265,5975,6790,6900,7335)
mylist <- list(before=intake.pre,after=intake.post)
mylist
mylist$before
d <- data.frame(intake.pre,intake.post)
d
d$intake.pre
intake.pre[5]
intake.pre[c(3,5,7)]
v <- c(3,5,7)
intake.pre[v]
intake.pre[1:5]
intake.pre[-c(3,5,7)]
intake.post[intake.pre > 7000]
intake.post[intake.pre > 7000 & intake.pre <= 8000]
intake.pre > 7000 & intake.pre <= 8000
d <- data.frame(intake.pre,intake.post)
d[5,1]
d[5,]
d[d$intake.pre>7000,]
sel <- d$intake.pre>7000
sel
d[sel,]
d[1:2,]
head(d)
energy
exp.lean <- energy$expend[energy$stature=="lean"]
exp.obese <- energy$expend[energy$stature=="obese"]
l <- split(energy$expend, energy$stature)
l
lapply(thuesen, mean, na.rm=T)
sapply(thuesen, mean, na.rm=T)
replicate(10,mean(rexp(20)))
m <- matrix(rnorm(12),4)
m
apply(m, 2, min)
tapply(energy$expend, energy$stature, median)
intake$post
sort(intake$post)
order(intake$post)
o <- order(intake$post)
intake$post[o]
intake$pre[o]
intake.sorted <- intake[o,]
save.image("ch1.RData")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
load("ch1.RData")
.foo <- dev.copy2eps
rm(dev.copy2eps)
ls()
dev.copy2eps <- .foo
rm(height, weight)
sink("myfile")
ls()
sink()
attach(thuesen)
blood.glucose
search()
detach()
search()
thue2 <- subset(thuesen,blood.glucose<7)
thue2
thue3 <- transform(thuesen,log.gluc=log(blood.glucose))
thue3
thue4 <- within(thuesen,{
   log.gluc <- log(blood.glucose)
   m <- mean(log.gluc)
   centered.log.gluc <- log.gluc - m
   rm(m)
})
thue4
d <- par(mar=c(5,4,4,2)+.1)
x <- runif(50,0,2)
y <- runif(50,0,2)
plot(x, y, main="Main title", sub="subtitle",
     xlab="x-label", ylab="y-label")
text(0.6,0.6,"text at (0.6,0.6)")
abline(h=.6,v=.6)
for (side in 1:4) mtext(-1:4,side=side,at=.7,line=-1:4)
mtext(paste("side",1:4), side=1:4, line=-1,font=2)
if (.make.epsf) dev.copy2eps(file="layout.ps")
par(d)
plot(x, y, type="n", xlab="", ylab="", axes=F)
points(x,y)
axis(1)
axis(2,at=seq(0.2,1.8,0.2))
box()
title(main="Main title", sub="subtitle",
    xlab="x-label", ylab="y-label")
set.seed(1234) #make it happen....
x <- rnorm(100)
hist(x,freq=F)
curve(dnorm(x),add=T)
h <- hist(x, plot=F)
ylim <- range(0, h$density, dnorm(0))
hist(x, freq=F, ylim=ylim)
curve(dnorm(x), add=T)
if (.make.epsf) dev.copy2eps(file="hist+norm.ps")
hist.with.normal <- function(x, xlab=deparse(substitute(x)),...)
{
    h <- hist(x, plot=F, ...)
    s <- sd(x)
    m <- mean(x)
    ylim <- range(0,h$density,dnorm(0,sd=s))
    hist(x, freq=F, ylim=ylim, xlab=xlab, ...)
    curve(dnorm(x,m,s), add=T)
}
hist.with.normal(rnorm(200))
y <- 12345
x <- y/2
while (abs(x*x-y) > 1e-10) x <- (x + y/x)/2
x
x^2
x <- y/2
repeat{
    x <- (x + y/x)/2
    if (abs(x*x-y) < 1e-10) break
}
x
x <- seq(0, 1,.05)
plot(x, x, ylab="y", type="l")
for ( j in 2:8 ) lines(x, x^j)
t.test(bmi, mu=22.5)$p.value
print
length(methods("print")) # quoted in text
thuesen2 <- read.table(
   system.file("rawdata","thuesen.txt",package="ISwR"), header=T)
thuesen2
levels(secretin$time)
system.file("rawdata", "thuesen.txt", package="ISwR")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
sample(1:40,5)
sample(c("H","T"), 10, replace=T)
sample(c("succ", "fail"), 10, replace=T, prob=c(0.9, 0.1))
1/prod(40:36)
prod(5:1)/prod(40:36)
1/choose(40,5)
x <- seq(-4,4,0.1)
plot(x,dnorm(x),type="l")
if (.make.epsf) dev.copy2eps(file="bellcurve.ps")
x <- 0:50
plot(x,dbinom(x,size=50,prob=.33),type="h")
if (.make.epsf) dev.copy2eps(file="binomdist.ps")
1-pnorm(160,mean=132,sd=13)
pbinom(16,size=20,prob=.5)
1-pbinom(15,size=20,prob=.5)
1-pbinom(15,20,.5)+pbinom(4,20,.5)
xbar <- 83
sigma <- 12
n <- 5
sem <- sigma/sqrt(n)
sem
xbar + sem * qnorm(0.025)
xbar + sem * qnorm(0.975)
set.seed(310367)
rnorm(10)
rnorm(10)
rnorm(10,mean=7,sd=5)
rbinom(10,size=20,prob=.5)
 ## no data sets used by exercises
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
x <- rnorm(50)
mean(x)
sd(x)
var(x)
median(x)
quantile(x)
pvec <- seq(0,1,0.1)
pvec
quantile(x,pvec)
attach(juul)
 mean(igf1)
mean(igf1,na.rm=T)
sum(!is.na(igf1))
summary(igf1)
summary(juul)
detach(juul)
juul$sex <- factor(juul$sex,labels=c("M","F"))
juul$menarche <- factor(juul$menarche,labels=c("No","Yes"))
juul$tanner <- factor(juul$tanner,
                      labels=c("I","II","III","IV","V"))
attach(juul)
summary(juul)
hist(x)
if (.make.epsf) dev.copy2eps(file="hist.ps")
mid.age <- c(2.5,7.5,13,16.5,17.5,19,22.5,44.5,70.5)
acc.count <- c(28,46,58,20,31,64,149,316,103)
age.acc <- rep(mid.age,acc.count)
brk <- c(0,5,10,16,17,18,20,25,60,80)
hist(age.acc,breaks=brk)
if (.make.epsf) dev.copy2eps(file="hist-acc-right.ps")
n <- length(x)
plot(sort(x),(1:n)/n,type="s",ylim=c(0,1))
if (.make.epsf) dev.copy2eps(file="empdist.ps")
qqnorm(x)
if (.make.epsf) dev.copy2eps(file="qqnorm.ps")
par(mfrow=c(1,2))
boxplot(IgM)
boxplot(log(IgM))
par(mfrow=c(1,1))
if (.make.epsf) dev.copy2eps(file="boxplot-IgM.ps")
attach(red.cell.folate)
tapply(folate,ventilation,mean)
tapply(folate,ventilation,sd)
tapply(folate,ventilation,length)
xbar <- tapply(folate, ventilation, mean)
s <- tapply(folate, ventilation, sd)
n <- tapply(folate, ventilation, length)
cbind(mean=xbar, std.dev=s, n=n)
tapply(igf1, tanner, mean)
tapply(igf1, tanner, mean, na.rm=T)
aggregate(juul[c("age","igf1")],
          list(sex=juul$sex), mean, na.rm=T)
aggregate(juul[c("age","igf1")], juul["sex"], mean, na.rm=T)
by(juul, juul["sex"], summary)
attach(energy)
expend.lean <- expend[stature=="lean"]
expend.obese <- expend[stature=="obese"]
par(mfrow=c(2,1))
hist(expend.lean,breaks=10,xlim=c(5,13),ylim=c(0,4),col="white")
hist(expend.obese,breaks=10,xlim=c(5,13),ylim=c(0,4),col="grey")
par(mfrow=c(1,1))
if (.make.epsf) dev.copy2eps(file="expend-hist-2on1.ps")
boxplot(expend ~ stature)
if (.make.epsf) dev.copy2eps(file="boxplots-expend-stat.ps")
boxplot(expend.lean,expend.obese)
opar <- par(mfrow=c(2,2), mex=0.8, mar=c(3,3,2,1)+.1)
stripchart(expend ~ stature)
stripchart(expend ~ stature, method="stack")
stripchart(expend ~ stature, method="jitter")
stripchart(expend ~ stature, method="jitter", jitter=.03)
par(opar)
if (.make.epsf) dev.copy2eps(file="stripcharts-expend-stat.ps")
caff.marital <- matrix(c(652,1537,598,242,36,46,38,21,218
,327,106,67),
nrow=3,byrow=T)
caff.marital
colnames(caff.marital) <- c("0","1-150","151-300",">300")
rownames(caff.marital) <- c("Married","Prev.married","Single")
caff.marital
names(dimnames(caff.marital)) <- c("marital","consumption")
caff.marital
as.data.frame(as.table(caff.marital))
table(sex)
table(sex,menarche)
table(menarche,tanner)
xtabs(~ tanner + sex, data=juul)
xtabs(~ dgn + diab + coma, data=stroke)
ftable(coma + diab ~ dgn, data=stroke)
t(caff.marital)
tanner.sex <- table(tanner,sex)
tanner.sex
margin.table(tanner.sex,1)
margin.table(tanner.sex,2)
prop.table(tanner.sex,1)
tanner.sex/sum(tanner.sex)
total.caff <- margin.table(caff.marital,2)
total.caff
barplot(total.caff, col="white")
if (.make.epsf) dev.copy2eps(file="simple-bar.ps")
par(mfrow=c(2,2))
barplot(caff.marital, col="white")
barplot(t(caff.marital), col="white")
barplot(t(caff.marital), col="white", beside=T)
barplot(prop.table(t(caff.marital),2), col="white", beside=T)
par(mfrow=c(1,1))
if (.make.epsf) dev.copy2eps(file="mat-4-bar.ps")
barplot(prop.table(t(caff.marital),2),beside=T,
legend.text=colnames(caff.marital),
col=c("white","grey80","grey50","black"))
if (.make.epsf) dev.copy2eps(file="pretty-bar.ps")
dotchart(t(caff.marital), lcolor="black")
if (.make.epsf) dev.copy2eps(file="dotchart.ps")
opar <- par(mfrow=c(2,2),mex=0.8, mar=c(1,1,2,1))
slices <- c("white","grey80","grey50","black")
pie(caff.marital["Married",], main="Married", col=slices)
pie(caff.marital["Prev.married",],
         main="Previously married", col=slices)
pie(caff.marital["Single",], main="Single", col=slices)
par(opar)
if (.make.epsf) dev.copy2eps(file="pie.ps")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
daily.intake <- c(5260,5470,5640,6180,6390,6515,
6805,7515,7515,8230,8770)
mean(daily.intake)
sd(daily.intake)
quantile(daily.intake)
t.test(daily.intake,mu=7725)
t.test(daily.intake,mu=7725)
wilcox.test(daily.intake, mu=7725)
attach(energy)
energy
t.test(expend~stature)
t.test(expend~stature, var.equal=T)
var.test(expend~stature)
wilcox.test(expend~stature)
attach(intake)
intake
post - pre
t.test(pre, post, paired=T)
t.test(pre, post) #WRONG!
wilcox.test(pre, post, paired=T)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
attach(thuesen)
lm(short.velocity~blood.glucose)
summary(lm(short.velocity~blood.glucose))
summary(lm(short.velocity~blood.glucose))
plot(blood.glucose,short.velocity)
abline(lm(short.velocity~blood.glucose))
if (.make.epsf) dev.copy2eps(file="velo-gluc-line.ps")
lm.velo <- lm(short.velocity~blood.glucose)
fitted(lm.velo)
resid(lm.velo)
options(error=expression(NULL))
plot(blood.glucose,short.velocity)
lines(blood.glucose,fitted(lm.velo))
options(error=NULL)
lines(blood.glucose[!is.na(short.velocity)],fitted(lm.velo))
cc <- complete.cases(thuesen)
options(na.action=na.exclude)
lm.velo <- lm(short.velocity~blood.glucose)
fitted(lm.velo)
segments(blood.glucose,fitted(lm.velo),
         blood.glucose,short.velocity)
if (.make.epsf) dev.copy2eps(file="velo-gluc-seg.ps")
plot(fitted(lm.velo),resid(lm.velo))
if (.make.epsf) dev.copy2eps(file="velo-gluc-resid.ps")
qqnorm(resid(lm.velo))
if (.make.epsf) dev.copy2eps(file="velo-gluc-qqnorm.ps")
predict(lm.velo)
predict(lm.velo,int="c")
predict(lm.velo,int="p")
pred.frame <- data.frame(blood.glucose=4:20)
pp <- predict(lm.velo, int="p", newdata=pred.frame)
pc <- predict(lm.velo, int="c", newdata=pred.frame)
plot(blood.glucose,short.velocity,
     ylim=range(short.velocity, pp, na.rm=T))
pred.gluc <- pred.frame$blood.glucose
matlines(pred.gluc, pc, lty=c(1,2,2), col="black")
matlines(pred.gluc, pp, lty=c(1,3,3), col="black")
if (.make.epsf) dev.copy2eps(file="velo-gluc-final.ps")
options(error=expression(NULL))
cor(blood.glucose,short.velocity)
options(error=NULL)
cor(blood.glucose,short.velocity,use="complete.obs")
cor(thuesen,use="complete.obs")
cor.test(blood.glucose,short.velocity)
cor.test(blood.glucose,short.velocity,method="spearman")
cor.test(blood.glucose,short.velocity,method="kendall")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
attach(red.cell.folate)
summary(red.cell.folate)
anova(lm(folate~ventilation))
attach(juul)
anova(lm(igf1~tanner))  ## WRONG!
juul$tanner <- factor(juul$tanner,
                      labels=c("I","II","III","IV","V"))
detach(juul)
attach(juul)
summary(tanner)
anova(lm(igf1~tanner))
summary(lm(folate~ventilation))
pairwise.t.test(folate, ventilation, p.adj="bonferroni")
pairwise.t.test(folate,ventilation)
oneway.test(folate~ventilation)
 pairwise.t.test(folate,ventilation,pool.sd=F)
xbar <- tapply(folate, ventilation, mean)
s <- tapply(folate, ventilation, sd)
n <- tapply(folate, ventilation, length)
sem <- s/sqrt(n)
stripchart(folate~ventilation, method="jitter",
   jitter=0.05, pch=16, vert=T)
arrows(1:3,xbar+sem,1:3,xbar-sem,angle=90,code=3,length=.1)
lines(1:3,xbar,pch=4,type="b",cex=2)
if (.make.epsf) dev.copy2eps(file="oneway.ps")
bartlett.test(folate~ventilation)
kruskal.test(folate~ventilation)
attach(heart.rate)
heart.rate
gl(9,1,36)
gl(4,9,36,labels=c(0,30,60,120))
anova(lm(hr~subj+time))
interaction.plot(time, subj, hr)
if (.make.epsf) dev.copy2eps(file="interaction-plot.ps")
interaction.plot(ordered(time),subj,hr)
friedman.test(hr~time|subj,data=heart.rate)
attach(thuesen)
lm.velo <- lm(short.velocity~blood.glucose)
anova(lm.velo)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
prop.test(39,215,.15)
binom.test(39,215,.15)
lewitt.machin.success <- c(9,4)
lewitt.machin.total <- c(12,13)
prop.test(lewitt.machin.success,lewitt.machin.total)
matrix(c(9,4,3,9),2)
lewitt.machin <- matrix(c(9,4,3,9),2)
fisher.test(lewitt.machin)
chisq.test(lewitt.machin)
caesar.shoe
caesar.shoe.yes <- caesar.shoe["Yes",]
caesar.shoe.total <- margin.table(caesar.shoe,2)
caesar.shoe.yes
caesar.shoe.total
prop.test(caesar.shoe.yes,caesar.shoe.total)
prop.trend.test(caesar.shoe.yes,caesar.shoe.total)
caff.marital <- matrix(c(652,1537,598,242,36,46,38,21,218
,327,106,67),
nrow=3,byrow=T)
colnames(caff.marital) <- c("0","1-150","151-300",">300")
rownames(caff.marital) <- c("Married","Prev.married","Single")
caff.marital
chisq.test(caff.marital)
chisq.test(caff.marital)$expected
chisq.test(caff.marital)$observed
E <- chisq.test(caff.marital)$expected
O <- chisq.test(caff.marital)$observed
(O-E)^2/E
attach(juul)
chisq.test(tanner,sex)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
curve(pt(x,25,ncp=3), from=0, to=6)
abline(v=qt(.975,25))
if (.make.epsf) dev.copy2eps(file="noncentral-t.ps")
pt(qt(.975,25),25,ncp=3)
power.t.test(delta=0.5, sd=2, sig.level = 0.01, power=0.9)
power.t.test(n=450, delta=0.5, sd=2, sig.level = 0.01)
power.t.test(delta=0.5, sd=2, sig.level = 0.01, power=0.9,
alt="one.sided")
power.t.test(delta=10, sd=10*sqrt(2), power=0.85, type="paired")
power.prop.test(power=.85,p1=.15,p2=.30)
### no data sets in exercises
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
age <- subset(juul, age >= 10 & age <= 16)$age
range(age)
agegr <- cut(age, seq(10,16,2), right=F, include.lowest=T)
length(age)
table(agegr)
agegr2 <- cut(age, seq(10,16,2), right=F)
table(agegr2)
q <- quantile(age, c(0, .25, .50, .75, 1))
q
ageQ <- cut(age, q, include.lowest=T)
table(ageQ)
levels(ageQ) <- c("1st", "2nd", "3rd", "4th")
levels(agegr) <- c("10-11", "12-13", "14-15")
pain <- c(0,3,2,2,1)
fpain <- factor(pain,levels=0:3,
       labels=c("none","mild","medium","severe"))
text.pain <-  c("none","severe", "medium", "medium", "mild")
factor(text.pain)
ftpain <- factor(text.pain)
ftpain2 <- factor(ftpain,
                  levels=c("none", "mild", "medium", "severe"))
ftpain3 <- ftpain2
levels(ftpain3) <- list(
        none="none",
        intermediate=c("mild","medium"),
        severe="severe")
ftpain3
ftpain4 <- ftpain2
levels(ftpain4) <- c("none","intermediate","intermediate","severe")
ftpain4
stroke <- read.csv2(
  system.file("rawdata","stroke.csv", package="ISwR"),
  na.strings=".")
names(stroke) <- tolower(names(stroke))
head(stroke)
stroke <- transform(stroke,
   died = as.Date(died, format="%d.%m.%Y"),
   dstr = as.Date(dstr, format="%d.%m.%Y"))
summary(stroke$died)
summary(stroke$dstr)
summary(stroke$died - stroke$dstr)
head(stroke$died - stroke$dstr)
o <- options(width=60) # minor cheat for visual purposes
stroke <- transform(stroke,
  end = pmin(died,  as.Date("1996-1-1"), na.rm = T),
  dead = !is.na(died) & died < as.Date("1996-1-1"))
head(stroke)
options(o); rm(o)
stroke <- transform(stroke,
  obstime = as.numeric(end - dstr, units="days")/365.25)
rawstroke <- read.csv2(
  system.file("rawdata","stroke.csv", package="ISwR"),
  na.strings=".")
ix <- c("DSTR", "DIED")
rawstroke[ix] <- lapply(rawstroke[ix],
                        as.Date, format="%d.%m.%Y")
head(rawstroke)
ix <- 6:9
rawstroke[ix] <- lapply(rawstroke[ix],
                        factor, levels=0:1, labels=c("No","Yes"))
strokesub <- ISwR::stroke[1:10,2:3]
strokesub
strokesub <- transform(strokesub,
  event = !is.na(died))
strokesub <- transform(strokesub,
 obstime = ifelse(event, died-dstr, as.Date("1996-1-1") - dstr))
strokesub
juulgrl <- subset(juul, sex==2, select=-c(testvol,sex))
juulboy <- subset(juul, sex==1, select=-c(menarche,sex))
juulgrl$sex <- factor("F")
juulgrl$testvol <- NA
juulboy$sex <- factor("M")
juulboy$menarche <- NA
juulall <- rbind(juulboy, juulgrl)
names(juulall)
levels(juulall$sex)
head(nickel)
head(ewrates)
nickel <- transform(nickel,
  agr = trunc(agein/5)*5,
  ygr = trunc((dob+agein-1)/5)*5+1)
mrg <- merge(nickel, ewrates,
  by.x=c("agr","ygr"), by.y=c("age","year"))
head(mrg,10)
head(alkfos)
a2 <- alkfos
names(a2) <- sub("c", "c.", names(a2))
names(a2)
a.long <- reshape(a2, varying=2:8, direction="long")
head(a.long)
tail(a.long)
o <- with(a.long, order(id, time))
head(a.long[o,], 10)
a.long2 <- na.omit(a.long)
attr(a.long2, "reshapeLong") <- NULL
a.wide2 <- reshape(a.long2, direction="wide", v.names="c",
                 idvar="id", timevar="time")
head(a.wide2)
l <- split(a.long$c, a.long$id)
l[1:3]
l2 <- lapply(l, function(x) x / x[1])
a.long$c.adj <- unsplit(l2, a.long$id)
subset(a.long, id==1)
a.long$c.adj <- ave(a.long$c, a.long$id,
    FUN = function(x) x / x[1])
all.equal(unsplit(l2, a.long$id),  a.long$c.adj)
l <- split(a.long, a.long$id)
l2 <- lapply(l, transform, c.adj = c / c[1])
a.long2 <- unsplit(l2, a.long$id)
all.equal(a.long2$c.adj,  a.long$c.adj)
head(nickel)
entry <- pmax(nickel$agein, 60)
exit <- pmin(nickel$ageout, 65)
valid <- (entry < exit)
entry <- entry[valid]
exit  <- exit[valid]
cens <- (nickel$ageout[valid] > 65)
nickel60 <- nickel[valid,]
nickel60$icd[cens] <- 0
nickel60$agein <- entry
nickel60$ageout <- exit
nickel60$agr <- 60
nickel60$ygr <- with(nickel60, trunc((dob+agein-1)/5)*5+1)
head(nickel60)
trim <- function(start)
{
  end   <- start + 5
  entry <- pmax(nickel$agein, start)
  exit  <- pmin(nickel$ageout, end)
  valid <- (entry < exit)
  cens  <- (nickel$ageout[valid] > end)
  result <- nickel[valid,]
  result$icd[cens] <- 0
  result$agein <- entry[valid]
  result$ageout <- exit[valid]
  result$agr <- start
  result$ygr <- with(result, trunc((dob+agein-1)/5)*5+1)
  result
}
head(trim(60))
nickel.expand <- do.call("rbind", lapply(seq(20,95,5), trim))
head(nickel.expand)
subset(nickel.expand, id==4)
nickel.expand <- merge(nickel.expand, ewrates,
  by.x=c("agr","ygr"), by.y=c("age","year"))
head(nickel.expand)
all.equal(nickel.expand, ISwR::nickel.expand)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
par(mex=0.5)
pairs(cystfibr, gap=0, cex.labels=0.9)
if (.make.epsf) dev.copy2eps(file="cyst-fibr.ps",height=4.5,width=4.49)
attach(cystfibr)
if (exists("age",.GlobalEnv,inh=F)) rm(age)
if (exists("height",.GlobalEnv,inh=F)) rm(height)
if (exists("weight",.GlobalEnv,inh=F)) rm(weight)
summary(lm(pemax~age+sex+height+weight+bmp+fev1+rv+frc+tlc))
1-25.5^2/var(pemax)
anova(lm(pemax~age+sex+height+weight+bmp+fev1+rv+frc+tlc))
955.4+155.0+632.3+2862.2+1549.1+561.9+194.6+92.4
7002.9/8
875.36/648.7
1-pf(1.349407,8,15)
## Not command output:
m1<-lm(pemax~age+sex+height+weight+bmp+fev1+rv+frc+tlc)
m2<-lm(pemax~age)
anova(m1,m2)
summary(lm(pemax~age+sex+height+weight+bmp+fev1+rv+frc+tlc))
summary(lm(pemax~age+sex+height+weight+bmp+fev1+rv+frc))
summary(lm(pemax~age+sex+height+weight+bmp+fev1+rv))
summary(lm(pemax~age+sex+height+weight+bmp+fev1))
summary(lm(pemax~age+sex+height+weight+bmp))
summary(lm(pemax~age+height+weight+bmp))
summary(lm(pemax~height+weight+bmp))
summary(lm(pemax~weight+bmp))
summary(lm(pemax~weight))
summary(lm(pemax~age+weight+height))
summary(lm(pemax~age+height))
summary(lm(pemax~age))
summary(lm(pemax~height))
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
attach(cystfibr)
summary(lm(pemax~height+I(height^2)))
pred.frame <- data.frame(height=seq(110,180,2))
lm.pemax.hq <- lm(pemax~height+I(height^2))
predict(lm.pemax.hq,interval="pred",newdata=pred.frame)
pp <- predict(lm.pemax.hq,newdata=pred.frame,interval="pred")
pc <- predict(lm.pemax.hq,newdata=pred.frame,interval="conf")
plot(height,pemax,ylim=c(0,200))
matlines(pred.frame$height,pp,lty=c(1,2,2),col="black")
matlines(pred.frame$height,pc,lty=c(1,3,3),col="black")
if (.make.epsf) dev.copy2eps(file="pemax-height-quad.ps")
x <- runif(20)
y <- 2*x+rnorm(20,0,0.3)
summary(lm(y~x))
summary(lm(y~x-1))
anova(lm(y~x))
anova(lm(y~x-1))
model.matrix(pemax~height+weight)
attach(red.cell.folate)
model.matrix(folate~ventilation)
attach(fake.trypsin)
summary(fake.trypsin)
anova(lm(trypsin~grpf))
anova(lm(trypsin~grp))
model1 <- lm(trypsin~grp)
model2 <- lm(trypsin~grpf)
anova(model1,model2)
anova(lm(trypsin~grp+grpf))
xbar.trypsin <- tapply(trypsin,grpf,mean)
stripchart(trypsin~grp, method="jitter",
   jitter=.1, vertical=T, pch=20)
lines(1:6,xbar.trypsin,type="b",pch=4,cex=2,lty=2)
abline(lm(trypsin~grp))
if (.make.epsf) dev.copy2eps(file="trypsin.ps")
n <- c(32,137, 38,44,16,4)
tryp.mean <- c(128,152,194,207,215,218)
tryp.sd <-c(50.9,58.5,49.3,66.3,60,14)
gr<-1:6
anova(lm(tryp.mean~gr+factor(gr),weights=n))
sum(tryp.sd^2*(n-1))
sum(n-1)
sum(tryp.sd^2*(n-1))/sum(n-1)
206698/3318.007 # F statistic for gr
1-pf(206698/3318.007,1,265) # p-value
4351/3318.007   # F statistic for factor(gr)
1-pf(4351/3318.007,4,265) # p-value
attach(coking)
anova(lm(time~width*temp))
tapply(time,list(width,temp),mean)
hellung
summary(hellung)
hellung$glucose <- factor(hellung$glucose, labels=c("Yes","No"))
summary(hellung)
attach(hellung)
plot(conc,diameter,pch=as.numeric(glucose))
locator <- function(n)list(x=4e5,y=26)
legend(locator(n=1),legend=c("glucose","no glucose"),pch=1:2)
if (.make.epsf) dev.copy2eps(file="hellung-raw.ps")
plot(conc,diameter,pch=as.numeric(glucose),log="x")
if (.make.epsf) dev.copy2eps(file="hellung-logx.ps")
plot(conc,diameter,pch=as.numeric(glucose),log="xy")
tethym.gluc <- hellung[glucose=="Yes",]
tethym.nogluc <- hellung[glucose=="No",]
lm.nogluc <- lm(log10(diameter)~ log10(conc),data=tethym.nogluc)
lm.gluc <- lm(log10(diameter)~ log10(conc),data=tethym.gluc)
abline(lm.nogluc)
abline(lm.gluc)
if (.make.epsf) dev.copy2eps(file="hellung-loglog-lines.ps")
summary(lm(log10(diameter)~ log10(conc), data=tethym.gluc))
summary(lm(log10(diameter)~ log10(conc), data=tethym.nogluc))
summary(lm(log10(diameter)~log10(conc)*glucose))
summary(lm(log10(diameter)~log10(conc)+glucose))
var.test(lm.gluc,lm.nogluc)
anova(lm(log10(diameter)~ log10(conc)*glucose))
anova(lm(log10(diameter)~glucose+log10(conc)))
anova(lm(log10(diameter)~log10(conc)+ glucose))
t.test(log10(diameter)~glucose)
attach(thuesen)
options(na.action="na.exclude")
lm.velo <- lm(short.velocity~blood.glucose)
opar <- par(mfrow=c(2,2), mex=0.6, mar=c(4,4,3,2)+.3)
plot(lm.velo, which=1:4)
par(opar)
if (.make.epsf) dev.copy2eps(file="regr-diag.ps")
opar <- par(mfrow=c(2,2), mex=0.6, mar=c(4,4,3,2)+.3)
plot(rstandard(lm.velo))
plot(rstudent(lm.velo))
plot(dffits(lm.velo),type="l")
matplot(dfbetas(lm.velo),type="l", col="black")
lines(sqrt(cooks.distance(lm.velo)), lwd=2)
par(opar)
if (.make.epsf) dev.copy2eps(file="regr-diag2.ps")
summary(lm(short.velocity~blood.glucose, subset=-13))
cookd <- cooks.distance(lm(pemax~height+weight))
cookd <- cookd/max(cookd)
cook.colors <- gray(1-sqrt(cookd))
plot(height,weight,bg=cook.colors,pch=21,cex=1.5)
points(height,weight,pch=1,cex=1.5)
if (.make.epsf) dev.copy2eps(file="cookd-cyst-fibr.ps")
attach(secher)
rst <- rstudent(lm(log10(bwt)~log10(ad)+log10(bpd)))
range(rst)
rst <- rst/3.71
plot(ad,bpd,log="xy",bg=gray(1-abs(rst)),
     pch=ifelse(rst>0,24,25), cex=1.5)
if (.make.epsf) dev.copy2eps(file="rstudent-secher.ps")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
no.yes <- c("No","Yes")
smoking <- gl(2,1,8,no.yes)
obesity <- gl(2,2,8,no.yes)
snoring <- gl(2,4,8,no.yes)
n.tot <- c(60,17,8,2,187,85,51,23)
n.hyp <- c(5,2,1,0,35,13,15,8)
data.frame(smoking,obesity,snoring,n.tot,n.hyp)
expand.grid(smoking=no.yes, obesity=no.yes, snoring=no.yes)
hyp.tbl <- cbind(n.hyp,n.tot-n.hyp)
hyp.tbl
glm(hyp.tbl~smoking+obesity+snoring,family=binomial("logit"))
glm(hyp.tbl~smoking+obesity+snoring,binomial)
prop.hyp <- n.hyp/n.tot
glm.hyp <- glm(prop.hyp~smoking+obesity+snoring,
               binomial,weights=n.tot)
glm(hyp.tbl~smoking+obesity+snoring, binomial("logit"))
glm.hyp <- glm(hyp.tbl~smoking+obesity+snoring,binomial)
summary(glm.hyp)
summary(glm(formula = hyp.tbl ~ smoking + obesity + snoring, family =
binomial))
glm.hyp <- glm(hyp.tbl~obesity+snoring,binomial)
summary(glm.hyp)
glm.hyp <- glm(hyp.tbl~smoking+obesity+snoring,binomial)
anova(glm.hyp, test="Chisq")
glm.hyp <- glm(hyp.tbl~snoring+obesity+smoking,binomial)
anova(glm.hyp, test="Chisq")
glm.hyp <- glm(hyp.tbl~obesity+snoring,binomial)
anova(glm.hyp, test="Chisq")
drop1(glm.hyp, test="Chisq")
caesar.shoe
shoe.score <- 1:6
shoe.score
summary(glm(t(caesar.shoe)~shoe.score,binomial))
anova(glm(t(caesar.shoe)~shoe.score,binomial))
caesar.shoe.yes <- caesar.shoe["Yes",]
caesar.shoe.no <- caesar.shoe["No",]
caesar.shoe.total <- caesar.shoe.yes+caesar.shoe.no
prop.trend.test(caesar.shoe.yes,caesar.shoe.total)
prop.test(caesar.shoe.yes,caesar.shoe.total)
confint(glm.hyp)
confint.default(glm.hyp)
library(MASS)
plot(profile(glm.hyp))
if (.make.epsf) dev.copy2eps(file="profile-hyp.ps")
exp(cbind(OR=coef(glm.hyp), confint(glm.hyp)))
juul$menarche <- factor(juul$menarche, labels=c("No","Yes"))
juul$tanner <- factor(juul$tanner)
juul.girl <- subset(juul,age>8 & age<20 &
                    complete.cases(menarche))
attach(juul.girl)
summary(glm(menarche~age,binomial))
summary(glm(menarche~age+tanner,binomial))
drop1(glm(menarche~age+tanner,binomial),test="Chisq")
predict(glm.hyp)
predict(glm.hyp, type="response")
glm.menarche <- glm(menarche~age, binomial)
Age <- seq(8,20,.1)
newages <- data.frame(age=Age)
predicted.probability <- predict(glm.menarche,
                                 newages,type="resp")
plot(predicted.probability ~ Age, type="l")
if (.make.epsf) dev.copy2eps(file="menarche-fit.ps")
fitted(glm.hyp)
prop.hyp
fitted(glm.hyp)*n.tot
data.frame(fit=fitted(glm.hyp)*n.tot,n.hyp,n.tot)
age.group <- cut(age,c(8,10,12,13,14,15,16,18,20))
tb <- table(age.group,menarche)
tb
rel.freq <- prop.table(tb,1)[,2]
rel.freq
points(rel.freq ~ c(9,11,12.5,13.5,14.5,15.5,17,19),pch=5)
if (.make.epsf) dev.copy2eps(file="menarche-fit+obs.ps")
age.gr <- cut(age,c(8,12,13,14,20))
summary(glm(menarche~age+age.gr,binomial))
anova(glm(menarche~age+age.gr,binomial))
1-pchisq(8.058,3)
anova(glm(menarche~age+I(age^2)+I(age^3)+age.gr,binomial))
glm.menarche <- glm(menarche~age+I(age^2)+I(age^3), binomial)
predicted.probability <-
    predict(glm.menarche, newages, type="resp")
plot(predicted.probability ~ Age, type="l")
points(rel.freq~c(9,11,12.5,13.5,14.5,15.5,17,19), pch=5)
if (.make.epsf) dev.copy2eps(file="menarche-cubic.ps")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
library(survival)
attach(melanom)
names(melanom)
Surv(days, status==1)
survfit(Surv(days,status==1)~1)
surv.all <- survfit(Surv(days,status==1)~1)
summary(surv.all)
plot(surv.all)
if (.make.epsf) dev.copy2eps(file="surv-all.ps")
surv.bysex <- survfit(Surv(days,status==1)~sex)
plot(surv.bysex)
if (.make.epsf) dev.copy2eps(file="surv-bysex.ps")
plot(surv.bysex, conf.int=T, col=c("black","gray"))
survdiff(Surv(days,status==1)~sex)
survdiff(Surv(days,status==1)~sex+strata(ulc))
summary(coxph(Surv(days,status==1)~sex))
summary(coxph(Surv(days,status==1)~sex+log(thick)+strata(ulc)))
plot(survfit(coxph(Surv(days,status==1)~
             log(thick)+sex+strata(ulc))))
if (.make.epsf) dev.copy2eps(file="surv-cox.ps")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
names(eba1977)
attach(eba1977)
fit <- glm(cases~city+age+offset(log(pop)), family=poisson)
summary(fit)
min(fitted(fit))
pchisq(deviance(fit), df.residual(fit), lower=F)
pchisq(23.45, 15, lower=F)
drop1(fit, test="Chisq")
fit2 <- glm(cases~(city=="Fredericia")+age+offset(log(pop)),
                 family=poisson)
anova(fit, fit2, test="Chisq")
drop1(fit2,  test="Chisq")
summary(fit2)
cf <- coefficients(summary(fit2))
est <- cf[,1]
s.e. <- cf[,2]
rr <- exp(cbind(est, est - s.e.*qnorm(.975), est
                     + s.e.*qnorm(.975) ))
colnames(rr) <- c("RateRatio", "CI.lo","CI.hi")
rr
exp(cbind(coef(fit2), confint(fit2)))
head(nickel.expand)
subset(nickel.expand, id==325)
nickel.expand <- within(nickel.expand,
     lung.cancer <- as.numeric(icd %in% c(162,163)))
attach(nickel.expand)
pyr <- tapply(ageout-agein,list(ygr,agr), sum)
print(round(pyr), na.print="-")
count <- tapply(lung.cancer, list(ygr, agr), sum)
print(count, na.print="-")
print(round(count/pyr*1000, 1), na.print="-")
expect.count <- tapply(lung/1e6*(ageout-agein),
                       list(ygr,agr), sum)
print(round(expect.count, 1), na.print="-")
expect.tot <- sum(lung/1e6*(ageout-agein))
expect.tot
count.tot <- sum(lung.cancer)
count.tot
count.tot/expect.tot
fit <- glm(lung.cancer ~ 1, poisson,
           offset = log((ageout-agein)*lung/1e6))
summary(fit)
exp(coef(fit))
tapply(lung.cancer, agr, sum)
tapply(lung.cancer, ygr, sum)
detach()
nickel.expand <- within(nickel.expand,{
    A <- factor(agr)
    Y <- factor(ygr)
    lv <- levels(A)
    lv[1:6] <- "< 50"
    lv[11:13] <- "70+"
    levels(A) <- lv
    lv <- levels(Y)
    lv[7:10] <- "1961ff"
    levels(Y) <- lv
    rm(lv)
})
attach(nickel.expand)
fit <- glm(lung.cancer ~ A + Y, poisson,
           offset=log((ageout-agein)*lung/1e6))
drop1(fit, test="Chisq")
fit <- glm(lung.cancer ~ Y - 1, poisson,
           offset=log((ageout-agein)*lung/1e6))
summary(fit)
round(exp(coef(fit)), 1)
expect.count <-  tapply(lung/1e6*(ageout-agein), Y, sum)
count <- tapply(lung.cancer, Y, sum)
cbind(count=count, expect=round(expect.count,1),
      SMR= round(count/expect.count, 1))
detach()
nickel.expand <- within(nickel.expand,{
    TFE <- cut(agein-age1st, c(0,20,30,40,50,100), right=F)
    AFE <- cut(age1st, c(0, 20, 27.5, 35, 100), right=F)
    YFE <- cut(dob + age1st, c(0, 1910, 1915, 1920, 1925),right=F)
    EXP <- cut(exposure, c(0, 0.5, 4.5, 8.5, 12.5, 25), right=F)
})
attach(nickel.expand)
fit <- glm(lung.cancer ~ TFE + AFE + YFE + EXP, poisson,
           offset=log((ageout-agein)*lung/1e6))
drop1(fit, test="Chisq")
summary(fit)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
t <- 0:10
y <- rnorm(11, mean=5*exp(-t/5), sd=.2)
plot(y ~ t)
if (.make.epsf) dev.copy2eps(file="nonlin-sim.ps")
nlsout <- nls(y ~ A*exp(-alpha*t), start=c(A=2, alpha=0.05))
summary(nlsout)
attach(subset(juul2, age<20 & age>5 & sex==1))
plot(height ~ age)
if (.make.epsf) dev.copy2eps(file="juul-a-h.ps")
plot(log(5.3-log(height))~age)
if (.make.epsf) dev.copy2eps(file="gomp-dif.ps")
lm(log(5.3-log(height))~age)
fit <- nls(height~alpha*exp(-beta*exp(-gamma*age)),
           start=c(alpha=exp(5.3),beta=exp(0.42),gamma=0.15))
summary(fit)
plot(age, height)
newage <- seq(5,20,length=500)
lines(newage, predict(fit,newdata=data.frame(age=newage)),lwd=2)
if (.make.epsf) dev.copy2eps(file="gompertz.ps")

fit <- nls(log(height)~log(alpha*exp(-beta*exp(-gamma*age))),
start=c(alpha=exp(5.3),beta=exp(.12),gamma=.12))
summary(fit)
plot(age, log(height))
lines(newage, predict(fit,newdata=data.frame(age=newage)),lwd=2)
if (.make.epsf) dev.copy2eps(file="log-gompertz.ps")
# count quoted in text, subtract 1 for SSD()
length(ls(pattern="SS.*", "package:stats"))-1
summary(nls(height~SSgompertz(age, Asym, b2, b3)))
cf <- coef(nls(height ~ SSgompertz(age, Asym, b2, b3)))
summary(nls(log(height) ~
               log(as.vector(SSgompertz(age,Asym, b2, b3))),
            start=as.list(cf)))
par(mfrow=c(3,1))
plot(profile(fit))
if (.make.epsf) dev.copy2eps(file="gomp-prof.ps")
confint(fit)
confint.default(fit)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
