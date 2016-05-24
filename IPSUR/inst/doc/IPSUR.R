### R code from vignette source 'IPSUR.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: IPSUR.Rnw:217-230
###################################################
###  IPSUR.R - Introduction to Probability and Statistics Using R
###  Copyright (C) 2011  G. Jay Kerns, <gkerns@ysu.edu>
###  This program is free software: you can redistribute it and/or modify
###  it under the terms of the GNU General Public License as published by
###  the Free Software Foundation, either version 3 of the License, or
###  (at your option) any later version.
###  This program is distributed in the hope that it will be useful,
###  but WITHOUT ANY WARRANTY; without even the implied warranty of
###  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###  GNU General Public License for more details.
###  You should have received a copy of the GNU General Public License
###  along with this program.  If not, see <http://www.gnu.org/licenses/>
###################################################


###################################################
### code chunk number 2: IPSUR.Rnw:233-240
###################################################
seed <- 42
set.seed(seed)
options(width = 75)
#library(random)
#i_seed <- randomNumbers(n = 624, col = 1, min = -1e+09, max = 1e+09)
#.Random.seed[2:626] <- as.integer(c(1, i_seed))
#save.image(file = "seed.RData")


###################################################
### code chunk number 3: IPSUR.Rnw:243-274
###################################################
options(useFancyQuotes = FALSE)
#library(prob)
library(RcmdrPlugin.IPSUR)
# Generate RcmdrTestDrive
n <- 168
# generate order 
order <- 1:n
# generate race 
race <- sample(c("White","AfAmer","Asian","Other"), size=n, prob=c(76,13,5,6), replace = TRUE)
race <- factor(race)
# generate gender and smoke 
tmp <- sample(4, size=n, prob=c(12,38,9,41), replace = TRUE) 
gender <- factor(ifelse(tmp < 3,"Male", "Female")) 
smoke <- factor(ifelse(tmp %in% c(1,3), "Yes", "No"))
# generate parking
parking <- rgeom(n, prob = 0.4) + 1
# generate salary 
m <- 17 + (as.numeric(gender)-1) 
s <- 1 + (2 - as.numeric(gender)) 
salary <- rnorm(n, mean = m, sd = s)
# simulate reduction 
x <- arima.sim(list(order=c(1,0,0), ar=.9), n=n) 
reduction <- as.numeric((20*x + order)/n + 5)
# simulate before and after
before <- rlogis(n, location = 68, scale = 3) 
m <- (as.numeric(smoke)-1)*2.5 
after <- before - rnorm(n, mean = m, sd=0.1)
RcmdrTestDrive <- data.frame(order = order, race = race, smoke = smoke, gender = gender, salary = salary, reduction = reduction, before = before, after = after, parking = parking)
# clean up
remove(list = names(RcmdrTestDrive))
remove(x, n, m, s, tmp)


###################################################
### code chunk number 4: IPSUR.Rnw:277-377
###################################################
plot.htest <- function (x, hypoth.or.conf = 'Hypoth',...) { 
require("HH") 
if (x$method == "1-sample proportions test with continuity correction" || x$method == "1-sample proportions test without continuity correction"){
mu <- x$null.value
obs.mean <- x$estimate
n <- NA
std.dev <- abs(obs.mean - mu)/sqrt(x$statistic)
deg.freedom <- NA
if(x$alternative == "two.sided"){
alpha.right <- (1 - attr(x$conf.int, "conf.level"))/2
Use.alpha.left <- TRUE
Use.alpha.right <- TRUE
} else if (x$alternative == "less") {
alpha.right <- 1 - attr(x$conf.int, "conf.level")
Use.alpha.left <- TRUE
Use.alpha.right <- FALSE
} else {
alpha.right <- 1 - attr(x$conf.int, "conf.level")
Use.alpha.left <- FALSE
Use.alpha.right <- TRUE
}
} else if (x$method == "One Sample z-test"){
mu <- x$null.value
obs.mean <- x$estimate
n <- x$parameter[1]
std.dev <- x$parameter[2]
deg.freedom <- NA
if(x$alternative == "two.sided"){
alpha.right <- (1 - attr(x$conf.int, "conf.level"))/2
Use.alpha.left <- TRUE
Use.alpha.right <- TRUE
} else if (x$alternative == "less") {
alpha.right <- 1 - attr(x$conf.int, "conf.level")
Use.alpha.left <- TRUE
Use.alpha.right <- FALSE
} else {
alpha.right <- 1 - attr(x$conf.int, "conf.level")
Use.alpha.left <- FALSE
Use.alpha.right <- TRUE
} 
} else if (x$method == "One Sample t-test" || x$method == "Paired t-test"){
mu <- x$null.value
obs.mean <- x$estimate
n <- x$parameter + 1
std.dev <- x$estimate/x$statistic*sqrt(n)
deg.freedom <- x$parameter
if(x$alternative == "two.sided"){
alpha.right <- (1 - attr(x$conf.int, "conf.level"))/2
Use.alpha.left <- TRUE
Use.alpha.right <- TRUE
} else if (x$alternative == "less") {
alpha.right <- 1 - attr(x$conf.int, "conf.level")
Use.alpha.left <- TRUE
Use.alpha.right <- FALSE
} else {
alpha.right <- 1 - attr(x$conf.int, "conf.level")
Use.alpha.left <- FALSE
Use.alpha.right <- TRUE
}
} else if (x$method == "Welch Two Sample t-test"){
mu <- x$null.value
obs.mean <- -diff(x$estimate)
n <- x$parameter + 2
std.dev <- obs.mean/x$statistic*sqrt(n)
deg.freedom <- x$parameter
if(x$alternative == "two.sided"){
alpha.right <- (1 - attr(x$conf.int, "conf.level"))/2
Use.alpha.left <- TRUE
Use.alpha.right <- TRUE
} else if (x$alternative == "less") {
alpha.right <- 1 - attr(x$conf.int, "conf.level")
Use.alpha.left <- TRUE
Use.alpha.right <- FALSE
} else {
alpha.right <- 1 - attr(x$conf.int, "conf.level")
Use.alpha.left <- FALSE
Use.alpha.right <- TRUE
} 
} else if (x$method == " Two Sample t-test"){
mu <- x$null.value
obs.mean <- -diff(x$estimate)
n <- x$parameter + 2
std.dev <- obs.mean/x$statistic*sqrt(n)
deg.freedom <- x$parameter
if(x$alternative == "two.sided"){
alpha.right <- (1 - attr(x$conf.int, "conf.level"))/2
Use.alpha.left <- TRUE
Use.alpha.right <- TRUE
} else if (x$alternative == "less") {
alpha.right <- 1 - attr(x$conf.int, "conf.level")
Use.alpha.left <- TRUE
Use.alpha.right <- FALSE
} else {
alpha.right <- 1 - attr(x$conf.int, "conf.level")
Use.alpha.left <- FALSE
Use.alpha.right <- TRUE
}
}
return(normal.and.t.dist(mu.H0 = mu, obs.mean = obs.mean, std.dev = std.dev, n = n, deg.freedom = deg.freedom, alpha.right = alpha.right, Use.obs.mean = TRUE, Use.alpha.left = Use.alpha.left, Use.alpha.right = Use.alpha.right, hypoth.or.conf = hypoth.or.conf))
}


###################################################
### code chunk number 5: IPSUR.Rnw:645-648 (eval = FALSE)
###################################################
## install.packages("IPSUR")
## library(IPSUR)
## read(IPSUR)


###################################################
### code chunk number 6: IPSUR.Rnw:871-872
###################################################
getOption("defaultPackages")


###################################################
### code chunk number 7: IPSUR.Rnw:1061-1064
###################################################
2 + 3       # add
4 * 5 / 6   # multiply and divide
7^8         # 7 to the 8th power


###################################################
### code chunk number 8: IPSUR.Rnw:1073-1079
###################################################
options(digits = 16)
10/3                 # see more digits
sqrt(2)              # square root
exp(1)               # Euler's constant, e
pi       
options(digits = 7)  # back to default


###################################################
### code chunk number 9: IPSUR.Rnw:1102-1104
###################################################
x <- 7*41/pi   # don't see the calculated value
x              # take a look


###################################################
### code chunk number 10: five
###################################################
sqrt(-1)              # isn't defined
sqrt(-1+0i)           # is defined
sqrt(as.complex(-1))  # same thing
(0 + 1i)^2            # should be -1
typeof((0 + 1i)^2)


###################################################
### code chunk number 11: IPSUR.Rnw:1167-1169
###################################################
x <- c(74, 31, 95, 61, 76, 34, 23, 54, 96)
x


###################################################
### code chunk number 12: IPSUR.Rnw:1198-1200
###################################################
seq(from = 1, to = 5)
seq(from = 2, by = -0.1, length.out = 4)


###################################################
### code chunk number 13: IPSUR.Rnw:1206-1207
###################################################
1:5


###################################################
### code chunk number 14: IPSUR.Rnw:1223-1227
###################################################
x[1]
x[2:4]
x[c(1,3,4,8)]
x[-c(1,3,4,8)]


###################################################
### code chunk number 15: IPSUR.Rnw:1233-1236
###################################################
LETTERS[1:5]
letters[-(6:24)]



###################################################
### code chunk number 16: IPSUR.Rnw:1246-1252
###################################################
x <- 1:5
sum(x)
length(x)
min(x)
mean(x)      # sample mean
sd(x)        # sample standard deviation


###################################################
### code chunk number 17: IPSUR.Rnw:1268-1269
###################################################
intersect


###################################################
### code chunk number 18: IPSUR.Rnw:1279-1280
###################################################
rev


###################################################
### code chunk number 19: IPSUR.Rnw:1287-1288
###################################################
methods(rev)


###################################################
### code chunk number 20: IPSUR.Rnw:1300-1301
###################################################
rev.default


###################################################
### code chunk number 21: IPSUR.Rnw:1310-1312
###################################################
wilcox.test
methods(wilcox.test)


###################################################
### code chunk number 22: IPSUR.Rnw:1331-1332
###################################################
exp


###################################################
### code chunk number 23: IPSUR.Rnw:1616-1618
###################################################
str(precip)
precip[1:4]


###################################################
### code chunk number 24: IPSUR.Rnw:1637-1638
###################################################
str(rivers)


###################################################
### code chunk number 25: IPSUR.Rnw:1656-1658
###################################################
str(discoveries)
discoveries[1:4]


###################################################
### code chunk number 26: IPSUR.Rnw:1697-1700 (eval = FALSE)
###################################################
## stripchart(precip, xlab="rainfall")
## stripchart(rivers, method="jitter", xlab="length")
## stripchart(discoveries, method="stack", xlab="number")


###################################################
### code chunk number 27: IPSUR.Rnw:1718-1723
###################################################
par(mfrow = c(1,3)) # 3 plots: 1 row, 3 columns
stripchart(precip, xlab="rainfall")
stripchart(rivers, method="jitter", xlab="length")
stripchart(discoveries, method="stack", xlab="number")
par(mfrow = c(1,1)) # back to normal


###################################################
### code chunk number 28: IPSUR.Rnw:1767-1769 (eval = FALSE)
###################################################
## hist(precip, main = "")
## hist(precip, freq = FALSE, main = "")


###################################################
### code chunk number 29: IPSUR.Rnw:1781-1785
###################################################
par(mfrow = c(1,2)) # 2 plots: 1 row, 2 columns
hist(precip, main = "")
hist(precip, freq = FALSE, main = "")
par(mfrow = c(1,1)) # back to normal


###################################################
### code chunk number 30: IPSUR.Rnw:1814-1816 (eval = FALSE)
###################################################
## hist(precip, breaks = 10, main = "")
## hist(precip, breaks = 200, main = "")


###################################################
### code chunk number 31: IPSUR.Rnw:1822-1826
###################################################
par(mfrow = c(1,2)) # 2 plots: 1 row, 2 columns
hist(precip, breaks = 10, main = "")
hist(precip, breaks = 200, main = "")
par(mfrow = c(1,1)) # back to normal


###################################################
### code chunk number 32: IPSUR.Rnw:1869-1871
###################################################
library(aplpack)
stem.leaf(UKDriverDeaths, depth = FALSE)


###################################################
### code chunk number 33: IPSUR.Rnw:1908-1910
###################################################
plot(LakeHuron, type = "h")
plot(LakeHuron, type = "p")


###################################################
### code chunk number 34: IPSUR.Rnw:1919-1923
###################################################
par(mfrow = c(2,1)) # 2 plots: 1 row, 2 columns
plot(LakeHuron, type = "h")
plot(LakeHuron, type = "p")
par(mfrow = c(1,1)) # back to normal


###################################################
### code chunk number 35: IPSUR.Rnw:1985-1986
###################################################
str(state.abb)


###################################################
### code chunk number 36: IPSUR.Rnw:2002-2004
###################################################
str(state.region)
state.region[1:5]


###################################################
### code chunk number 37: IPSUR.Rnw:2032-2036
###################################################
Tbl <- table(state.division)
Tbl               # frequencies
Tbl/sum(Tbl)      # relative frequencies
prop.table(Tbl)   # same thing


###################################################
### code chunk number 38: IPSUR.Rnw:2056-2058 (eval = FALSE)
###################################################
## barplot(table(state.region), cex.names = 0.50)
## barplot(prop.table(table(state.region)), cex.names = 0.50)


###################################################
### code chunk number 39: IPSUR.Rnw:2077-2081
###################################################
par(mfrow = c(1,2)) # 2 plots: 1 row, 2 columns
barplot(table(state.region), cex.names = 0.50)
barplot(prop.table(table(state.region)), cex.names = 0.50)
par(mfrow = c(1,1))


###################################################
### code chunk number 40: IPSUR.Rnw:2115-2117 (eval = FALSE)
###################################################
## library(qcc)
## pareto.chart(table(state.division), ylab="Frequency")


###################################################
### code chunk number 41: IPSUR.Rnw:2123-2125
###################################################
library(qcc)
pareto.chart(table(state.division), ylab="Frequency")


###################################################
### code chunk number 42: IPSUR.Rnw:2147-2149 (eval = FALSE)
###################################################
## x <- table(state.region)
## dotchart(as.vector(x), labels = names(x))


###################################################
### code chunk number 43: IPSUR.Rnw:2155-2157
###################################################
x <- table(state.region)
dotchart(as.vector(x), labels = names(x))


###################################################
### code chunk number 44: IPSUR.Rnw:2193-2196
###################################################
x <- 5:9
y <- (x < 7.3)
y


###################################################
### code chunk number 45: IPSUR.Rnw:2211-2212
###################################################
!y


###################################################
### code chunk number 46: IPSUR.Rnw:2229-2232
###################################################
x <- c(3, 7, NA, 4, 7)
y <- c(5, NA, 1, 2, 2)
x + y


###################################################
### code chunk number 47: IPSUR.Rnw:2244-2246
###################################################
sum(x)
sum(x, na.rm = TRUE)


###################################################
### code chunk number 48: IPSUR.Rnw:2257-2260
###################################################
is.na(x)
z <- x[!is.na(x)]
sum(z)


###################################################
### code chunk number 49: IPSUR.Rnw:2357-2359
###################################################
library(aplpack)
stem.leaf(faithful$eruptions)


###################################################
### code chunk number 50: IPSUR.Rnw:2725-2728
###################################################
library(e1071)
skewness(discoveries)
2*sqrt(6/length(discoveries))


###################################################
### code chunk number 51: IPSUR.Rnw:2735-2737
###################################################
kurtosis(UKDriverDeaths)
4*sqrt(6/length(UKDriverDeaths))


###################################################
### code chunk number 52: IPSUR.Rnw:2819-2820
###################################################
stem.leaf(rivers)


###################################################
### code chunk number 53: IPSUR.Rnw:2844-2845
###################################################
stem.leaf(precip)


###################################################
### code chunk number 54: IPSUR.Rnw:2957-2958
###################################################
boxplot.stats(rivers)$out


###################################################
### code chunk number 55: IPSUR.Rnw:2964-2965
###################################################
boxplot.stats(rivers, coef = 3)$out


###################################################
### code chunk number 56: IPSUR.Rnw:3016-3019
###################################################
x <- 5:8
y <- letters[3:6]
A <- data.frame(v1 = x, v2 = y)


###################################################
### code chunk number 57: IPSUR.Rnw:3039-3042
###################################################
A[3,]
A[1, ]
A[ ,2]


###################################################
### code chunk number 58: IPSUR.Rnw:3057-3059
###################################################
names(A)
A$v1


###################################################
### code chunk number 59: IPSUR.Rnw:3210-3212 (eval = FALSE)
###################################################
## library(lattice)
## xyplot()


###################################################
### code chunk number 60: IPSUR.Rnw:3292-3294 (eval = FALSE)
###################################################
## library(lattice)
## bwplot(~weight | feed, data = chickwts)


###################################################
### code chunk number 61: IPSUR.Rnw:3300-3302
###################################################
library(lattice)
print(bwplot(~ weight | feed, data = chickwts))


###################################################
### code chunk number 62: IPSUR.Rnw:3317-3318 (eval = FALSE)
###################################################
## histogram(~age | education, data = infert)


###################################################
### code chunk number 63: IPSUR.Rnw:3324-3326
###################################################
library(lattice)
print(histogram(~age | education, data = infert))


###################################################
### code chunk number 64: IPSUR.Rnw:3339-3340 (eval = FALSE)
###################################################
## xyplot(Petal.Length ~ Petal.Width | Species, data = iris)


###################################################
### code chunk number 65: IPSUR.Rnw:3346-3348
###################################################
library(lattice)
print(xyplot(Petal.Length ~ Petal.Width | Species, data = iris))


###################################################
### code chunk number 66: IPSUR.Rnw:3361-3362 (eval = FALSE)
###################################################
## coplot(conc ~ uptake | Type * Treatment, data = CO2)


###################################################
### code chunk number 67: IPSUR.Rnw:3368-3370
###################################################
library(lattice)
print(coplot(conc ~ uptake | Type * Treatment, data = CO2))


###################################################
### code chunk number 68: IPSUR.Rnw:3401-3403
###################################################
attach(RcmdrTestDrive)
names(RcmdrTestDrive)


###################################################
### code chunk number 69: "Find summary statistics"
###################################################
summary(RcmdrTestDrive)


###################################################
### code chunk number 70: IPSUR.Rnw:3450-3451
###################################################
table(race)


###################################################
### code chunk number 71: IPSUR.Rnw:3461-3462
###################################################
barplot(table(RcmdrTestDrive$race), main="", xlab="race", ylab="Frequency", legend.text=FALSE, col=NULL) 


###################################################
### code chunk number 72: IPSUR.Rnw:3498-3500
###################################################
x <- tapply(salary, list(gender = gender), mean)
x


###################################################
### code chunk number 73: IPSUR.Rnw:3506-3507
###################################################
by(salary, gender, mean, na.rm = TRUE)


###################################################
### code chunk number 74: IPSUR.Rnw:3527-3528
###################################################
x[which(x==max(x))]


###################################################
### code chunk number 75: IPSUR.Rnw:3536-3538
###################################################
y <- tapply(salary, list(gender = gender), sd)
y


###################################################
### code chunk number 76: IPSUR.Rnw:3552-3553
###################################################
boxplot(salary~gender, xlab="salary", ylab="gender", main="", notch=FALSE, varwidth=TRUE, horizontal=TRUE, data=RcmdrTestDrive) 


###################################################
### code chunk number 77: IPSUR.Rnw:3589-3590
###################################################
x = sort(reduction)


###################################################
### code chunk number 78: IPSUR.Rnw:3593-3597
###################################################
x[137]
IQR(x)
fivenum(x)
fivenum(x)[4] - fivenum(x)[2]


###################################################
### code chunk number 79: IPSUR.Rnw:3609-3610
###################################################
boxplot(reduction, xlab="reduction", main="", notch=FALSE, varwidth=TRUE, horizontal=TRUE, data=RcmdrTestDrive) 


###################################################
### code chunk number 80: IPSUR.Rnw:3614-3619
###################################################
temp <- fivenum(x)
inF <- 1.5 * (temp[4] - temp[2]) + temp[4]
outF <- 3 * (temp[4] - temp[2]) + temp[4]
which(x > inF)
which(x > outF)


###################################################
### code chunk number 81: IPSUR.Rnw:3660-3662
###################################################
c(mean(before), median(before))
c(mean(after), median(after))


###################################################
### code chunk number 82: IPSUR.Rnw:3680-3681
###################################################
boxplot(before, xlab="before", main="", notch=FALSE, varwidth=TRUE, horizontal=TRUE, data=RcmdrTestDrive) 


###################################################
### code chunk number 83: IPSUR.Rnw:3696-3697
###################################################
boxplot(after, xlab="after", notch=FALSE, varwidth=TRUE, horizontal=TRUE, data=RcmdrTestDrive) 


###################################################
### code chunk number 84: IPSUR.Rnw:3718-3721
###################################################
sd(before)
mad(after)
IQR(after)/1.349


###################################################
### code chunk number 85: IPSUR.Rnw:3736-3739
###################################################
library(e1071)
skewness(before)
kurtosis(before)


###################################################
### code chunk number 86: IPSUR.Rnw:3758-3760
###################################################
skewness(after)
kurtosis(after)


###################################################
### code chunk number 87: IPSUR.Rnw:3775-3776
###################################################
hist(before, xlab="before", data=RcmdrTestDrive) 


###################################################
### code chunk number 88: IPSUR.Rnw:3781-3782
###################################################
hist(after, xlab="after", data=RcmdrTestDrive) 


###################################################
### code chunk number 89: IPSUR.Rnw:3836-3844
###################################################
require(diagram)
par(mex = 0.2, cex = 0.5)
openplotmat(frame.plot=TRUE)
straightarrow(from = c(0.46,0.74), to = c(0.53,0.71), arr.pos = 1)
straightarrow(from = c(0.3,0.65), to = c(0.3,0.51), arr.pos = 1)
textellipse(mid = c(0.74,0.55), box.col = grey(0.95), radx = 0.24, rady = 0.22, lab = c(expression(bold(underline(DETERMINISTIC))), expression(2*H[2]+O[2] %->% H[2]*O), "3 + 4 = 7"), cex = 2 )
textrect(mid = c(0.3, 0.75), radx = 0.15, rady = 0.1, lab = c("Experiments"), cex = 2 )
textellipse(mid = c(0.29,0.25), box.col = grey(0.95), radx = 0.27, rady = 0.22, lab = c(expression(bold(underline(RANDOM))), "toss coin, roll die", "count ants on sidewalk", "measure rainfall" ), cex = 2 )


###################################################
### code chunk number 90: IPSUR.Rnw:3884-3886
###################################################
S <- data.frame(lands = c("down","up","side"))
S


###################################################
### code chunk number 91: IPSUR.Rnw:3916-3918
###################################################
library(prob)
tosscoin(1) 


###################################################
### code chunk number 92: IPSUR.Rnw:3925-3926
###################################################
tosscoin(3) 


###################################################
### code chunk number 93: IPSUR.Rnw:3931-3932
###################################################
rolldie(1) 


###################################################
### code chunk number 94: IPSUR.Rnw:3944-3945
###################################################
head(cards()) 


###################################################
### code chunk number 95: IPSUR.Rnw:4016-4017
###################################################
urnsamples(1:3, size = 2, replace = TRUE, ordered = TRUE)


###################################################
### code chunk number 96: IPSUR.Rnw:4034-4035
###################################################
urnsamples(1:3, size = 2, replace = FALSE, ordered = TRUE)


###################################################
### code chunk number 97: IPSUR.Rnw:4053-4054
###################################################
urnsamples(1:3, size = 2, replace = FALSE, ordered = FALSE) 


###################################################
### code chunk number 98: IPSUR.Rnw:4068-4069
###################################################
urnsamples(1:3, size = 2, replace = TRUE, ordered = FALSE) 


###################################################
### code chunk number 99: IPSUR.Rnw:4131-4134
###################################################
S <- tosscoin(2, makespace = TRUE) 
S[1:3, ] 
S[c(2,4), ] 


###################################################
### code chunk number 100: IPSUR.Rnw:4141-4142
###################################################
S <- cards() 


###################################################
### code chunk number 101: IPSUR.Rnw:4145-4147
###################################################
subset(S, suit == "Heart") 
subset(S, rank %in% 7:9)


###################################################
### code chunk number 102: IPSUR.Rnw:4153-4154
###################################################
subset(rolldie(3), X1+X2+X3 > 16) 


###################################################
### code chunk number 103: IPSUR.Rnw:4174-4177
###################################################
x <- 1:10 
y <- 8:12 
y %in% x


###################################################
### code chunk number 104: IPSUR.Rnw:4194-4195
###################################################
isin(x,y) 


###################################################
### code chunk number 105: IPSUR.Rnw:4203-4205
###################################################
x <- 1:10 
y <- c(3,3,7) 


###################################################
### code chunk number 106: IPSUR.Rnw:4208-4210
###################################################
all(y %in% x)
isin(x,y) 


###################################################
### code chunk number 107: IPSUR.Rnw:4223-4225
###################################################
isin(x, c(3,4,5), ordered = TRUE) 
isin(x, c(3,5,4), ordered = TRUE) 


###################################################
### code chunk number 108: IPSUR.Rnw:4234-4236
###################################################
S <- rolldie(4) 
subset(S, isin(S, c(2,2,6), ordered = TRUE)) 


###################################################
### code chunk number 109: IPSUR.Rnw:4268-4271
###################################################
S = cards() 
A = subset(S, suit == "Heart") 
B = subset(S, rank %in% 7:9)


###################################################
### code chunk number 110: IPSUR.Rnw:4276-4280
###################################################
union(A,B) 
intersect(A,B) 
setdiff(A,B) 
setdiff(B,A) 


###################################################
### code chunk number 111: IPSUR.Rnw:4546-4549
###################################################
outcomes <- rolldie(1) 
p <- rep(1/6, times = 6) 
probspace(outcomes, probs = p) 


###################################################
### code chunk number 112: IPSUR.Rnw:4558-4559
###################################################
probspace(1:6, probs = p) 


###################################################
### code chunk number 113: IPSUR.Rnw:4567-4568
###################################################
probspace(1:6) 


###################################################
### code chunk number 114: IPSUR.Rnw:4577-4578
###################################################
rolldie(1, makespace = TRUE)


###################################################
### code chunk number 115: IPSUR.Rnw:4607-4608
###################################################
probspace(tosscoin(1), probs = c(0.70, 0.30)) 


###################################################
### code chunk number 116: IPSUR.Rnw:4853-4856
###################################################
S <- cards(makespace = TRUE) 
A <- subset(S, suit == "Heart") 
B <- subset(S, rank %in% 7:9)


###################################################
### code chunk number 117: IPSUR.Rnw:4861-4862
###################################################
prob(A) 


###################################################
### code chunk number 118: IPSUR.Rnw:4867-4868
###################################################
prob(S, suit == "Heart") 


###################################################
### code chunk number 119: IPSUR.Rnw:5067-5071
###################################################
nsamp(n=3, k=2, replace = TRUE, ordered = TRUE) 
nsamp(n=3, k=2, replace = FALSE, ordered = TRUE) 
nsamp(n=3, k=2, replace = FALSE, ordered = FALSE) 
nsamp(n=3, k=2, replace = TRUE, ordered = FALSE) 


###################################################
### code chunk number 120: IPSUR.Rnw:5106-5109
###################################################
n <- c(11,7,31) 
k <- c(3,4,3) 
r <- c(FALSE,FALSE,TRUE) 


###################################################
### code chunk number 121: IPSUR.Rnw:5112-5113
###################################################
x <- nsamp(n, k, rep = r, ord = TRUE) 


###################################################
### code chunk number 122: IPSUR.Rnw:5125-5126
###################################################
prod(x) 


###################################################
### code chunk number 123: IPSUR.Rnw:5131-5132
###################################################
(11*10*9)*(7*6*5*4)*313 


###################################################
### code chunk number 124: IPSUR.Rnw:5137-5138
###################################################
prod(9:11)*prod(4:7)*313 


###################################################
### code chunk number 125: IPSUR.Rnw:5143-5144
###################################################
prod(factorial(c(11,7))/factorial(c(8,3)))*313 


###################################################
### code chunk number 126: IPSUR.Rnw:5219-5224
###################################################
g <- Vectorize(pbirthday.ipsur)
plot(1:50, g(1:50), xlab = "Number of people in room", ylab = "Prob(at least one match)")
abline(h = 0.5)
abline(v = 23, lty = 2)
remove(g)


###################################################
### code chunk number 127: IPSUR.Rnw:5357-5360
###################################################
library(prob)
S <- rolldie(2, makespace = TRUE)  # assumes ELM
head(S)                            #  first few rows


###################################################
### code chunk number 128: IPSUR.Rnw:5365-5367
###################################################
A <- subset(S, X1 == X2)
B <- subset(S, X1 + X2 >= 8)


###################################################
### code chunk number 129: IPSUR.Rnw:5375-5377
###################################################
prob(A, given = B)
prob(B, given = A)


###################################################
### code chunk number 130: IPSUR.Rnw:5385-5387
###################################################
prob(S, X1==X2, given = (X1 + X2 >= 8) )
prob(S, X1+X2 >= 8, given = (X1==X2) )


###################################################
### code chunk number 131: IPSUR.Rnw:5471-5475
###################################################
library(prob)
L <- cards()
M <- urnsamples(L, size = 2)
N <- probspace(M)


###################################################
### code chunk number 132: IPSUR.Rnw:5487-5488
###################################################
prob(N, all(rank == "A"))


###################################################
### code chunk number 133: IPSUR.Rnw:5523-5527
###################################################
library(prob)
L <- rep(c("red","green"), times = c(7,3))
M <- urnsamples(L, size = 3, replace = FALSE, ordered = TRUE)
N <- probspace(M)


###################################################
### code chunk number 134: IPSUR.Rnw:5543-5544
###################################################
prob(N, isrep(N, "red", 3))


###################################################
### code chunk number 135: IPSUR.Rnw:5551-5552
###################################################
prob(N, isrep(N, "red", 2))


###################################################
### code chunk number 136: IPSUR.Rnw:5561-5562
###################################################
prob(N, isin(N, c("red","green","red"), ordered = TRUE))


###################################################
### code chunk number 137: IPSUR.Rnw:5571-5572
###################################################
prob(N, isin(N, c("red","green","red")))


###################################################
### code chunk number 138: IPSUR.Rnw:5603-5606
###################################################
.Table <- xtabs(~smoke+gender, data=RcmdrTestDrive)
addmargins(.Table) # Table with Marginal Distributions
remove(.Table)


###################################################
### code chunk number 139: IPSUR.Rnw:5724-5727
###################################################
S <- tosscoin(10, makespace = TRUE)
A <- subset(S, isrep(S, vals = "T", nrep = 10))
1 - prob(A)


###################################################
### code chunk number 140: IPSUR.Rnw:5765-5766
###################################################
iidspace(c("H","T"), ntrials = 3, probs = c(0.7, 0.3)) 


###################################################
### code chunk number 141: IPSUR.Rnw:5973-5977
###################################################
prior <- c(0.6, 0.3, 0.1)
like <- c(0.003, 0.007, 0.010)
post <- prior * like
post / sum(post)


###################################################
### code chunk number 142: IPSUR.Rnw:5997-6000
###################################################
newprior <- post
post <- newprior * like^7
post / sum(post)


###################################################
### code chunk number 143: IPSUR.Rnw:6021-6023
###################################################
fastpost <- prior * like^8
fastpost / sum(fastpost)


###################################################
### code chunk number 144: IPSUR.Rnw:6117-6119
###################################################
S <- rolldie(3, nsides = 4, makespace = TRUE) 
S <- addrv(S, U = X1-X2+X3) 


###################################################
### code chunk number 145: IPSUR.Rnw:6126-6127
###################################################
head(S)


###################################################
### code chunk number 146: IPSUR.Rnw:6133-6134
###################################################
prob(S, U > 6) 


###################################################
### code chunk number 147: IPSUR.Rnw:6152-6155
###################################################
S <- addrv(S, FUN = max, invars = c("X1","X2","X3"), name = "V") 
S <- addrv(S, FUN = sum, invars = c("X1","X2","X3"), name = "W") 
head(S) 


###################################################
### code chunk number 148: IPSUR.Rnw:6183-6184
###################################################
marginal(S, vars = "V") 


###################################################
### code chunk number 149: IPSUR.Rnw:6193-6194
###################################################
marginal(S, vars = c("V", "W")) 


###################################################
### code chunk number 150: IPSUR.Rnw:6210-6211
###################################################
rnorm(1)


###################################################
### code chunk number 151: IPSUR.Rnw:6375-6377
###################################################
x <- c(0,1,2,3)
f <- c(1/8, 3/8, 3/8, 1/8)


###################################################
### code chunk number 152: IPSUR.Rnw:6386-6388
###################################################
mu <- sum(x * f)
mu


###################################################
### code chunk number 153: IPSUR.Rnw:6397-6401
###################################################
sigma2 <- sum((x-mu)^2 * f)
sigma2
sigma <- sqrt(sigma2)
sigma


###################################################
### code chunk number 154: IPSUR.Rnw:6408-6410
###################################################
F = cumsum(f)
F


###################################################
### code chunk number 155: IPSUR.Rnw:6420-6423
###################################################
library(distrEx)
X <- DiscreteDistribution(supp = 0:3, prob = c(1,3,3,1)/8)
E(X); var(X); sd(X)


###################################################
### code chunk number 156: IPSUR.Rnw:6586-6589
###################################################
A <- data.frame(Pr=dbinom(0:4, size = 4, prob = 0.5))
rownames(A) <- 0:4 
A


###################################################
### code chunk number 157: IPSUR.Rnw:6613-6615
###################################################
pbinom(9, size=12, prob=1/6) - pbinom(6, size=12, prob=1/6)
diff(pbinom(c(6,9), size = 12, prob = 1/6))  # same thing


###################################################
### code chunk number 158: IPSUR.Rnw:6662-6667
###################################################
plot(0, xlim = c(-1.2, 4.2), ylim = c(-0.04, 1.04), type = "n", xlab = "number of successes", ylab = "cumulative probability")
abline(h = c(0,1), lty = 2, col = "grey")
lines(stepfun(0:3, pbinom(-1:3, size = 3, prob = 0.5)), verticals = FALSE, do.p = FALSE)
points(0:3, pbinom(0:3, size = 3, prob = 0.5), pch = 16, cex = 1.2)
points(0:3, pbinom(-1:2, size = 3, prob = 0.5), pch = 1, cex = 1.2)


###################################################
### code chunk number 159: IPSUR.Rnw:6689-6692
###################################################
library(distr)
X <- Binom(size = 3, prob = 1/2)
X


###################################################
### code chunk number 160: IPSUR.Rnw:6702-6704
###################################################
d(X)(1)   # pmf of X evaluated at x = 1
p(X)(2)   # cdf of X evaluated at x = 2


###################################################
### code chunk number 161: IPSUR.Rnw:6716-6717
###################################################
plot(X, cex = 0.2)


###################################################
### code chunk number 162: IPSUR.Rnw:6926-6930
###################################################
X <- Binom(size = 3, prob = 0.45)
library(distrEx)
E(X)
E(3*X + 4)


###################################################
### code chunk number 163: IPSUR.Rnw:6942-6944
###################################################
var(X)
sd(X)


###################################################
### code chunk number 164: IPSUR.Rnw:6993-6995
###################################################
x <- c(4, 7, 9, 11, 12)
ecdf(x)


###################################################
### code chunk number 165: IPSUR.Rnw:7008-7009 (eval = FALSE)
###################################################
## plot(ecdf(x))


###################################################
### code chunk number 166: IPSUR.Rnw:7015-7016
###################################################
plot(ecdf(x))


###################################################
### code chunk number 167: IPSUR.Rnw:7034-7037
###################################################
epdf <- function(x) function(t){sum(x %in% t)/length(x)}
x <- c(0,0,1)
epdf(x)(0)       # should be 2/3


###################################################
### code chunk number 168: IPSUR.Rnw:7045-7047
###################################################
x <- c(0,0,1)
sample(x, size = 7, replace = TRUE)


###################################################
### code chunk number 169: IPSUR.Rnw:7122-7123
###################################################
dhyper(3, m = 17, n = 233, k = 5)


###################################################
### code chunk number 170: IPSUR.Rnw:7133-7136
###################################################
A <- data.frame(Pr=dhyper(0:4, m = 17, n = 233, k = 5))
rownames(A) <- 0:4 
A


###################################################
### code chunk number 171: IPSUR.Rnw:7150-7151
###################################################
dhyper(5, m = 17, n = 233, k = 5)


###################################################
### code chunk number 172: IPSUR.Rnw:7171-7172
###################################################
phyper(2, m = 17, n = 233, k = 5)


###################################################
### code chunk number 173: IPSUR.Rnw:7189-7190
###################################################
phyper(1, m = 17, n = 233, k = 5, lower.tail = FALSE)


###################################################
### code chunk number 174: IPSUR.Rnw:7244-7245
###################################################
rhyper(10, m = 17, n = 233, k = 5)


###################################################
### code chunk number 175: IPSUR.Rnw:7320-7321
###################################################
pgeom(4, prob = 0.812, lower.tail = FALSE)


###################################################
### code chunk number 176: IPSUR.Rnw:7372-7373
###################################################
dnbinom(5, size = 7, prob = 0.5)


###################################################
### code chunk number 177: IPSUR.Rnw:7483-7484
###################################################
diff(ppois(c(47, 50), lambda = 50))


###################################################
### code chunk number 178: IPSUR.Rnw:7598-7605
###################################################
xmin <- qbinom(.0005, size=31 , prob=0.447) 
xmax <- qbinom(.9995, size=31 , prob=0.447) 
.x <- xmin:xmax 
plot(.x, dbinom(.x, size=31, prob=0.447), xlab="Number of Successes", ylab="Probability Mass",    main="Binomial Dist'n: Trials = 31, Prob of success = 0.447", type="h") 
points(.x, dbinom(.x, size=31, prob=0.447), pch=16) 
abline( h = 0, lty = 2, col = "grey" ) 
remove(.x, xmin, xmax)


###################################################
### code chunk number 179: IPSUR.Rnw:7613-7622
###################################################
xmin <- qbinom(.0005, size=31 , prob=0.447) 
xmax <- qbinom(.9995, size=31 , prob=0.447) 
.x <- xmin:xmax 
plot( stepfun(.x, pbinom((xmin-1):xmax, size=31, prob=0.447)), verticals=FALSE, do.p=FALSE, xlab="Number of Successes", ylab="Cumulative Probability", main="Binomial Dist'n: Trials = 31, Prob of success = 0.447") 
points( .x, pbinom(xmin:xmax, size=31, prob=0.447), pch = 16, cex=1.2 ) 
points( .x, pbinom((xmin-1):(xmax-1), size=31, prob=0.447), pch = 1,    cex=1.2 ) 
abline( h = 1, lty = 2, col = "grey" ) 
abline( h = 0, lty = 2, col = "grey" ) 
remove(.x, xmin, xmax) 


###################################################
### code chunk number 180: IPSUR.Rnw:7629-7630
###################################################
dbinom(17, size = 31, prob = 0.447)


###################################################
### code chunk number 181: IPSUR.Rnw:7636-7637
###################################################
pbinom(13, size = 31, prob = 0.447)


###################################################
### code chunk number 182: IPSUR.Rnw:7643-7644
###################################################
pbinom(11, size = 31, prob = 0.447, lower.tail = FALSE)


###################################################
### code chunk number 183: IPSUR.Rnw:7650-7651
###################################################
pbinom(14, size = 31, prob = 0.447, lower.tail = FALSE)


###################################################
### code chunk number 184: IPSUR.Rnw:7657-7659
###################################################
sum(dbinom(16:19, size = 31, prob = 0.447))
diff(pbinom(c(19,15), size = 31, prob = 0.447, lower.tail = FALSE))


###################################################
### code chunk number 185: IPSUR.Rnw:7665-7668
###################################################
library(distrEx)
X = Binom(size = 31, prob = 0.447)
E(X)


###################################################
### code chunk number 186: IPSUR.Rnw:7674-7675
###################################################
var(X)


###################################################
### code chunk number 187: IPSUR.Rnw:7681-7682
###################################################
sd(X)


###################################################
### code chunk number 188: IPSUR.Rnw:7688-7689
###################################################
E(4*X + 51.324)


###################################################
### code chunk number 189: IPSUR.Rnw:7694-7695
###################################################
rnorm(1)


###################################################
### code chunk number 190: IPSUR.Rnw:7946-7948
###################################################
f <- function(x) 3*x^2
integrate(f, lower = 0.14, upper = 0.71)


###################################################
### code chunk number 191: IPSUR.Rnw:7963-7965
###################################################
g <- function(x) 3/x^3
integrate(g, lower = 1, upper = Inf)


###################################################
### code chunk number 192: IPSUR.Rnw:7980-7984
###################################################
library(distr)
f <- function(x) 3*x^2
X <- AbscontDistribution(d = f, low1 = 0, up1 = 1)
p(X)(0.71) - p(X)(0.14)


###################################################
### code chunk number 193: IPSUR.Rnw:7991-7995
###################################################
library(distrEx)
E(X)
var(X)
3/80


###################################################
### code chunk number 194: IPSUR.Rnw:8110-8111
###################################################
pnorm(1:3)-pnorm(-(1:3))


###################################################
### code chunk number 195: IPSUR.Rnw:8150-8152
###################################################
g <- function(x) pnorm(x, mean = 100, sd = 15) - 0.99
uniroot(g, interval = c(130, 145))


###################################################
### code chunk number 196: IPSUR.Rnw:8154-8155
###################################################
temp <- round(uniroot(g, interval = c(130, 145))$root, 4)


###################################################
### code chunk number 197: IPSUR.Rnw:8225-8226
###################################################
qnorm(0.99, mean = 100, sd = 15)


###################################################
### code chunk number 198: IPSUR.Rnw:8236-8237
###################################################
qnorm(c(0.025, 0.01, 0.005), lower.tail = FALSE)


###################################################
### code chunk number 199: IPSUR.Rnw:8409-8413
###################################################
library(distr)
X <- Norm(mean = 0, sd = 1)
Y <- 4 - 3*X
Y


###################################################
### code chunk number 200: IPSUR.Rnw:8425-8427
###################################################
Y <- exp(X)
Y


###################################################
### code chunk number 201: IPSUR.Rnw:8452-8454
###################################################
W <- sin(exp(X) + 27)
W


###################################################
### code chunk number 202: IPSUR.Rnw:8465-8468
###################################################
p(W)(0.5)
W <- sin(exp(X) + 27)
p(W)(0.5)


###################################################
### code chunk number 203: IPSUR.Rnw:8590-8593 (eval = FALSE)
###################################################
## curve(dchisq(x, df = 3), from = 0, to = 20, ylab = "y")
## ind <- c(4, 5, 10, 15)
## for (i in ind) curve(dchisq(x, df = i), 0, 20, add = TRUE)


###################################################
### code chunk number 204: IPSUR.Rnw:8599-8602
###################################################
curve(dchisq(x, df = 3), from = 0, to = 20, ylab = "y")
ind <- c(4, 5, 10, 15)
for (i in ind) curve(dchisq(x, df = i), 0, 20, add = TRUE)


###################################################
### code chunk number 205: IPSUR.Rnw:8762-8764
###################################################
library(actuar)
mgamma(1:4, shape = 13, rate = 1)


###################################################
### code chunk number 206: IPSUR.Rnw:8769-8770
###################################################
plot(function(x){mgfgamma(x, shape = 13, rate = 1)}, from=-0.1, to=0.1, ylab = "gamma mgf")


###################################################
### code chunk number 207: IPSUR.Rnw:8776-8777
###################################################
plot(function(x){mgfgamma(x, shape = 13, rate = 1)}, from=-0.1, to=0.1, ylab = "gamma mgf")


###################################################
### code chunk number 208: IPSUR.Rnw:8808-8809
###################################################
rnorm(1)


###################################################
### code chunk number 209: IPSUR.Rnw:8835-8836
###################################################
pnorm(2.64, lower.tail = FALSE)


###################################################
### code chunk number 210: IPSUR.Rnw:8842-8843
###################################################
pnorm(0.87) - 1/2


###################################################
### code chunk number 211: IPSUR.Rnw:8849-8850
###################################################
2 * pnorm(-1.39)


###################################################
### code chunk number 212: IPSUR.Rnw:9164-9168
###################################################
S <- rolldie(2, makespace = TRUE)
S <- addrv(S, FUN = max, invars = c("X1","X2"), name = "U")
S <- addrv(S, FUN = sum, invars = c("X1","X2"), name = "V")
head(S)


###################################################
### code chunk number 213: IPSUR.Rnw:9182-9184
###################################################
UV <- marginal(S, vars = c("U", "V"))
head(UV)


###################################################
### code chunk number 214: IPSUR.Rnw:9192-9193
###################################################
xtabs(round(probs,3) ~ U + V, data = UV)


###################################################
### code chunk number 215: IPSUR.Rnw:9200-9202
###################################################
marginal(UV, vars = "U")
head(marginal(UV, vars = "V"))


###################################################
### code chunk number 216: IPSUR.Rnw:9210-9213
###################################################
temp <- xtabs(probs ~ U + V, data = UV)
rowSums(temp)
colSums(temp)


###################################################
### code chunk number 217: IPSUR.Rnw:9304-9308
###################################################
Eu <- sum(S$U*S$probs)
Ev <- sum(S$V*S$probs)
Euv <- sum(S$U*S$V*S$probs)
Euv - Eu * Ev


###################################################
### code chunk number 218: IPSUR.Rnw:9710-9715 (eval = FALSE)
###################################################
## library(mvtnorm)
## x <- y <- seq(from = -3, to = 3, length.out = 30)
## f <- function(x,y) dmvnorm(cbind(x,y), mean = c(0,0), sigma = diag(2))
## z <- outer(x, y, FUN = f)
## persp(x, y, z, theta = -30, phi = 30, ticktype = "detailed")


###################################################
### code chunk number 219: IPSUR.Rnw:9724-9729
###################################################
library(mvtnorm)
x <- y <- seq(from = -3, to = 3, length.out = 30)
f <- function(x,y) dmvnorm(cbind(x,y), mean = c(0,0), sigma = diag(2))
z <- outer(x, y, FUN = f)
persp(x, y, z, theta = -30, phi = 30, ticktype = "detailed")


###################################################
### code chunk number 220: IPSUR.Rnw:10045-10052
###################################################
library(combinat)
tmp <- t(xsimplex(3, 6))
p <- apply(tmp, MARGIN = 1, FUN = dmultinom, prob = c(36,27,37))
library(prob)
S <- probspace(tmp, probs = p)
ProbTable <- xtabs(probs ~ X1 + X2, data = S)
round(ProbTable, 3)


###################################################
### code chunk number 221: IPSUR.Rnw:10085-10087
###################################################
library(lattice)
print(cloud(probs ~ X1 + X2, data = S, type = c("p","h"), lwd = 2, pch = 16, cex = 1.5), screen = list(z = 15, x = -70))


###################################################
### code chunk number 222: IPSUR.Rnw:10337-10340 (eval = FALSE)
###################################################
## curve(dt(x, df = 30), from = -3, to = 3, lwd = 3, ylab = "y")
## ind <- c(1, 2, 3, 5, 10)
## for (i in ind) curve(dt(x, df = i), -3, 3, add = TRUE)


###################################################
### code chunk number 223: IPSUR.Rnw:10346-10349
###################################################
curve(dt(x, df = 30), from = -3, to = 3, lwd = 3, ylab = "y")
ind <- c(1, 2, 3, 5, 10)
for (i in ind) curve(dt(x, df = i), -3, 3, add = TRUE)


###################################################
### code chunk number 224: IPSUR.Rnw:10364-10365
###################################################
qt(0.01, df = 23, lower.tail = FALSE)


###################################################
### code chunk number 225: IPSUR.Rnw:10446-10448 (eval = FALSE)
###################################################
## library(TeachingDemos)
## example(clt.examp)


###################################################
### code chunk number 226: IPSUR.Rnw:10453-10455 (eval = FALSE)
###################################################
## library(distrTeach)
## example(illustrateCLT)


###################################################
### code chunk number 227: IPSUR.Rnw:10634-10635
###################################################
iqrs <- replicate(100, IQR(rnorm(100)))


###################################################
### code chunk number 228: IPSUR.Rnw:10640-10641
###################################################
mean(iqrs)    # close to 1


###################################################
### code chunk number 229: IPSUR.Rnw:10646-10647
###################################################
sd(iqrs)


###################################################
### code chunk number 230: IPSUR.Rnw:10655-10656
###################################################
hist(iqrs, breaks = 20)


###################################################
### code chunk number 231: IPSUR.Rnw:10668-10669
###################################################
mads <- replicate(100, mad(rnorm(100)))


###################################################
### code chunk number 232: IPSUR.Rnw:10674-10675
###################################################
mean(mads)    # close to 1.349


###################################################
### code chunk number 233: IPSUR.Rnw:10680-10681
###################################################
sd(mads)


###################################################
### code chunk number 234: IPSUR.Rnw:10689-10690
###################################################
hist(mads, breaks = 20)


###################################################
### code chunk number 235: IPSUR.Rnw:10707-10710
###################################################
k = 1
n = sample(10:30, size=10, replace = TRUE)
mu = round(rnorm(10, mean = 20))


###################################################
### code chunk number 236: IPSUR.Rnw:10821-10822
###################################################
pnorm(43.1, mean = 37, sd = 9, lower.tail = FALSE)


###################################################
### code chunk number 237: IPSUR.Rnw:10926-10934
###################################################
heights = rep(0, 16)
for (j in 7:15) heights[j] <- dhyper(3, m = 7, n = j - 7, k = 4)
plot(6:15, heights[6:15], pch = 16, cex = 1.5, xlab = "number of fish in pond", ylab = "Likelihood")
abline(h = 0)
lines(6:15, heights[6:15], type = "h", lwd = 2, lty = 3)
text(9, heights[9]/6, bquote(hat(F)==.(9)), cex = 2, pos = 4)
lines(9, heights[9], type = "h", lwd = 2)
points(9, 0, pch = 4, lwd = 3, cex = 2) 


###################################################
### code chunk number 238: IPSUR.Rnw:10985-10988 (eval = FALSE)
###################################################
## curve(x^5*(1-x)^2, from = 0, to = 1, xlab = "p", ylab = "L(p)")
## curve(x^4*(1-x)^3, from = 0, to = 1, add = TRUE)
## curve(x^3*(1-x)^4, 0, 1, add = TRUE)


###################################################
### code chunk number 239: IPSUR.Rnw:10994-10997
###################################################
curve(x^5*(1-x)^2, 0, 1, xlab = "p", ylab = "L(p)")
curve(x^4*(1-x)^3, 0, 1, add = TRUE)
curve(x^3*(1-x)^4, 0, 1, add = TRUE)


###################################################
### code chunk number 240: IPSUR.Rnw:11046-11059
###################################################
dat <- rbinom(27, size = 1, prob = 0.3)
like <- function(x){
r <- 1
for (k in 1:27){ r <- r*dbinom(dat[k], size = 1, prob = x)}
return(r)
}
curve(like, from = 0, to = 1, xlab = "parameter space", ylab = "Likelihood", lwd = 3, col = "blue")
abline(h = 0, lwd = 1, lty = 3, col = "grey")
mle <- mean(dat)
mleobj <- like(mle)
lines(mle, mleobj, type = "h", lwd = 2, lty = 3, col = "red")
points(mle, 0, pch = 4, lwd = 2, cex = 2, col = "red")
text(mle, mleobj/6, substitute(hat(theta)==a, list(a=round(mle, 4))), cex = 2, pos = 4)


###################################################
### code chunk number 241: IPSUR.Rnw:11180-11183
###################################################
x <- mtcars$am
L <- function(p,x) prod(dbinom(x, size = 1, prob = p))
optimize(L, interval = c(0,1), x = x, maximum = TRUE)


###################################################
### code chunk number 242: IPSUR.Rnw:11186-11187
###################################################
A <- optimize(L, interval = c(0,1), x = x, maximum = TRUE)


###################################################
### code chunk number 243: IPSUR.Rnw:11206-11208
###################################################
minuslogL <- function(p,x) -sum(dbinom(x, size = 1, prob = p, log = TRUE))
optimize(minuslogL, interval = c(0,1), x = x)


###################################################
### code chunk number 244: IPSUR.Rnw:11231-11234
###################################################
minuslogL <- function(mu, sigma2){
  -sum(dnorm(x, mean = mu, sd = sqrt(sigma2), log = TRUE))
}


###################################################
### code chunk number 245: IPSUR.Rnw:11245-11249
###################################################
x <- PlantGrowth$weight
library(stats4)
MaxLikeEst <- mle(minuslogL, start = list(mu = 5, sigma2 = 0.5))
summary(MaxLikeEst)


###################################################
### code chunk number 246: IPSUR.Rnw:11259-11262
###################################################
mean(x)
var(x)*29/30
sd(x)/sqrt(30)


###################################################
### code chunk number 247: IPSUR.Rnw:11340-11343
###################################################
set.seed(seed + 1)
library(TeachingDemos)
ci.examp()


###################################################
### code chunk number 248: IPSUR.Rnw:11399-11401
###################################################
library(aplpack)
with(PlantGrowth, stem.leaf(weight))


###################################################
### code chunk number 249: IPSUR.Rnw:11413-11416
###################################################
dim(PlantGrowth)   # sample size is first entry
with(PlantGrowth, mean(weight))
qnorm(0.975)


###################################################
### code chunk number 250: IPSUR.Rnw:11431-11433
###################################################
library(TeachingDemos)
plot(z.test(PlantGrowth$weight, stdev = 0.70), "Conf")


###################################################
### code chunk number 251: IPSUR.Rnw:11538-11541
###################################################
library(TeachingDemos)
temp <- with(PlantGrowth, z.test(weight, stdev = 0.7))
temp


###################################################
### code chunk number 252: IPSUR.Rnw:11549-11551 (eval = FALSE)
###################################################
## library(IPSUR)
## plot(temp, "Conf")


###################################################
### code chunk number 253: IPSUR.Rnw:11706-11709
###################################################
library(Hmisc)
binconf(x = 7, n = 25, method = "asymptotic")
binconf(x = 7, n = 25, method = "wilson")


###################################################
### code chunk number 254: IPSUR.Rnw:11717-11719
###################################################
tab <- xtabs(~gender, data = RcmdrTestDrive)
prop.test(rbind(tab), conf.level = 0.95, correct = FALSE)


###################################################
### code chunk number 255: IPSUR.Rnw:11722-11725
###################################################
A <- as.data.frame(Titanic)
library(reshape)
B <- with(A, untable(A, Freq))


###################################################
### code chunk number 256: IPSUR.Rnw:11931-11932
###################################################
dhyper(0, m = 26, n = 26, k = 5)


###################################################
### code chunk number 257: IPSUR.Rnw:12058-12059
###################################################
- qnorm(0.99)


###################################################
### code chunk number 258: IPSUR.Rnw:12067-12070
###################################################
A <- as.data.frame(UCBAdmissions)
head(A)
xtabs(Freq ~ Admit, data = A)


###################################################
### code chunk number 259: IPSUR.Rnw:12075-12077
###################################################
phat <- 1755/(1755 + 2771)
(phat - 0.4)/sqrt(0.4 * 0.6/(1755 + 2771)) 


###################################################
### code chunk number 260: IPSUR.Rnw:12099-12100
###################################################
-qnorm(0.95)


###################################################
### code chunk number 261: IPSUR.Rnw:12150-12151
###################################################
pnorm(-1.680919)


###################################################
### code chunk number 262: IPSUR.Rnw:12183-12184
###################################################
prop.test(1755, 1755 + 2771, p = 0.4, alternative = "less", conf.level = 0.99, correct = FALSE)


###################################################
### code chunk number 263: IPSUR.Rnw:12189-12193 (eval = FALSE)
###################################################
## library(IPSUR)
## library(HH)
## temp <- prop.test(1755, 1755 + 2771, p = 0.4, alternative = "less", conf.level = 0.99, correct = FALSE)
## plot(temp, 'Hypoth')


###################################################
### code chunk number 264: IPSUR.Rnw:12199-12201
###################################################
library(HH)
plot(prop.test(1755, 1755 + 2771, p = 0.4, alternative = "less", conf.level = 0.99, correct = FALSE), 'Hypoth')


###################################################
### code chunk number 265: IPSUR.Rnw:12301-12304
###################################################
x <- rnorm(37, mean = 2, sd = 3)
library(TeachingDemos)
z.test(x, mu = 1, sd = 3, conf.level = 0.90)


###################################################
### code chunk number 266: IPSUR.Rnw:12310-12312
###################################################
library(HH)
plot(prop.test(1755, 1755 + 2771, p = 0.4, alternative = "less", conf.level = 0.99, correct = FALSE), 'Hypoth')


###################################################
### code chunk number 267: IPSUR.Rnw:12330-12332
###################################################
x <- rnorm(13, mean = 2, sd = 3)
t.test(x, mu = 0, conf.level = 0.90, alternative = "greater")


###################################################
### code chunk number 268: IPSUR.Rnw:12365-12367
###################################################
library(TeachingDemos)
sigma.test(women$height, sigma = 8)


###################################################
### code chunk number 269: IPSUR.Rnw:12441-12442
###################################################
t.test(extra ~ group, data = sleep, paired = TRUE)


###################################################
### code chunk number 270: IPSUR.Rnw:12454-12455
###################################################
ks.test(randu$x, "punif")


###################################################
### code chunk number 271: IPSUR.Rnw:12464-12465
###################################################
shapiro.test(women$height)


###################################################
### code chunk number 272: IPSUR.Rnw:12476-12477
###################################################
with(chickwts, by(weight, feed, shapiro.test))


###################################################
### code chunk number 273: IPSUR.Rnw:12482-12483
###################################################
temp <- lm(weight ~ feed, data = chickwts)


###################################################
### code chunk number 274: IPSUR.Rnw:12488-12489
###################################################
anova(temp)


###################################################
### code chunk number 275: IPSUR.Rnw:12497-12506
###################################################
y1 <- rnorm(300, mean = c(2,8,22))
plot(y1, xlim = c(-1,25), ylim = c(0,0.45) , type = "n")
f <- function(x){dnorm(x, mean = 2)}
curve(f, from = -1, to = 5, add = TRUE, lwd = 2)
f <- function(x){dnorm(x, mean = 8)}
curve(f, from = 5, to = 11, add = TRUE, lwd = 2)
f <- function(x){dnorm(x, mean = 22)}
curve(f, from = 19, to = 25, add = TRUE, lwd = 2)
rug(y1)


###################################################
### code chunk number 276: IPSUR.Rnw:12522-12530
###################################################
y2 <- rnorm(300, mean = c(4,4.1,4.3))
hist(y2, 30, prob = TRUE)
f <- function(x){dnorm(x, mean = 4)/3}
curve(f, add = TRUE, lwd = 2)
f <- function(x){dnorm(x, mean = 4.1)/3}
curve(f, add = TRUE, lwd = 2)
f <- function(x){dnorm(x, mean = 4.3)/3}
curve(f, add = TRUE, lwd = 2)


###################################################
### code chunk number 277: IPSUR.Rnw:12540-12546
###################################################
library(HH)
old.omd <- par(omd = c(.05,.88, .05,1))
F.setup(df1 = 5, df2 = 30)
F.curve(df1 = 5, df2 = 30, col='blue')
F.observed(3, df1 = 5, df2 = 30)
par(old.omd)


###################################################
### code chunk number 278: IPSUR.Rnw:12587-12589
###################################################
library(TeachingDemos)
power.examp()


###################################################
### code chunk number 279: IPSUR.Rnw:12744-12761
###################################################
 # open window
plot(c(0,5), c(0,6.5), type = "n", xlab="x", ylab="y")
## the x- and y-axes
abline(h=0, v=0, col = "gray60")
# regression line
abline(a = 2.5, b = 0.5, lwd = 2)
# normal curves
x <- 600:3000/600
y <- dnorm(x, mean = 3, sd = 0.5)
lines(y + 1.0, x)
lines(y + 2.5, x + 0.75)
lines(y + 4.0, x + 1.5)
# pretty it up
abline(v = c(1, 2.5, 4), lty = 2, col = "grey")
segments(1,3, 1+dnorm(0,0,0.5),3, lty = 2, col = "gray")
segments(2.5, 3.75, 2.5+dnorm(0,0,0.5), 3.75, lty = 2, col = "gray")
segments(4,4.5, 4+dnorm(0,0,0.5),4.5, lty = 2, col = "gray")


###################################################
### code chunk number 280: IPSUR.Rnw:12777-12778
###################################################
head(cars)


###################################################
### code chunk number 281: IPSUR.Rnw:12790-12791
###################################################
plot(dist ~ speed, data = cars)


###################################################
### code chunk number 282: IPSUR.Rnw:12805-12806 (eval = FALSE)
###################################################
## plot(dist ~ speed, data = cars)


###################################################
### code chunk number 283: IPSUR.Rnw:12885-12886
###################################################
cars.lm <- lm(dist ~ speed, data = cars)


###################################################
### code chunk number 284: IPSUR.Rnw:12904-12905
###################################################
coef(cars.lm)


###################################################
### code chunk number 285: IPSUR.Rnw:12921-12923
###################################################
plot(dist ~ speed, data = cars, pch = 16)
abline(coef(cars.lm))


###################################################
### code chunk number 286: IPSUR.Rnw:12932-12934 (eval = FALSE)
###################################################
## plot(dist ~ speed, data = cars, pch = 16)
## abline(coef(cars))


###################################################
### code chunk number 287: IPSUR.Rnw:13006-13007
###################################################
cars[5, ]


###################################################
### code chunk number 288: IPSUR.Rnw:13039-13040
###################################################
fitted(cars.lm)[1:5]


###################################################
### code chunk number 289: IPSUR.Rnw:13053-13054
###################################################
predict(cars.lm, newdata = data.frame(speed = c(6, 8, 21)))


###################################################
### code chunk number 290: IPSUR.Rnw:13089-13090
###################################################
residuals(cars.lm)[1:5]


###################################################
### code chunk number 291: IPSUR.Rnw:13111-13113
###################################################
carsumry <- summary(cars.lm)
carsumry$sigma


###################################################
### code chunk number 292: IPSUR.Rnw:13170-13171
###################################################
summary(cars.lm)


###################################################
### code chunk number 293: IPSUR.Rnw:13173-13175
###################################################
A <- round(summary(cars.lm)$coef, 3)
B <- round(confint(cars.lm), 3)


###################################################
### code chunk number 294: IPSUR.Rnw:13189-13190
###################################################
confint(cars.lm)


###################################################
### code chunk number 295: IPSUR.Rnw:13261-13262
###################################################
new <- data.frame(speed = c(5, 6, 21))


###################################################
### code chunk number 296: IPSUR.Rnw:13268-13269
###################################################
predict(cars.lm, newdata = new, interval = "confidence")


###################################################
### code chunk number 297: IPSUR.Rnw:13272-13273
###################################################
carsCI <- round(predict(cars.lm, newdata = new, interval = "confidence"), 2)


###################################################
### code chunk number 298: IPSUR.Rnw:13278-13279
###################################################
predict(cars.lm, newdata = new, interval = "prediction")


###################################################
### code chunk number 299: IPSUR.Rnw:13282-13283
###################################################
carsPI <- round(predict(cars.lm, newdata = new, interval = "prediction"), 2)


###################################################
### code chunk number 300: IPSUR.Rnw:13336-13338
###################################################
library(HH)
print(ci.plot(cars.lm))


###################################################
### code chunk number 301: IPSUR.Rnw:13348-13350 (eval = FALSE)
###################################################
## library(HH)
## ci.plot(cars.lm)


###################################################
### code chunk number 302: IPSUR.Rnw:13394-13395
###################################################
summary(cars.lm)


###################################################
### code chunk number 303: IPSUR.Rnw:13397-13399
###################################################
A <- round(summary(cars.lm)$coef, 3)
B <- round(confint(cars.lm), 3)


###################################################
### code chunk number 304: IPSUR.Rnw:13484-13485
###################################################
anova(cars.lm)


###################################################
### code chunk number 305: IPSUR.Rnw:13498-13499
###################################################
carsumry$r.squared


###################################################
### code chunk number 306: IPSUR.Rnw:13511-13512
###################################################
sqrt(carsumry$r.squared)


###################################################
### code chunk number 307: IPSUR.Rnw:13553-13554
###################################################
anova(cars.lm)


###################################################
### code chunk number 308: IPSUR.Rnw:13610-13611
###################################################
plot(cars.lm, which = 2)


###################################################
### code chunk number 309: IPSUR.Rnw:13647-13648
###################################################
shapiro.test(residuals(cars.lm))


###################################################
### code chunk number 310: IPSUR.Rnw:13698-13699
###################################################
plot(cars.lm, which = 3)


###################################################
### code chunk number 311: IPSUR.Rnw:13733-13735
###################################################
library(lmtest)
bptest(cars.lm)


###################################################
### code chunk number 312: IPSUR.Rnw:13767-13768
###################################################
plot(cars.lm, which = 1)


###################################################
### code chunk number 313: IPSUR.Rnw:13797-13799
###################################################
library(lmtest)
dwtest(cars.lm, alternative = "two.sided")


###################################################
### code chunk number 314: IPSUR.Rnw:13941-13943
###################################################
sres <- rstandard(cars.lm)
sres[1:5]


###################################################
### code chunk number 315: IPSUR.Rnw:13949-13950
###################################################
sres[which(abs(sres) > 2)]


###################################################
### code chunk number 316: IPSUR.Rnw:13958-13960
###################################################
sdelres <- rstudent(cars.lm)
sdelres[1:5]


###################################################
### code chunk number 317: IPSUR.Rnw:13967-13969
###################################################
t0.005 <- qt(0.005, df = 47, lower.tail = FALSE)
sdelres[which(abs(sdelres) > t0.005)]


###################################################
### code chunk number 318: IPSUR.Rnw:13976-13979
###################################################
leverage <- hatvalues(cars.lm)
leverage[1:5]
leverage[which(leverage > 4/50)]


###################################################
### code chunk number 319: IPSUR.Rnw:14020-14022
###################################################
dfb <- dfbetas(cars.lm)
head(dfb)


###################################################
### code chunk number 320: IPSUR.Rnw:14042-14044
###################################################
dff <- dffits(cars.lm)
dff[1:5]


###################################################
### code chunk number 321: IPSUR.Rnw:14080-14082
###################################################
cooksD <- cooks.distance(cars.lm)
cooksD[1:5]


###################################################
### code chunk number 322: IPSUR.Rnw:14090-14091
###################################################
plot(cars.lm, which = 4)


###################################################
### code chunk number 323: IPSUR.Rnw:14108-14110
###################################################
F0.50 <- qf(0.5, df1 = 2, df2 = 48)
cooksD[which(cooksD > F0.50)]


###################################################
### code chunk number 324: IPSUR.Rnw:14123-14124 (eval = FALSE)
###################################################
## influence.measures(cars.lm)


###################################################
### code chunk number 325: IPSUR.Rnw:14134-14137 (eval = FALSE)
###################################################
## par(mfrow = c(2,2))
## plot(cars.lm)
## par(mfrow = c(1,1))


###################################################
### code chunk number 326: IPSUR.Rnw:14148-14151
###################################################
par(mfrow = c(2,2))
plot(cars.lm)
par(mfrow = c(1,1))


###################################################
### code chunk number 327: IPSUR.Rnw:14173-14175 (eval = FALSE)
###################################################
## plot(cars.lm, which = 5)   # std'd resids vs lev plot
## identify(leverage, sres, n = 4)   # identify 4 points


###################################################
### code chunk number 328: IPSUR.Rnw:14297-14298
###################################################
head(trees)


###################################################
### code chunk number 329: IPSUR.Rnw:14311-14313
###################################################
library(lattice)
print(splom(trees))


###################################################
### code chunk number 330: IPSUR.Rnw:14322-14324 (eval = FALSE)
###################################################
## library(lattice)
## splom(trees)


###################################################
### code chunk number 331: IPSUR.Rnw:14374-14378 (eval = FALSE)
###################################################
## library(scatterplot3d)
## s3d <- with(trees, scatterplot3d(Girth, Height, Volume, pch = 16, highlight.3d = TRUE, angle = 60))
## fit <- lm(Volume ~ Girth + Height, data = trees)
## s3d$plane3d(fit)


###################################################
### code chunk number 332: IPSUR.Rnw:14384-14388
###################################################
library(scatterplot3d)
s3d <- with(trees, scatterplot3d(Girth, Height, Volume, pch = 16, highlight.3d = TRUE, angle = 60))
fit <- lm(Volume ~ Girth + Height, data = trees)
s3d$plane3d(fit)


###################################################
### code chunk number 333: IPSUR.Rnw:14461-14463
###################################################
trees.lm <- lm(Volume ~ Girth + Height, data = trees)
trees.lm


###################################################
### code chunk number 334: IPSUR.Rnw:14476-14477
###################################################
head(model.matrix(trees.lm))


###################################################
### code chunk number 335: IPSUR.Rnw:14546-14547
###################################################
fitted(trees.lm)[1:5]


###################################################
### code chunk number 336: IPSUR.Rnw:14563-14564
###################################################
new <- data.frame(Girth = c(9.1, 11.6, 12.5), Height = c(69, 74, 87))


###################################################
### code chunk number 337: IPSUR.Rnw:14569-14570
###################################################
new


###################################################
### code chunk number 338: IPSUR.Rnw:14575-14576
###################################################
predict(trees.lm, newdata = new)


###################################################
### code chunk number 339: IPSUR.Rnw:14579-14580
###################################################
treesFIT <- round(predict(trees.lm, newdata = new), 1)


###################################################
### code chunk number 340: IPSUR.Rnw:14636-14637
###################################################
residuals(trees.lm)[1:5]


###################################################
### code chunk number 341: IPSUR.Rnw:14646-14648
###################################################
treesumry <- summary(trees.lm)
treesumry$sigma


###################################################
### code chunk number 342: IPSUR.Rnw:14704-14705
###################################################
confint(trees.lm)


###################################################
### code chunk number 343: IPSUR.Rnw:14708-14709
###################################################
treesPAR <- round(confint(trees.lm), 1)


###################################################
### code chunk number 344: IPSUR.Rnw:14755-14756
###################################################
new <- data.frame(Girth = c(9.1, 11.6, 12.5), Height = c(69, 74, 87))


###################################################
### code chunk number 345: IPSUR.Rnw:14761-14762
###################################################
predict(trees.lm, newdata = new, interval = "confidence")


###################################################
### code chunk number 346: IPSUR.Rnw:14765-14766
###################################################
treesCI <- round(predict(trees.lm, newdata = new, interval = "confidence"), 1)


###################################################
### code chunk number 347: IPSUR.Rnw:14771-14772
###################################################
predict(trees.lm, newdata = new, interval = "prediction")


###################################################
### code chunk number 348: IPSUR.Rnw:14775-14776
###################################################
treesPI <- round(predict(trees.lm, newdata = new, interval = "prediction"), 1)


###################################################
### code chunk number 349: IPSUR.Rnw:14871-14873
###################################################
treesumry$r.squared
treesumry$adj.r.squared


###################################################
### code chunk number 350: IPSUR.Rnw:14905-14906
###################################################
treesumry$fstatistic


###################################################
### code chunk number 351: IPSUR.Rnw:14956-14957
###################################################
treesumry


###################################################
### code chunk number 352: IPSUR.Rnw:14986-14987
###################################################
plot(Volume ~ Girth, data = trees)


###################################################
### code chunk number 353: IPSUR.Rnw:15073-15075
###################################################
treesquad.lm <- lm(Volume ~ scale(Girth) + I(scale(Girth)^2), data = trees)
summary(treesquad.lm)


###################################################
### code chunk number 354: IPSUR.Rnw:15089-15091 (eval = FALSE)
###################################################
## plot(Volume ~ scale(Girth), data = trees)
## lines(fitted(treesquad.lm) ~ scale(Girth), data = trees)


###################################################
### code chunk number 355: IPSUR.Rnw:15103-15105
###################################################
plot(Volume ~ scale(Girth), data = trees)
lines(fitted(treesquad.lm) ~ scale(Girth), data = trees)


###################################################
### code chunk number 356: IPSUR.Rnw:15126-15128
###################################################
new <- data.frame(Girth = c(9.1, 11.6, 12.5))
predict(treesquad.lm, newdata = new, interval = "prediction")


###################################################
### code chunk number 357: IPSUR.Rnw:15140-15141
###################################################
summary(lm(Volume ~ Girth + I(Girth^2), data = trees))


###################################################
### code chunk number 358: IPSUR.Rnw:15222-15224
###################################################
treesint.lm <- lm(Volume ~ Girth + Height + Girth:Height, data = trees)
summary(treesint.lm)


###################################################
### code chunk number 359: IPSUR.Rnw:15236-15239
###################################################
confint(treesint.lm)
new <- data.frame(Girth = c(9.1, 11.6, 12.5), Height = c(69, 74, 87))
predict(treesint.lm, newdata = new, interval = "prediction")


###################################################
### code chunk number 360: IPSUR.Rnw:15287-15289
###################################################
trees$Tall <- cut(trees$Height, breaks = c(-Inf, 76, Inf), labels = c("no","yes"))
trees$Tall[1:5]


###################################################
### code chunk number 361: IPSUR.Rnw:15329-15330
###################################################
class(trees$Tall)


###################################################
### code chunk number 362: IPSUR.Rnw:15338-15340
###################################################
treesdummy.lm <- lm(Volume ~ Girth + Tall, data = trees)
summary(treesdummy.lm)


###################################################
### code chunk number 363: IPSUR.Rnw:15388-15396 (eval = FALSE)
###################################################
## treesTall <- split(trees, trees$Tall)
## treesTall[["yes"]]$Fit <- predict(treesdummy.lm, treesTall[["yes"]])
## treesTall[["no"]]$Fit <- predict(treesdummy.lm, treesTall[["no"]])
## plot(Volume ~ Girth, data = trees, type = "n")
## points(Volume ~ Girth, data = treesTall[["yes"]], pch = 1)
## points(Volume ~ Girth, data = treesTall[["no"]], pch = 2)
## lines(Fit ~ Girth, data = treesTall[["yes"]])
## lines(Fit ~ Girth, data = treesTall[["no"]])


###################################################
### code chunk number 364: IPSUR.Rnw:15402-15410
###################################################
treesTall <- split(trees, trees$Tall)
treesTall[["yes"]]$Fit <- predict(treesdummy.lm, treesTall[["yes"]])
treesTall[["no"]]$Fit <- predict(treesdummy.lm, treesTall[["no"]])
plot(Volume ~ Girth, data = trees, type = "n")
points(Volume ~ Girth, data = treesTall[["yes"]], pch = 1)
points(Volume ~ Girth, data = treesTall[["no"]], pch = 2)
lines(Fit ~ Girth, data = treesTall[["yes"]])
lines(Fit ~ Girth, data = treesTall[["no"]])


###################################################
### code chunk number 365: IPSUR.Rnw:15503-15505
###################################################
treesfull.lm <- lm(Volume ~ Girth + I(Girth^2) + Height + I(Height^2), data = trees)
summary(treesfull.lm)


###################################################
### code chunk number 366: IPSUR.Rnw:15520-15521
###################################################
treesreduced.lm <- lm(Volume ~ -1 + Girth + I(Girth^2), data = trees)


###################################################
### code chunk number 367: IPSUR.Rnw:15528-15529
###################################################
anova(treesreduced.lm, treesfull.lm)


###################################################
### code chunk number 368: IPSUR.Rnw:15539-15541
###################################################
treesreduced2.lm <- lm(Volume ~ Girth + I(Girth^2) + Height, data = trees)
anova(treesreduced2.lm, treesfull.lm)


###################################################
### code chunk number 369: IPSUR.Rnw:15636-15638
###################################################
treesNonlin.lm <- lm(log(Volume) ~ log(Girth) + log(Height), data = trees)
summary(treesNonlin.lm)


###################################################
### code chunk number 370: IPSUR.Rnw:15650-15651
###################################################
exp(confint(treesNonlin.lm))


###################################################
### code chunk number 371: IPSUR.Rnw:15659-15661
###################################################
new <- data.frame(Girth = c(9.1, 11.6, 12.5), Height = c(69, 74, 87))
exp(predict(treesNonlin.lm, newdata = new, interval = "confidence"))


###################################################
### code chunk number 372: IPSUR.Rnw:15700-15708
###################################################
# fake data 
set.seed(1) 
x <- seq(from = 0, to = 1000, length.out = 200) 
y <- 1 + 2*(sin((2*pi*x/360) - 3))^2 + rnorm(200, sd = 2)
plot(x, y)
acc.nls <- nls(y ~ a + b*(sin((2*pi*x/360) - c))^2, start = list(a = 0.9, b = 2.3, c = 2.9))
summary(acc.nls)
#plot(x, fitted(acc.nls))


###################################################
### code chunk number 373: IPSUR.Rnw:15909-15912
###################################################
srs <- rnorm(25, mean = 3)
resamps <- replicate(1000, sample(srs, 25, TRUE), simplify = FALSE)
xbarstar <- sapply(resamps, mean, simplify = TRUE)


###################################################
### code chunk number 374: IPSUR.Rnw:15918-15920
###################################################
hist(xbarstar, breaks = 40, prob = TRUE)
curve(dnorm(x, 3, 0.2), add = TRUE)


###################################################
### code chunk number 375: IPSUR.Rnw:15943-15945 (eval = FALSE)
###################################################
## hist(xbarstar, breaks = 40, prob = TRUE)
## curve(dnorm(x, 3, 0.2), add = TRUE)  # overlay true normal density


###################################################
### code chunk number 376: IPSUR.Rnw:15959-15962
###################################################
mean(xbarstar)
mean(srs)
mean(xbarstar) - mean(srs)


###################################################
### code chunk number 377: IPSUR.Rnw:15981-15982
###################################################
sd(xbarstar)


###################################################
### code chunk number 378: IPSUR.Rnw:16012-16015
###################################################
resamps <- replicate(1000, sample(rivers, 141, TRUE), simplify = FALSE)
medstar <- sapply(resamps, median, simplify = TRUE)
sd(medstar)


###################################################
### code chunk number 379: IPSUR.Rnw:16021-16022
###################################################
hist(medstar, breaks = 40, prob = TRUE)


###################################################
### code chunk number 380: IPSUR.Rnw:16035-16036 (eval = FALSE)
###################################################
## hist(medstar, breaks = 40, prob = TRUE)


###################################################
### code chunk number 381: IPSUR.Rnw:16039-16042
###################################################
median(rivers)
mean(medstar)
mean(medstar) - median(rivers)


###################################################
### code chunk number 382: IPSUR.Rnw:16067-16070
###################################################
library(boot)
mean_fun <- function(x, indices) mean(x[indices])
boot(data = srs, statistic = mean_fun, R = 1000)


###################################################
### code chunk number 383: IPSUR.Rnw:16075-16077
###################################################
median_fun <- function(x, indices) median(x[indices])
boot(data = rivers, statistic = median_fun, R = 1000)


###################################################
### code chunk number 384: IPSUR.Rnw:16144-16149
###################################################
btsamps <- replicate(2000, sample(stack.loss, 21, TRUE), simplify = FALSE)
thetast <- sapply(btsamps, median, simplify = TRUE)
mean(thetast)
median(stack.loss)
quantile(thetast, c(0.025, 0.975))


###################################################
### code chunk number 385: IPSUR.Rnw:16159-16163
###################################################
library(boot)
med_fun <- function(x, ind) median(x[ind])
med_boot <- boot(stack.loss, med_fun, R = 2000)
boot.ci(med_boot, type = c("perc", "norm", "bca"))


###################################################
### code chunk number 386: IPSUR.Rnw:16312-16314
###################################################
library(coin)
oneway_test(len ~ supp, data = ToothGrowth)


###################################################
### code chunk number 387: IPSUR.Rnw:16326-16327
###################################################
t.test(len ~ supp, data = ToothGrowth, alt = "greater", var.equal = TRUE)


###################################################
### code chunk number 388: IPSUR.Rnw:16330-16332
###################################################
A <- show(oneway_test(len ~ supp, data = ToothGrowth))
B <- t.test(len ~ supp, data = ToothGrowth, alt = "greater", var.equal = TRUE)


###################################################
### code chunk number 389: IPSUR.Rnw:16405-16408 (eval = FALSE)
###################################################
## install.packages("IPSUR", repos="http://R-Forge.R-project.org")
## library(IPSUR)
## read(IPSUR)


###################################################
### code chunk number 390: IPSUR.Rnw:16419-16422 (eval = FALSE)
###################################################
## install.packages("IPSUR", repos="http://R-Forge.R-project.org")
## library(IPSUR)
## read(IPSUR)


###################################################
### code chunk number 391: IPSUR.Rnw:16433-16436 (eval = FALSE)
###################################################
## install.packages("IPSUR", repos="http://R-Forge.R-project.org")
## library(IPSUR)
## read(IPSUR)


###################################################
### code chunk number 392: IPSUR.Rnw:16448-16449
###################################################
sessionInfo()


###################################################
### code chunk number 393: IPSUR.Rnw:16967-16968
###################################################
x <- c(3, 5, 9)


###################################################
### code chunk number 394: IPSUR.Rnw:16975-16976
###################################################
y <- c(3, "5", TRUE)


###################################################
### code chunk number 395: IPSUR.Rnw:16995-16996
###################################################
matrix(letters[1:6], nrow = 2, ncol = 3)


###################################################
### code chunk number 396: IPSUR.Rnw:17003-17004
###################################################
matrix(letters[1:6], nrow = 2, ncol = 3, byrow = TRUE)


###################################################
### code chunk number 397: IPSUR.Rnw:17012-17013
###################################################
matrix(c(1,"2",NA, FALSE), nrow = 2, ncol = 3)


###################################################
### code chunk number 398: IPSUR.Rnw:17023-17027
###################################################
A <- matrix(1:6, 2, 3)
B <- matrix(2:7, 2, 3)
A + B
A * B


###################################################
### code chunk number 399: IPSUR.Rnw:17038-17040
###################################################
try(A * B)     # an error
A %*% t(B)     # this is alright


###################################################
### code chunk number 400: IPSUR.Rnw:17046-17047
###################################################
solve(A %*% t(B))     # input matrix must be square


###################################################
### code chunk number 401: IPSUR.Rnw:17054-17055
###################################################
array(LETTERS[1:24], dim = c(3,4,2))


###################################################
### code chunk number 402: IPSUR.Rnw:17070-17075
###################################################
x <- c(1.3, 5.2, 6)
y <- letters[1:3]
z <- c(TRUE, FALSE, TRUE)
A <- data.frame(x, y, z)
A


###################################################
### code chunk number 403: IPSUR.Rnw:17083-17085
###################################################
names(A) <- c("Fred","Mary","Sue")
A


###################################################
### code chunk number 404: IPSUR.Rnw:17112-17114
###################################################
A <- as.data.frame(Titanic)
head(A)


###################################################
### code chunk number 405: IPSUR.Rnw:17132-17135
###################################################
library(reshape)
B <- with(A, untable(A, Freq))
head(B)


###################################################
### code chunk number 406: IPSUR.Rnw:17155-17158
###################################################
C <- B[, -5]
rownames(C) <- 1:dim(C)[1]
head(C)


###################################################
### code chunk number 407: IPSUR.Rnw:17171-17175
###################################################
tab <- matrix(1:6, nrow = 2, ncol = 3)
rownames(tab) <- c('first', 'second')
colnames(tab) <- c('A', 'B', 'C')
tab  # Counts


###################################################
### code chunk number 408: IPSUR.Rnw:17186-17192
###################################################
p <- c("milk","tea")
g <- c("milk","tea")
catgs <- expand.grid(poured = p, guessed = g)
cnts <- c(3, 1, 1, 3)
D <- cbind(catgs, count = cnts)
xtabs(count ~ poured + guessed, data = D)


###################################################
### code chunk number 409: IPSUR.Rnw:17276-17278 (eval = FALSE)
###################################################
## library(foreign)
## read.spss("foo.sav")


###################################################
### code chunk number 410: IPSUR.Rnw:17326-17328
###################################################
Tmp <- Puromycin[order(Puromycin$conc), ]
head(Tmp)


###################################################
### code chunk number 411: IPSUR.Rnw:17333-17334 (eval = FALSE)
###################################################
## with(Puromycin, Puromycin[order(conc), ])


###################################################
### code chunk number 412: IPSUR.Rnw:17341-17342 (eval = FALSE)
###################################################
## with(Puromycin, Puromycin[order(state, conc), ])


###################################################
### code chunk number 413: IPSUR.Rnw:17348-17350
###################################################
Tmp <- with(Puromycin, Puromycin[order(-conc), ])
head(Tmp)


###################################################
### code chunk number 414: IPSUR.Rnw:17358-17360
###################################################
Tmp <- with(Puromycin, Puromycin[order(-xtfrm(state)), ])
head(Tmp)


###################################################
### code chunk number 415: IPSUR.Rnw:18172-18174 (eval = FALSE)
###################################################
## library(odfWeave)
## odfWeave(file = "infile.odt", dest = "outfile.odt")


###################################################
### code chunk number 416: IPSUR.Rnw:18237-18239
###################################################
library(Hmisc)
summary(cbind(Sepal.Length, Sepal.Width) ~ Species, data = iris)


###################################################
### code chunk number 417: IPSUR.Rnw:18260-18261
###################################################
set.seed(095259)


###################################################
### code chunk number 418: IPSUR.Rnw:18310-18312
###################################################
options(digits = 16)
runif(1)


###################################################
### code chunk number 419: IPSUR.Rnw:18762-18767
###################################################
rm(.Random.seed)
# try(dir.create("../../data"), silent = TRUE)
# save.image(file = "../../data/IPSUR.RData")
# tools::resaveRdaFiles('../../data', compress = 'xz')
# Stangle(file="IPSUR.Rnw", output="../IPSUR.R", annotate=TRUE)


