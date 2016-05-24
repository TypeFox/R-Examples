
set.seed(290875)
library("party")
library("survival")

### get rid of the NAMESPACE
attach(asNamespace("party"))

### check nominal level printing
set.seed(290875)
x <- gl(5, 50)
df <- data.frame(y = c(rnorm(50, 0), rnorm(50, 1), rnorm(50, 2), rnorm(50, 3), rnorm(50, 4)), 
                 x = x, z = rnorm(250))
ctree(y ~ x, data = df)

### check asymptotic vs. MonteCarlo, especially categorical splits after
### MonteCarlo resampling
a <- ctree(y ~ x + z, data = df, control = ctree_control(stump = TRUE))
b <- ctree(y ~ x + z, data = df,
           control =  ctree_control(testtype = "Monte", stump = TRUE))
stopifnot(isequal(a@tree$psplit, b@tree$psplit))  
stopifnot(isequal(a@tree$criterion$statistic, b@tree$criterion$statistic))

### we did not check the hyper parameters
try(cforest_control(minsplit = -1))
try(cforest_control(ntree = -1))
try(cforest_control(maxdepth = -1))
try(cforest_control(nresample = 10))

### NA handling for factors and in random forest
### more than one (ordinal) response variable
xo <- ordered(x)
x[sample(1:length(x), 10)] <- NA
cforest(y + xo ~ x + z, data = df, 
        control = cforest_unbiased(ntree = 50))

### make sure minsplit is OK in the presence of missing values
### spotted by Han Lee <Han.Lee@GeodeCapital.com>
load("t1.RData")
tr <- try(ctree(p ~., data = t1))
stopifnot(!inherits(tr, "try-error"))

### make sure number of surrogate splits exceeds number of inputs by 1
### spotted by Henric Nilsson <henric.nilsson@phadia.com>
airq <- subset(airquality, !is.na(Ozone))
tr <- try(ctree(Ozone ~ Wind, data = airq,
          controls = ctree_control(maxsurrogate = 3)))
stopifnot(inherits(tr, "try-error"))

### ctree() used only the first of a multivariate response
### spotted by Henric Nilsson <henric.nilsson@phadia.com>
airq <- subset(airquality, complete.cases(Ozone, Solar.R))
airOzoSol1 <- ctree(Ozone + Solar.R ~ Wind + Temp + Month + Day,
                    data = airq)
airOzoSol2 <- ctree(Solar.R + Ozone ~ Wind + Temp + Month + Day,
                    data = airq)
stopifnot(isequal(airOzoSol1@where, airOzoSol2@where))

### one variable with all values missing
dat <- data.frame(y = rnorm(100), x1 = runif(100), x2 = rep(NA, 100))
ctree(y ~ x1 + x2, data = dat)

### one factor with only one level
dat$x2 <- factor(rep(0, 100))
try(ctree(y ~ x1 + x2, data = dat))

### weights for sampling without replacement for cforest
### spotted by Carolin Strobl <carolin.strol@stat.uni-muenchen.de>
airq <- subset(airquality, !is.na(Ozone))
cctrl <- cforest_control(replace = FALSE, fraction = 0.5)
n <- nrow(airq)
w <- double(n)


if (FALSE) {
### forest objects have weights remove in 0.9-13

### case weights
x <- runif(w)
w[x > 0.5] <- 1
w[x > 0.9] <- 2

rf <- cforest(Ozone ~ .,data = airq, weights = w, control = cctrl)
rfw <- sapply(rf@ensemble, function(x) x[[2]])
stopifnot(all(colSums(rfw) == ceiling(sum(w) / 2)))
stopifnot(max(abs(rfw[w == 0,])) == 0)

### real weights
w <- runif(n)
w[1:10] <- 0
rf <- cforest(Ozone ~ .,data = airq, weights = w, control = cctrl)
rfw <- sapply(rf@ensemble, function(x) x[[2]])
stopifnot(all(colSums(rfw) == ceiling(sum(w > 0) / 2)))
stopifnot(max(abs(rfw[w == 0,])) == 0)
}

### cforest with multivariate response
df <- data.frame(y1 = rnorm(100), y2 = rnorm(100), x1 = runif(100), x2 = runif(100))
df$y1[df$x1 < 0.5] <- df$y1[df$x1 < 0.5] + 1
cf <- cforest(y1 + y2 ~ x1 + x2, data = df)
pr <- predict(cf)
stopifnot(length(pr) == nrow(df) || lengthl(pr[[1]]) != 2)

### varimp with ordered response
### spotted by Max Kuhn <Max.Kuhn@pfizer.com>
data("mammoexp", package = "TH.data")
test <- cforest(ME ~ ., data = mammoexp, control = cforest_unbiased(ntree = 50))
stopifnot(sum(abs(varimp(test))) > 0)

### missing values in factors lead to segfaults on 64 bit systems
### spotted by Carolin Strobl <carolin.strobl@lme.de>
y <- rnorm(100)
x <- gl(2, 50)
z <- gl(2, 50)[sample(1:100)]
y <- y + (x == "1") * 3
xNA <- x
xNA[1:2] <- NA
ctree(y ~ xNA )


y <- rnorm(100)
x <- y + rnorm(100, sd = 0.1)

tmp <- data.frame(x, y)

x[sample(1:100)[1:10]] <- NA

ct1 <- ctree(y ~ x, data = tmp)
ct2 <- ctree(y ~ x, data = tmp[complete.cases(tmp),])
w <- as.double(complete.cases(tmp))
ct3 <- ctree(y ~ x, data = tmp, weights = w)

xx <- data.frame(x = rnorm(100))
t1 <- max(abs(predict(ct2, newdata = xx) - predict(ct3, newdata = xx))) == 0
t2 <- nterminal(ct1@tree) == nterminal(ct2@tree)
t3 <- nterminal(ct3@tree) == nterminal(ct1@tree)
t4 <- all.equal(ct2@tree$psplit, ct1@tree$psplit)
stopifnot(t1 && t2 && t3 && t4)

y <- rnorm(100)
x <- cut(y, c(-Inf, -1, 0, 1, Inf))

tmp <- data.frame(x, y)

x[sample(1:100)[1:10]] <- NA

ct1 <- ctree(y ~ x, data = tmp)
ct2 <- ctree(y ~ x, data = tmp[complete.cases(tmp),])
w <- as.double(complete.cases(tmp))
ct3 <- ctree(y ~ x, data = tmp, weights = w)

stopifnot(all.equal(ct2@tree$psplit, ct1@tree$psplit))
stopifnot(all.equal(ct2@tree$psplit, ct3@tree$psplit))

### predictions for obs with zero weights
### spotted by Mark Difford <mark_difford@yahoo.co.uk>
airq <- subset(airquality, !is.na(Ozone))
w <- rep(1, nrow(airq))
w[1:5] <- 0

ctw <- ctree(Ozone ~ ., data = airq, weights = w)
stopifnot(all.equal(predict(ctw)[1:5], predict(ctw, newdata = airq)[1:5]))
rfw <- cforest(Ozone ~ ., data = airq, weights = w)
stopifnot(all.equal(predict(rfw)[1:5], predict(rfw, newdata = airq)[1:5]))

### more surrogate splits than available requested
### spotted by Henric Nilsson <henric.nilsson@sorch.se>
airq <- data.frame(airq,
                    x1 = factor(ifelse(runif(nrow(airq)) < 0.5, 0, 1)),
                    x2 = factor(ifelse(runif(nrow(airq)) < 0.5, 0, 1)),
                    x3 = factor(ifelse(runif(nrow(airq)) < 0.5, 0, 1)))

foo <- function(nm) 
    ctree(Ozone ~ ., data = airq,
          controls = ctree_control(maxsurrogate = nm))
foo(4)
try(foo(5))
try(foo(6))

### variance = 0 due to constant variables
### spotted by Sebastian Wietzke <Sebastian.Wietzke@axa.de>
v <- rep(0,20)
w <- rep(0,20)
x <- 1:20
y <- rep(1,20)
z <- c(4,5,8,2,6,1,3,6,8,2,5,8,9,3,5,8,9,4,6,8)
tmp <- ctree(z ~ v+w+x+y,controls = ctree_control(mincriterion = 0.80,
             minsplit = 2, minbucket = 1, testtype = "Univariate", teststat = "quad"))
stopifnot(all(tmp@tree$criterion$criterion[c(1,2,4)] == 0))

### optimal split in last observation lead to selection of suboptimal split
data("GlaucomaM", package = "TH.data")
tmp <- subset(GlaucomaM, vari <= 0.059)
weights <- rep(1.0, nrow(tmp))
stopifnot(all.equal(Split(tmp$vasg, tmp$Class, weights, 
                    ctree_control()@splitctrl)[[1]], 0.066))

### model.matrix.survReg was missing from modeltools
data("GBSG2", package = "TH.data")
nloglik <- function(x) -logLik(x)
GBSG2$time <- GBSG2$time/365
mobGBSG2 <- mob(Surv(time, cens) ~ horTh + pnodes | progrec + menostat +
  estrec + menostat + age + tsize + tgrade, data = GBSG2, model = survReg,
  control = mob_control(objfun = nloglik, minsplit = 40))
plot(mobGBSG2, terminal = node_scatterplot, tp_args = list(yscale = c(-0.1, 11)))

### factors were evaluated for surrogate splits
data("Ozone", package = "mlbench")
Ozone$V2 <- ordered(Ozone$V2)
Ozone <- subset(Ozone, !is.na(V4))
rf <- cforest(V4 ~ ., data = Ozone, control = cforest_unbiased(maxsurrogate = 7))

### scores for response
### spotted and fixed by Silke Janitza <janitza@ibe.med.uni-muenchen.de>
tmp <- data.frame(y = gl(3, 10, ordered = TRUE), x = gl(3, 10, ordered = TRUE))
ct <- ctree(y ~ x, data = tmp, scores = list(y = c(0, 10, 11), x = c(1, 2, 5)))
stopifnot(isTRUE(all.equal(ct@responses@scores, list(y = c(0, 10, 11)))))

### deal with empty levels for teststat = "quad" by
### removing elements of the teststatistic with zero variance  
### reported by Wei-Yin Loh <loh@stat.wisc.edu>
tdata <-
structure(list(ytrain = structure(c(3L, 7L, 3L, 2L, 1L, 6L, 2L, 
1L, 1L, 2L, 1L, 2L, 3L, 3L, 2L, 1L, 2L, 6L, 2L, 4L, 6L, 1L, 2L, 
3L, 7L, 6L, 4L, 6L, 2L, 2L, 1L, 2L, 6L, 1L, 7L, 1L, 3L, 6L, 2L, 
1L, 7L, 2L, 7L, 2L, 3L, 2L, 1L, 1L, 3L, 1L, 6L, 2L, 2L, 2L, 2L, 
2L, 1L, 1L, 6L, 6L, 7L, 2L, 2L, 2L, 2L, 2L, 1L, 3L, 6L, 5L, 1L, 
1L, 4L, 7L, 2L, 3L, 3L, 3L, 1L, 8L, 1L, 6L, 2L, 8L, 3L, 4L, 6L, 
2L, 7L, 3L, 6L, 6L, 1L, 1L, 2L, 6L, 3L, 3L, 1L, 2L, 3L, 1L, 2L, 
7L, 2L, 3L, 6L, 2L, 5L, 2L, 2L, 2L, 1L, 3L, 3L, 7L, 3L, 2L, 3L, 
3L, 1L, 6L, 1L, 1L, 1L, 7L, 1L, 3L, 7L, 6L, 1L, 3L, 3L, 6L, 4L, 
2L, 3L, 2L, 8L, 3L, 4L, 2L, 2L, 2L, 3L, 2L, 2L, 2L, 3L, 4L, 6L, 
4L, 8L, 2L, 2L, 3L, 3L, 2L, 3L, 6L, 2L, 1L, 2L, 2L, 7L, 2L, 1L, 
1L, 7L, 2L, 7L, 6L, 6L, 6L), .Label = c("0", "1", "2", "3", "4", 
"5", "6", "7"), class = "factor"), landmass = c(5L, 3L, 4L, 6L, 
3L, 4L, 1L, 2L, 2L, 6L, 3L, 1L, 5L, 5L, 1L, 3L, 1L, 4L, 1L, 5L, 
4L, 2L, 1L, 5L, 3L, 4L, 5L, 4L, 4L, 1L, 4L, 1L, 4L, 2L, 5L, 2L, 
4L, 4L, 6L, 1L, 1L, 3L, 3L, 3L, 4L, 1L, 1L, 2L, 4L, 1L, 4L, 4L, 
3L, 2L, 6L, 3L, 3L, 2L, 4L, 4L, 3L, 3L, 3L, 3L, 1L, 6L, 1L, 4L, 
4L, 2L, 1L, 1L, 5L, 3L, 3L, 6L, 5L, 5L, 3L, 5L, 3L, 4L, 1L, 5L, 
5L, 5L, 4L, 6L, 5L, 5L, 4L, 4L, 3L, 3L, 4L, 4L, 5L, 5L, 3L, 6L, 
4L, 1L, 6L, 5L, 1L, 4L, 4L, 6L, 5L, 3L, 1L, 6L, 1L, 4L, 4L, 5L, 
5L, 3L, 5L, 5L, 2L, 6L, 2L, 2L, 6L, 3L, 1L, 5L, 3L, 4L, 4L, 5L, 
4L, 4L, 5L, 6L, 4L, 4L, 5L, 5L, 5L, 1L, 1L, 1L, 4L, 2L, 3L, 3L, 
5L, 5L, 4L, 5L, 4L, 6L, 2L, 4L, 5L, 1L, 5L, 4L, 3L, 2L, 1L, 1L, 
5L, 6L, 3L, 2L, 5L, 6L, 3L, 4L, 4L, 4L), zone = c(1L, 1L, 1L, 
3L, 1L, 2L, 4L, 3L, 3L, 2L, 1L, 4L, 1L, 1L, 4L, 1L, 4L, 1L, 4L, 
1L, 2L, 3L, 4L, 1L, 1L, 4L, 1L, 2L, 1L, 4L, 4L, 4L, 1L, 3L, 1L, 
4L, 2L, 2L, 3L, 4L, 4L, 1L, 1L, 1L, 1L, 4L, 4L, 3L, 1L, 4L, 1L, 
1L, 4L, 3L, 2L, 1L, 1L, 4L, 2L, 4L, 1L, 1L, 4L, 1L, 4L, 1L, 4L, 
4L, 4L, 4L, 4L, 4L, 1L, 1L, 4L, 2L, 1L, 1L, 4L, 1L, 1L, 4L, 4L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 4L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 
1L, 4L, 4L, 1L, 1L, 4L, 4L, 2L, 2L, 1L, 1L, 4L, 2L, 4L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 4L, 2L, 3L, 3L, 1L, 1L, 4L, 1L, 1L, 2L, 1L, 
1L, 4L, 4L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 4L, 4L, 4L, 1L, 4L, 1L, 
1L, 1L, 1L, 2L, 1L, 1L, 2L, 4L, 1L, 1L, 4L, 1L, 1L, 4L, 3L, 4L, 
4L, 1L, 2L, 1L, 4L, 1L, 3L, 1L, 2L, 2L, 2L), area = c(648L, 29L, 
2388L, 0L, 0L, 1247L, 0L, 2777L, 2777L, 7690L, 84L, 19L, 1L, 
143L, 0L, 31L, 23L, 113L, 0L, 47L, 600L, 8512L, 0L, 6L, 111L, 
274L, 678L, 28L, 474L, 9976L, 4L, 0L, 623L, 757L, 9561L, 1139L, 
2L, 342L, 0L, 51L, 115L, 9L, 128L, 43L, 22L, 0L, 49L, 284L, 1001L, 
21L, 28L, 1222L, 1L, 12L, 18L, 337L, 547L, 91L, 268L, 10L, 108L, 
249L, 0L, 132L, 0L, 0L, 109L, 246L, 36L, 215L, 28L, 112L, 1L, 
93L, 103L, 1904L, 1648L, 435L, 70L, 21L, 301L, 323L, 11L, 372L, 
98L, 181L, 583L, 0L, 236L, 10L, 30L, 111L, 0L, 3L, 587L, 118L, 
333L, 0L, 0L, 0L, 1031L, 1973L, 1L, 1566L, 0L, 447L, 783L, 0L, 
140L, 41L, 0L, 268L, 128L, 1267L, 925L, 121L, 195L, 324L, 212L, 
804L, 76L, 463L, 407L, 1285L, 300L, 313L, 9L, 11L, 237L, 26L, 
0L, 2150L, 196L, 72L, 1L, 30L, 637L, 1221L, 99L, 288L, 66L, 0L, 
0L, 0L, 2506L, 63L, 450L, 41L, 185L, 36L, 945L, 514L, 57L, 1L, 
5L, 164L, 781L, 0L, 84L, 236L, 245L, 178L, 0L, 9363L, 22402L, 
15L, 0L, 912L, 333L, 3L, 256L, 905L, 753L, 391L), population = c(16L, 
3L, 20L, 0L, 0L, 7L, 0L, 28L, 28L, 15L, 8L, 0L, 0L, 90L, 0L, 
10L, 0L, 3L, 0L, 1L, 1L, 119L, 0L, 0L, 9L, 7L, 35L, 4L, 8L, 24L, 
0L, 0L, 2L, 11L, 1008L, 28L, 0L, 2L, 0L, 2L, 10L, 1L, 15L, 5L, 
0L, 0L, 6L, 8L, 47L, 5L, 0L, 31L, 0L, 0L, 1L, 5L, 54L, 0L, 1L, 
1L, 17L, 61L, 0L, 10L, 0L, 0L, 8L, 6L, 1L, 1L, 6L, 4L, 5L, 11L, 
0L, 157L, 39L, 14L, 3L, 4L, 57L, 7L, 2L, 118L, 2L, 6L, 17L, 0L, 
3L, 3L, 1L, 1L, 0L, 0L, 9L, 6L, 13L, 0L, 0L, 0L, 2L, 77L, 0L, 
2L, 0L, 20L, 12L, 0L, 16L, 14L, 0L, 2L, 3L, 5L, 56L, 18L, 9L, 
4L, 1L, 84L, 2L, 3L, 3L, 14L, 48L, 36L, 3L, 0L, 22L, 5L, 0L, 
9L, 6L, 3L, 3L, 0L, 5L, 29L, 39L, 2L, 15L, 0L, 0L, 0L, 20L, 0L, 
8L, 6L, 10L, 18L, 18L, 49L, 2L, 0L, 1L, 7L, 45L, 0L, 1L, 13L, 
56L, 3L, 0L, 231L, 274L, 0L, 0L, 15L, 60L, 0L, 22L, 28L, 6L, 
8L), language = structure(c(10L, 6L, 8L, 1L, 6L, 10L, 1L, 2L, 
2L, 1L, 4L, 1L, 8L, 6L, 1L, 6L, 1L, 3L, 1L, 10L, 10L, 6L, 1L, 
10L, 5L, 3L, 10L, 10L, 3L, 1L, 6L, 1L, 10L, 2L, 7L, 2L, 3L, 10L, 
1L, 2L, 2L, 6L, 5L, 6L, 3L, 1L, 2L, 2L, 8L, 2L, 10L, 10L, 6L, 
1L, 1L, 9L, 3L, 3L, 10L, 1L, 4L, 4L, 1L, 6L, 1L, 1L, 2L, 3L, 
6L, 1L, 3L, 2L, 7L, 9L, 6L, 10L, 6L, 8L, 1L, 10L, 6L, 3L, 1L, 
9L, 8L, 10L, 10L, 1L, 10L, 8L, 10L, 10L, 4L, 4L, 10L, 10L, 10L, 
10L, 10L, 10L, 8L, 2L, 10L, 10L, 1L, 8L, 10L, 10L, 10L, 6L, 6L, 
1L, 2L, 3L, 10L, 10L, 8L, 6L, 8L, 6L, 2L, 1L, 2L, 2L, 10L, 5L, 
2L, 8L, 6L, 10L, 6L, 8L, 3L, 1L, 7L, 1L, 10L, 6L, 10L, 8L, 10L, 
1L, 1L, 1L, 8L, 6L, 6L, 4L, 8L, 7L, 10L, 10L, 3L, 10L, 1L, 8L, 
9L, 1L, 8L, 10L, 1L, 2L, 1L, 1L, 5L, 6L, 6L, 2L, 10L, 1L, 6L, 
10L, 10L, 10L), .Label = c("1", "2", "3", "4", "5", "6", "7", 
"8", "9", "10"), class = "factor"), bars = c(0L, 0L, 2L, 0L, 
3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 3L, 0L, 0L, 0L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 2L, 1L, 0L, 1L, 0L, 0L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
0L, 0L, 0L, 0L, 3L, 3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 3L, 
1L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 3L, 3L, 0L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 2L, 0L, 
0L, 3L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 
0L, 0L, 0L, 1L, 0L, 0L, 0L, 3L, 0L, 0L, 0L, 0L, 3L, 3L, 0L, 0L, 
3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 5L, 0L, 0L, 0L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
0L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 3L, 0L), stripes = c(3L, 0L, 
0L, 0L, 0L, 2L, 1L, 3L, 3L, 0L, 3L, 3L, 0L, 0L, 0L, 0L, 2L, 0L, 
0L, 0L, 5L, 0L, 0L, 0L, 3L, 2L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 2L, 
0L, 3L, 0L, 0L, 0L, 5L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 3L, 3L, 
3L, 3L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 5L, 3L, 3L, 1L, 9L, 0L, 0L, 
0L, 0L, 2L, 0L, 0L, 3L, 0L, 3L, 0L, 2L, 3L, 3L, 0L, 2L, 0L, 0L, 
0L, 0L, 3L, 0L, 5L, 0L, 3L, 2L, 0L, 11L, 2L, 3L, 2L, 3L, 14L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 5L, 3L, 0L, 3L, 1L, 0L, 3L, 
3L, 0L, 5L, 3L, 0L, 2L, 0L, 0L, 0L, 3L, 0L, 0L, 2L, 5L, 0L, 0L, 
0L, 3L, 0L, 0L, 3L, 2L, 0L, 0L, 3L, 0L, 3L, 0L, 0L, 0L, 0L, 3L, 
5L, 0L, 0L, 3L, 0L, 0L, 5L, 5L, 0L, 0L, 0L, 0L, 0L, 3L, 6L, 0L, 
9L, 0L, 13L, 0L, 0L, 0L, 3L, 0L, 0L, 3L, 0L, 0L, 7L), colours = c(5L, 
3L, 3L, 5L, 3L, 3L, 3L, 2L, 3L, 3L, 2L, 3L, 2L, 2L, 3L, 3L, 8L, 
2L, 6L, 4L, 3L, 4L, 6L, 4L, 5L, 3L, 3L, 3L, 3L, 2L, 5L, 6L, 5L, 
3L, 2L, 3L, 2L, 3L, 4L, 3L, 3L, 3L, 3L, 2L, 4L, 6L, 3L, 3L, 4L, 
2L, 4L, 3L, 3L, 6L, 7L, 2L, 3L, 3L, 3L, 4L, 3L, 3L, 3L, 2L, 3L, 
7L, 2L, 3L, 4L, 5L, 2L, 2L, 6L, 3L, 3L, 2L, 3L, 4L, 3L, 2L, 3L, 
3L, 3L, 2L, 4L, 2L, 4L, 4L, 3L, 4L, 4L, 3L, 3L, 3L, 3L, 3L, 4L, 
3L, 3L, 3L, 2L, 4L, 2L, 3L, 7L, 2L, 5L, 3L, 3L, 3L, 3L, 3L, 2L, 
3L, 2L, 3L, 4L, 3L, 3L, 2L, 3L, 4L, 6L, 2L, 4L, 2L, 3L, 2L, 7L, 
4L, 4L, 2L, 3L, 3L, 2L, 4L, 2L, 5L, 4L, 4L, 4L, 5L, 4L, 4L, 4L, 
4L, 2L, 2L, 4L, 3L, 4L, 3L, 4L, 2L, 3L, 2L, 2L, 6L, 4L, 5L, 3L, 
3L, 6L, 3L, 2L, 4L, 4L, 7L, 2L, 3L, 4L, 4L, 4L, 5L), red = c(1L, 
1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 
1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
0L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 
1L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 
1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 
0L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 
1L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 1L, 
1L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), green = c(1L, 
0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 
1L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 
0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 
0L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 
1L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 
1L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 
1L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
1L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 
1L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 
1L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 
0L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 1L), blue = c(0L, 
0L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 
0L, 1L, 0L, 1L, 1L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 
1L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 
1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 
1L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 
0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 
0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 0L, 
0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 
1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L), gold = c(1L, 
1L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 
0L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 
0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 
0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 1L, 
1L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
0L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 
0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 
1L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 
1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 
1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 1L), white = c(1L, 
0L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 
0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 1L, 
1L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 
1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 
1L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 
1L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
0L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 
1L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 1L), black = c(1L, 
1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 
0L, 1L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 
0L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 
0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 
0L, 0L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 1L), orange = c(0L, 
0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 
1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L), mainhue = structure(c(5L, 
7L, 5L, 2L, 4L, 7L, 8L, 2L, 2L, 2L, 7L, 2L, 7L, 5L, 2L, 4L, 2L, 
5L, 7L, 6L, 2L, 5L, 2L, 4L, 7L, 7L, 7L, 7L, 4L, 7L, 4L, 2L, 4L, 
7L, 7L, 4L, 5L, 7L, 2L, 2L, 2L, 8L, 8L, 7L, 2L, 5L, 2L, 4L, 1L, 
2L, 5L, 5L, 8L, 2L, 2L, 8L, 8L, 8L, 5L, 7L, 4L, 1L, 8L, 2L, 4L, 
2L, 2L, 4L, 4L, 5L, 1L, 2L, 2L, 7L, 2L, 7L, 7L, 7L, 8L, 8L, 8L, 
8L, 5L, 8L, 1L, 7L, 7L, 7L, 7L, 7L, 2L, 7L, 7L, 7L, 7L, 7L, 7L, 
7L, 7L, 2L, 5L, 5L, 2L, 7L, 2L, 7L, 4L, 2L, 3L, 7L, 8L, 2L, 2L, 
6L, 5L, 2L, 7L, 7L, 7L, 5L, 7L, 1L, 7L, 7L, 2L, 8L, 7L, 3L, 7L, 
7L, 5L, 5L, 5L, 5L, 8L, 5L, 2L, 6L, 8L, 7L, 4L, 5L, 2L, 5L, 7L, 
7L, 2L, 7L, 7L, 7L, 5L, 7L, 5L, 7L, 7L, 7L, 7L, 2L, 5L, 4L, 7L, 
8L, 8L, 8L, 7L, 7L, 4L, 7L, 7L, 7L, 7L, 5L, 5L, 5L), .Label = c("black", 
"blue", "brown", "gold", "green", "orange", "red", "white"), class = "factor"), 
    circles = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 1L, 0L, 0L, 1L, 0L, 1L, 4L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 
    0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 
    0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L), crosses = c(0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
    1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 2L, 1L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 
    0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 
    0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L), saltires = c(0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 
    0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L), quarters = c(0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
    0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 
    0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
    0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 
    0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 
    0L, 0L, 0L, 0L), sunstars = c(1L, 1L, 1L, 0L, 0L, 1L, 0L, 
    0L, 1L, 6L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 22L, 
    0L, 0L, 1L, 1L, 14L, 3L, 1L, 0L, 1L, 4L, 1L, 1L, 5L, 0L, 
    4L, 1L, 15L, 0L, 1L, 0L, 0L, 0L, 1L, 10L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 7L, 
    0L, 0L, 0L, 1L, 0L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 1L, 
    0L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    1L, 1L, 0L, 0L, 1L, 1L, 0L, 4L, 1L, 0L, 1L, 1L, 1L, 2L, 0L, 
    6L, 4L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 2L, 5L, 1L, 0L, 4L, 
    0L, 1L, 0L, 2L, 0L, 2L, 0L, 1L, 0L, 5L, 5L, 1L, 0L, 0L, 1L, 
    0L, 2L, 0L, 0L, 0L, 1L, 0L, 0L, 2L, 1L, 0L, 0L, 1L, 0L, 0L, 
    1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 50L, 1L, 0L, 0L, 7L, 1L, 
    5L, 1L, 0L, 0L, 1L), crescent = c(0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
    1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L), triangle = c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
    0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 
    0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 
    0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L), icon = c(1L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 1L, 0L, 
    1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 
    0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
    0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
    1L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 
    0L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
    1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L), animate = c(0L, 
    1L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 
    1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 
    0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 
    0L, 1L, 0L, 0L, 0L, 1L, 1L, 1L), text = c(0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
    0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
    0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
    0L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L), topleft = structure(c(1L, 6L, 4L, 2L, 
    2L, 6L, 7L, 2L, 2L, 7L, 6L, 2L, 7L, 4L, 2L, 1L, 6L, 4L, 7L, 
    5L, 2L, 4L, 7L, 7L, 7L, 6L, 2L, 7L, 4L, 6L, 6L, 7L, 2L, 2L, 
    6L, 3L, 4L, 6L, 7L, 2L, 2L, 7L, 7L, 6L, 7L, 4L, 2L, 3L, 6L, 
    2L, 4L, 4L, 7L, 7L, 7L, 7L, 2L, 2L, 4L, 6L, 1L, 1L, 7L, 2L, 
    6L, 6L, 2L, 6L, 6L, 1L, 1L, 2L, 7L, 6L, 2L, 6L, 4L, 6L, 4L, 
    2L, 4L, 6L, 3L, 7L, 1L, 6L, 1L, 6L, 6L, 6L, 4L, 2L, 2L, 6L, 
    7L, 1L, 2L, 6L, 7L, 2L, 4L, 4L, 2L, 6L, 7L, 6L, 4L, 2L, 2L, 
    6L, 7L, 7L, 2L, 5L, 4L, 2L, 6L, 6L, 6L, 7L, 7L, 6L, 6L, 6L, 
    2L, 7L, 6L, 7L, 2L, 6L, 4L, 4L, 4L, 4L, 6L, 2L, 2L, 5L, 7L, 
    6L, 3L, 4L, 2L, 2L, 6L, 4L, 2L, 6L, 6L, 2L, 4L, 6L, 6L, 7L, 
    7L, 6L, 6L, 7L, 6L, 1L, 7L, 7L, 7L, 2L, 6L, 1L, 3L, 3L, 6L, 
    2L, 2L, 4L, 4L, 4L), .Label = c("black", "blue", "gold", 
    "green", "orange", "red", "white"), class = "factor"), botright = structure(c(5L, 
    7L, 8L, 7L, 7L, 1L, 2L, 2L, 2L, 2L, 7L, 2L, 7L, 5L, 2L, 7L, 
    7L, 5L, 7L, 7L, 2L, 5L, 2L, 4L, 7L, 5L, 7L, 8L, 4L, 7L, 5L, 
    2L, 4L, 7L, 7L, 7L, 5L, 7L, 2L, 2L, 2L, 8L, 7L, 7L, 5L, 5L, 
    2L, 7L, 1L, 2L, 7L, 7L, 8L, 2L, 2L, 8L, 7L, 7L, 2L, 5L, 4L, 
    4L, 7L, 2L, 7L, 7L, 2L, 5L, 5L, 5L, 7L, 2L, 2L, 5L, 2L, 8L, 
    7L, 1L, 6L, 2L, 7L, 5L, 4L, 8L, 5L, 7L, 5L, 2L, 7L, 7L, 2L, 
    7L, 7L, 2L, 5L, 5L, 8L, 7L, 7L, 2L, 5L, 7L, 2L, 7L, 2L, 7L, 
    4L, 2L, 2L, 2L, 8L, 2L, 2L, 5L, 5L, 2L, 1L, 7L, 5L, 5L, 8L, 
    1L, 2L, 7L, 7L, 7L, 7L, 3L, 7L, 5L, 5L, 5L, 7L, 2L, 8L, 5L, 
    2L, 2L, 8L, 1L, 4L, 7L, 2L, 5L, 1L, 5L, 2L, 7L, 1L, 7L, 2L, 
    7L, 5L, 7L, 8L, 7L, 7L, 2L, 1L, 7L, 7L, 8L, 8L, 7L, 7L, 5L, 
    8L, 7L, 7L, 7L, 7L, 5L, 3L, 5L), .Label = c("black", "blue", 
    "brown", "gold", "green", "orange", "red", "white"), class = "factor")), .Names = c("ytrain", 
"landmass", "zone", "area", "population", "language", "bars", 
"stripes", "colours", "red", "green", "blue", "gold", "white", 
"black", "orange", "mainhue", "circles", "crosses", "saltires", 
"quarters", "sunstars", "crescent", "triangle", "icon", "animate", 
"text", "topleft", "botright"), row.names = c(NA, -174L), class = "data.frame")
tdata$language <-  factor(tdata$language)
tdata$ytrain <- factor(tdata$ytrain)

library("coin")

m <- ctree(ytrain ~ language, data = subset(tdata, language != "8"), 
    control = ctree_control(testtype = "Univariate", maxdepth = 1L))
it <- independence_test(ytrain ~ language, data = subset(tdata, language != "8"), 
                        teststat = "quad")
stopifnot(isTRUE(all.equal(m@tree$criterion$statistic, 
                           statistic(it), check.attributes = FALSE)))

### easier example
levels(tdata$language) <- c(1, 1, 1, 1, 1, 1, 2, 8, 1, 1)
levels(tdata$ytrain) <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 6)
m <- ctree(ytrain ~ language, data = subset(tdata, language != "8"),
    control = ctree_control(testtype = "Univariate", maxdepth = 1L))
it <- independence_test(ytrain ~ language, data = subset(tdata, language != "8"), 
                        teststat = "quad")
stopifnot(isTRUE(all.equal(m@tree$criterion$statistic, 
                           statistic(it), check.attributes = FALSE)))

## the whole exercise manually
Y <- model.matrix(~ language - 1, data = subset(tdata, language != "8"))
X <- model.matrix(~ ytrain -1, data = subset(tdata, language != "8"))
w <- rep(1, nrow(X))

lin <- coin:::LinearStatistic(Y, X, weights = w)
expcov <- coin:::ExpectCovarLinearStatistic(Y, X, weights = w)

tmp <- new("LinStatExpectCovar", ncol(Y), ncol(X))
tmp@linearstatistic <- lin
tmp@expectation <- expcov@expectation
tmp@covariance <- expcov@covariance

a <- .Call("R_linexpcovReduce", tmp)

u <- matrix(tmp@linearstatistic - tmp@expectation, nc = 1)
d <- tmp@dimension
u <- matrix(tmp@linearstatistic - tmp@expectation, nc = 1)[1:d,,drop = FALSE]
S <- coin:::MPinv(matrix(as.vector(tmp@covariance[1:d^2]), ncol = d))

stat <- t(u) %*% S$MPinv %*% u
stopifnot(isTRUE(all.equal(stat[1,1], statistic(it), 
                           check.attributes = FALSE)))

x <- matrix(as.vector(tmp@covariance[1:d^2]), ncol = d)
s <- svd(x)

m <- new("svd_mem", 18L)
m@p <- as.integer(d)

s2 <- .Call("R_svd", x, m)

stopifnot(max(abs(s$d - m@s[1:d])) < sqrt(.Machine$double.eps))
stopifnot(max(abs(s$v - t(matrix(m@v[1:d^2], nrow = d)))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(s$u - matrix(m@u[1:d^2], nrow = d))) < sqrt(.Machine$double.eps))

s2 <- .Call("R_svd", tmp@covariance, m)

stopifnot(max(abs(s$d - m@s[1:d])) < sqrt(.Machine$double.eps))
stopifnot(max(abs(s$v - t(matrix(m@v[1:d^2], nrow = d)))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(s$u - matrix(m@u[1:d^2], nrow = d))) < sqrt(.Machine$double.eps))

a <- .Call("R_MPinv", tmp@covariance, sqrt(.Machine$double.eps), m)

stat <- t(u) %*% matrix(a@MPinv[1:d^2], ncol = d) %*% u  
stopifnot(isTRUE(all.equal(stat[1,1], statistic(it), 
                           check.attributes = FALSE)))
