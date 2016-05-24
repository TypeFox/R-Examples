
set.seed(290875)
library("party")
library("coin")

"hohnloser" <-
structure(list(EF = as.integer(c(11, 11, 12, 13, 13, 13, 15,
17, 20, 20, 20, 20, 20, 21, 22, 22, 22, 22, 23, 24, 24, 24, 24,
24, 24, 24, 25, 25, 26, 26, 26, 27, 28, 30, 30, 31, 31, 32, 33,
33, 33, 33, 34, 34, 34, 34, 36, 37, 38, 38, 38, 39, 40, 41, 41,
41, 43, 43, 43, 44, 44, 49, 50, 51, 51, 51, 52, 52, 52, 56, 56,
56, 57, 57, 58, 58, 58, 59, 60, 60, 61, 64, 64, 64, 64, 65, 70,
70, 72, 75, 77, 77, 80, 93)), month = as.integer(c(1, 5, 14,
2, 10, 39, 16, 17, 1, 1, 1, 8, 29, 22, 1, 3, 11, 15, 13, 1, 1,
3, 5, 7, 11, 33, 3, 16, 1, 13, 23, 20, 12, 1, 1, 18, 20, 23,
9, 12, 17, 21, 1, 5, 14, 38, 6, 1, 3, 12, 18, 8, 19, 3, 10, 15,
19, 31, 33, 23, 24, 5, 13, 4, 21, 28, 3, 16, 37, 1, 3, 33, 23,
29, 5, 9, 36, 19, 1, 10, 7, 1, 6, 7, 14, 6, 5, 23, 36, 30, 10,
20, 7, 22)), cens = as.integer(c(0, 1, 0, 1, 0, 0, 1, 0, 1, 1,
1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0,
0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1,
0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
))), .Names = c("EF", "month", "cens"), class = "data.frame", row.names =
c("1",
"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
"14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
"25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35",
"36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46",
"47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57",
"58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68",
"69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79",
"80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90",
"91", "92", "93", "94"))


### get rid of the NAMESPACE
attach(asNamespace("party"))

### 
###
###    Regression tests for cutpoint search
###    
###    functions defined in file `./src/Splits.c'    

### tests for function C_Split
x <- rnorm(100)
y <- rnorm(100)
weights <- rep(1, length(x))
splitctrl <- new("SplitControl")
split <- Split(x, y, weights, splitctrl)
mydata <- data.frame(y, x)
ms <- show(maxstat_test(y ~ x, data = mydata, distribution = approximate(10)))
stopifnot(isequal(split[[1]], ms$estimate[[1]]))
stopifnot(isequal(split[[2]], ms$statistic))
stopifnot(isequal(max(split[[3]]), ms$statistic))

### Hohnloser data
ms <-  show(maxstat_test(Surv(month, cens) ~ EF, data = hohnloser,
distribution = approximate(10)))
splitctrl <- new("SplitControl")
splitctrl@minprob <- 0.1
splitctrl@minsplit <- as.integer(5)

split <- Split(hohnloser$EF, logrank_trafo(Surv(hohnloser$month, hohnloser$cens)),
               rep(1, nrow(hohnloser)), splitctrl)
stopifnot(isequal(split[[1]], ms$estimate[[1]]))
stopifnot(isequal(split[[2]], ms$statistic))
stopifnot(isequal(max(split[[3]]), ms$statistic))

### categorical splits
n <- 100
xf <- gl(5, 100/5)
yf <- gl(4, 100/4)[sample(1:length(xf))]
weights <- rep(1, length(xf))
splitctrl <- new("SplitControl")
splitctrl@minprob <- 0.1
splitctrl@minsplit <- as.integer(5)
split <- Split(xf, yf, weights, splitctrl)
split

### Check if the statistic used for selecting the split is
### correct: For the ranks of a continuous response the statistic
### needs to be equal to the standardized Wilcoxon statistic

y <- rnorm(100) + c(rep(0, 25), rep(1, 25), rep(0, 25), rep(1, 25))
x <- gl(4, 25)
weights <- rep(1, length(y))
split <- Split(x, rank(y), weights, splitctrl)
levelset <- levels(x)[split[[4]] == 1]
tstat <- split[[2]]
p <- wilcox.test(y ~ I(x %in% levelset),corr = FALSE,
                alternative = "less")$p.value
stopifnot(isequal(round(abs(qnorm(p)), 6), round(tstat, 6)))

y <- rnorm(100) + c(rep(0, 25), rep(1, 25), rep(0, 25), rep(1, 25))
x <- rnorm(100)
weights <- rep(1, length(y))
split <- Split(x, rank(y), weights, splitctrl)
tstat <- split[[2]]
p <- wilcox.test(y ~ I(x <= split[[1]]), corr = FALSE,
                alternative = "less")$p.value
stopifnot(isequal(round(abs(qnorm(p)), 6), round(tstat, 6)))
