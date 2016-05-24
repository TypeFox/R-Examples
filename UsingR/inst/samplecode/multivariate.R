
## @knitr 
head(state.x77, n=3)


## @knitr 
state.x77[1,2]               # first row, second column
state.x77[1, 3:4]            # first row, third and fourth columns


## @knitr 
state.x77[state.x77[,"Population"] < 500, 1:6]


## @knitr 
colnames(state.x77)


## @knitr 
rbind(1:5, c(1, 1, 2, 3, 5))            # 2 by 5 matrix without names


## @knitr 
m <- cbind(1:3, c(1.1, 1.2, 1.3), c(1, 1, 2)) # a 3 by 3 matrix 
colnames(m) <- c("x", "y", "z")               # or cbind(x=..., ...)
m


## @knitr 
m[1:6]                                  # first 6 entries of m


## @knitr 
matrix(1:10, nrow=2)
matrix(1:10, nrow=2, byrow=TRUE)


## @knitr 
DF <- Cars93[1:3, 1:5]
DF


## @knitr 
DF[ , "Price"]


## @knitr 
DF[ , "Price", drop=FALSE]


## @knitr 
DF["Price"]


## @knitr 
require(MASS)
dim(UScereal)                           # rows, columns
length(levels(UScereal$mfr))            # number of manufacturers
length(levels(UScereal["vitamins"]))    # vitamin categories
sum(UScereal[, "sugars"] > 10)          # sugar levels above 10
mean(UScereal[UScereal$fat > 5, "calories"]) # conditional mean
mean(UScereal[UScereal$fat <= 5, "calories"]) # conditional mean
mean(UScereal[UScereal["shelf"] == 2, "calories"])


## @knitr 
l <- lm(mpg ~ wt, data=mtcars)
length(l)                               # no. components
names(l)
l[["residuals"]]                        # numeric, named data


## @knitr 
d <- data.frame(a=1, b="two")
names(d) <- c("A", "B")
dimnames(d) <- list("1", c("eh", "bee"))
colnames(d) <- c("EH", "BEE")
d <- setNames(d, c("ahh", "buh"))


## @knitr eval=FALSE
## data.frame(nm1 = vec1, nm2=vec2, ...)


## @knitr 
d <- as.data.frame(state.x77)
class(d)


## @knitr 
d1 <- data.matrix(d)
class(d1)


## @knitr 
d <- mtcars[1:5,]                       # first 5 rows
mean(d$mpg) - sd(d$mpg) <= d$mpg & d$mpg <= mean(d$mpg) + sd(d$mpg)


## @knitr 
with(d, mean(mpg) - sd(mpg) <= mpg & mpg <= mean(mpg) + sd(mpg))


## @knitr 
d <- Cars93[1:3, 1:4]                   # first 3 rows, 4 columns
d[1,1] <- d[3,4] <- NA                  # set two values to NA
d


## @knitr 
d[1:2, 4] <- round(d[1:2, 4])           # round numeric values


## @knitr 
d[3,c(2,4)] <- list("A3", 30)           # warning


## @knitr 
levels(d$Model) <- c(levels(d$Model), c("A3", "A4", "A6"))
d[3,c(2,4)] <- list("A3", 30)    


## @knitr 
d[4, ] <- list("Audi", "A4", "Midsize", 35)


## @knitr 
d <- rbind(d, list("Audi", "A6", "Large", 45))


## @knitr 
d[, 5] <- d$Min.Price * 1.3             # price in Euros


## @knitr 
d$Min.Price.Euro <- d$Min.Price * 1.3


## @knitr 
names(d) <- tolower(names(d))


## @knitr 
names(d)[3] <- "car type"


## @knitr 
aq <- airquality[1:5, ]                 # shorten
aq
subset(aq, select = Ozone:Wind)         # range of names
subset(aq, select = -c(Month, Day))     # same result
subset(aq, subset = !is.na(Ozone), select=Ozone:Wind) # drop a row


## @knitr 
DF <- data.frame(a=c(NA, 1, 2), b=c("one", NA, "three"))
subset(DF, !is.na(a))                   # drop first, keep second
subset(DF, complete.cases(DF))           # drop first, second


## @knitr 
d$min.price.euro <- d$min.price * 1.5


## @knitr 
d <- within(d, {min.price.euro = min.price * 1.5})


## @knitr 
d <- transform(d, min.price.euro = min.price * 1.5)


## @knitr echo=FALSE
speed <- head(reshape(morley, v.names="Speed", timevar="Expt", idvar="Run", direction="wide"))
speed <- speed[,-1]
rownames(speed) <- 1:6


## @knitr 
speed


## @knitr 
m <- reshape(speed, varying=names(speed)[1:5], direction="long")
head(m)                                 # first 6 rows only


## @knitr 
speed$Run <- LETTERS[1:6]
m <- reshape(speed, varying=names(speed)[1:5], direction="long")
head(m)


## @knitr 
reshape(m, v.names="Speed", timevar="time", idvar="Run", 
        direction="wide")


## @knitr cache=FALSE
domestic <- "
The Avengers,	   	623357910
The Dark Knight Rises,	448139099
The Hunger Games,     	408010692
Skyfall,	   	304360277
The Hobbit,	     	303003568
"
foreign <- "
The Avengers,           1511.8
Skyfall,		1108.6
The Dark Knight Rises,	1084.4
The Hobbit,        	1017.0
Ice Age,  	         877.2
"


## @knitr cache=FALSE
df.domestic <- read.csv(textConnection(domestic), header=FALSE)
names(df.domestic) <- c("Name", "Domestic")
df.foreign <- read.csv(textConnection(foreign), header=FALSE)
names(df.foreign) <- c("name", "foreign")


## @knitr cache=FALSE
merge(df.domestic, df.foreign, by.x="Name", by.y="name", all=FALSE)


## @knitr cache=FALSE
merge(df.domestic, df.foreign, by.x="Name", by.y="name", all=TRUE)


## @knitr cache=FALSE
merge(df.domestic, df.foreign, by.x="Name", by.y="name", all.x=TRUE)


## @knitr 
babies$id <- as.character(babies$id)    # change type of variable


## @knitr 
with(babies, gestation[gestation > 45 * 7])


## @knitr 
babies <- within(babies, {
  gestation[gestation == 999] <- NA
  wt[wt == 999] <- NA
  wt1[wt1 == 999] <- NA
  dwt[dwt == 999] <- NA
  ht[ht == 99] <- NA
  dht[dht == 99] <- NA
})


## @knitr 
babies$smoke <- factor(babies$smoke)
levels(babies$smoke) <- list("never"=0, "smokes now"=1, 
                             "Until current pregnancy"=2, 
                             "once did, not now"=3)


## @knitr 
babies$number <- factor(babies$number)
levels(babies$number) <- list("never"=0, "1-4"=1, "5-9"=2,
                              "10-14"=3, "15-19"=4, "20-29"=5,
                              "30-39"=6, "40-60"=7, "60+"=8,
                              "smoke, don't know"=9, "unknown"=98)


## @knitr 
require(lubridate)
(x <- ymd("1961-01-01"))


## @knitr 
x + 10 * days(1)


## @knitr 
babies$date <- x + (babies$date - 1096) * days(1)


## @knitr 
bmi <- function(wt, ht) (wt/2.2) / (ht*2.54/100)^2
babies <- transform(babies, bmi = bmi(wt1, ht),
                    dbmi = bmi(dwt, dht))


## @knitr 
subset(babies, abs(dbmi - bmi) > 14, 
       select=c(date, gestation, wt, race))


## @knitr 
mean(fat$neck) / mean(fat$wrist) 
mean(fat$neck/fat$wrist)  


## @knitr 
with(fat, mean(neck) / mean(wrist))
with(fat, mean(neck / wrist))


## @knitr 
subset(Cars93, Origin == "non-USA" & Cylinders == 4 & Max.Price <= 15)


## @knitr 
x <- 1:5
x - 3


## @knitr 
x <- 1:5
res <- integer(length(x))               # allocate temporary storage
for(i in 1:length(x)) {                 # iterate over indices
  res[i] <- x[i] - 3                    # compute f, do bookkeeping
}
res  


## @knitr 
vmedian <- Vectorize(median)            # returns a function
vmedian(homedata)


## @knitr 
collection <- c(4, 9, 16)
Map(sqrt, collection)


## @knitr 
sqrt(collection)


## @knitr 
sapply(collection, sqrt)


## @knitr 
lst <- with(ToothGrowth, split(len, supp)) 
sapply(lst, mean)


## @knitr 
with(ToothGrowth, 
     tapply(len, supp, mean)            # (X, INDEX, FUN)
     )


## @knitr 
with(ToothGrowth, 
     tapply(len, list(supp, dose), mean)            # (X, INDEX, FUN)
     )


## @knitr 
aggregate(len ~ supp, data=ToothGrowth, mean)


## @knitr 
mean_formula <- function(formula, data) {
  out <- aggregate(formula, data=data, mean)
  xtabs(formula, out)
}


## @knitr 
mean_formula(len ~ supp, data=ToothGrowth)


## @knitr echo=FALSE
mean.formula <- function(...) mean_formula(...)


## @knitr 
mean(len ~ supp, data=ToothGrowth)


## @knitr 
lst <- with(ToothGrowth, split(len, supp))
sapply(lst, summary)


## @knitr 
sapply(mtcars, mean)


## @knitr 
m <- mtcars[1:4, 1:3]
m[1,1] <- m[2,2] <- NA
sapply(m, mean)
sapply(m, mean, na.rm=TRUE)


## @knitr 
m <- rbind(c(1,2), c(3,4))
sqrt(m)


## @knitr 
(m <- replicate(5, rnorm(3)))


## @knitr 
rowSums(m)                              # add along rows
colSums(m)                              # add down columns


## @knitr 
sapply(m, sum)


## @knitr 
apply(m, 1, mean)                       # rowMeans alternative
apply(m, 2, mean)                       # colMeans alternative


## @knitr 
c(sum(m[1,]), sum(m[2,]), sum(m[3,]))


## @knitr 
apply(m, 2, summary)


## @knitr 
xbars <- apply(m, 2, mean)
centers <- sweep(m, 2, xbars, FUN="-")  # "-" is default
centers


## @knitr 
sds <- apply(m, 2, sd)
z_scores <- sweep(centers, 2, sds, FUN="/")
z_scores


## @knitr 
min(3, 4)                              # 3
min(c(1,4), c(2,3))                    # not c(1,3) as maybe desired


## @knitr 
Map(min, c(1,4), c(2,3))


## @knitr 
mapply(min, c(1,4), c(2,3))


## @knitr 
our_sweep <- function(col, center) col - center
mapply(our_sweep, as.data.frame(m), apply(m, 2, mean))


## @knitr 
body <- Animals$body; brain <- Animals$brain 
do.call(cor, Map(rank, list(body, brain)))


## @knitr 
do.call(cor, Map(rank, setNames(Animals, NULL)))


## @knitr 
m <- Cars93[1:2, 1:15]                  # 15 columns
Filter(is.factor, m)                    # 6 are factors


## @knitr 
Reduce("+", 1:4)


## @knitr 
Reduce(function(x,y) ifelse(x > y, x, y), c(6, 4, 7, 9, 10, 3))


## @knitr 
## gcd by Euclidean algorithm, a, b integers
gcd <- function(a, b) {while (b != 0) {t = b; b = a %% b; a = t}; a}
## scm (lcm is R function name in graphics package)
scm <- function(a, b) (a * b) / gcd(a, b)


## @knitr 
scm(3, 5)                               # no common primes
scm(3, 6)                               # common prime


## @knitr 
Reduce(scm, 1:20)      # smallest number divisible by 1, 2, ..., 20


## @knitr eval=FALSE
## sapply(wellbeing[,-(1:2)], function(y) {
##        cor(wellbeing[,2], y, use="complete.obs")
##      })


## @knitr 
library(LearnEDA)
l <- with(beatles, split(time,  album))
sapply(l, length)


## @knitr 
sapply(mtcars, sd)
Vectorize(sd)(mtcars)


## @knitr 
sapply(Filter(is.numeric, Cars93), sd)


## @knitr 
sapply(Filter(is.numeric, Cars93), sd, na.rm=TRUE)


## @knitr 
teams <- split(batting, batting$teamID)
team_avg <- sapply(teams, function(DF) with(DF, sum(H) / sum(AB)))
sort(team_avg)


## @knitr 
players <- split(batting, batting$playerID)
traded_players <- Filter(function(x) nrow(x) > 1, players)
names(traded_players)


## @knitr 
d <- c("1,2","3,4","5,6")
strsplit(d,",")


## @knitr 
sapply(d, function(x) x[1])


## @knitr 
d <- data.frame(a=1:3, b=c(1, NA, 3), c=c("one", "two", NA))
Filter(function(x) !any(is.na(x)), d)


## @knitr 
f <- function(nm, x) sprintf("Variable %s has class %s", nm, class(x)[1])
our_func <- function(DF) mapply(f, names(DF), DF)
our_func(mtcars[1:3])



## @knitr 
fruits <- c("Bananas", "Oranges", "Avocados", "Celeries?")
sapply(fruits, function(x) 
   paste(x, "are fruit number", which(fruits==x)))


## @knitr 
sapply(seq_along(fruits), function(i) paste(fruits[i], "are fruit number", i))


## @knitr 
mapply(paste, fruits, "are fruit number", seq_along(fruits))


## @knitr 
require("gdata")                        # must be installed
f <- "http://www.eia.gov/petroleum/gasdiesel/xls/pswrgvwall.xls"
gas_prices <- read.xls(f, sheet=2, skip=2)
gas_prices <- setNames(gas_prices[,1:2], c("Date", "Weekly_US"))


## @knitr gas_price_graph, eval=FALSE
## gas_prices$Date <- as.Date(substr(gas_prices$Date, 1, 10),
##                            format="%b %d%Y")
## plot(Weekly_US ~ Date, gas_prices, type="l")


## @knitr echo=FALSE, out.width=singlewide
gas_prices$Date <- as.Date(substr(gas_prices$Date, 1, 10), 
                           format="%b %d%Y")
plot(Weekly_US ~ Date, gas_prices, type="l")


## @knitr 
key <- "0AoaQTPQhRgkqdEthU0ZZeThtcWtvcWpZUThiX2JUMGc"
f <- paste("https://docs.google.com/spreadsheet/pub?key=", 
           key,
           "&single=true&gid=0&output=csv", sep="")
require(RCurl)
read.csv(textConnection(getURL(f)), header=TRUE)


## @knitr quandl, eval=FALSE
## require(Quandl)
## ch_0014 <- Quandl("WORLDBANK/CHN_SP_POP_0014_TO_ZS")
## ch_1564 <- Quandl("WORLDBANK/CHN_SP_POP_1564_TO_ZS")
## ch_65up <- Quandl("WORLDBANK/CHN_SP_POP_65UP_TO_ZS")
## ch_all <- Reduce(function(x,y) merge(x, y, by="Date"),
##                  list(ch_0014, ch_1564, ch_65up))
## names(ch_all) <- c("Date", "[0,14]", "[15,64]", "[65,)")


## @knitr echo=FALSE
ch_0014 <- suppressWarnings(Quandl("WORLDBANK/CHN_SP_POP_0014_TO_ZS"))
ch_1564 <- suppressWarnings(Quandl("WORLDBANK/CHN_SP_POP_1564_TO_ZS"))
ch_65up <- suppressWarnings(Quandl("WORLDBANK/CHN_SP_POP_65UP_TO_ZS"))
ch_all <- Reduce(function(x,y) merge(x, y, by="Date"), list(ch_0014, ch_1564, ch_65up))
names(ch_all) <- c("Date", "[0,14]", "[15,64]", "[65,)")


## @knitr chinese_demographics, eval=FALSE
## heights <- t(ch_all[,-1])
## colnames(heights) <- format(ch_all[,'Date'], format="%Y")
## barplot(heights, main="Proportion of [0-14], [15-64], [65,)")


## @knitr 
require(RJSONIO)
f <- "http://www.quandl.com/api/v1/datasets/PRAGUESE/PX.json"
out <- fromJSON(f)
out$column_names                        # names
out$data[1]                             # one from 1000s of values


## @knitr portugese_stocks, eval=FALSE
## pluck <- function(l, key) l[[key]]      # pluck from a list
## px <- data.frame(Date = as.Date(sapply(out$data, pluck, key=1)),
##                  index = sapply(out$data, pluck, key=2),
##                  perc_change =  sapply(out$data, pluck, key=3))
## plot(index ~ Date, data=px, type="l", main="Portugese stock index")


## @knitr echo=FALSE, out.width=doublewide
heights <- t(ch_all[,-1])
colnames(heights) <- format(ch_all[,'Date'], format="%Y")
barplot(heights, main="Proportion of [0-14], [15-64], [65,)")
pluck <- function(l, key) l[[key]]      # pluck from a list
px <- data.frame(Date = as.Date(sapply(out$data, pluck, key=1)),
                 index = sapply(out$data, pluck, key=2),
                 perc_change =  sapply(out$data, pluck, key=3))
plot(index ~ Date, data=px, type="l", main="Portugese stock index")


## @knitr 
require(XML)
## fit in 80 characters
url_base = "http://en.wikipedia.org/wiki/"
ch <- "List_of_highest-grossing_films_in_China"
us_can <- "List_of_highest-grossing_films_in_Canada_and_the_United_States"
##
china_all  <- readHTMLTable(paste(url_base, ch, sep=""))[[1]]
us_can_all <- readHTMLTable(paste(url_base, us_can, sep=""))[[2]]


## @knitr cache=FALSE
in_common <- merge(china_all, us_can_all, by="Title")
## tidy up
elide <- function(x, n=20) 
  ifelse(nchar(x) < n, x, sprintf("%s...", substr(x,0,n)))
rownames(in_common) <- sapply(as.character(in_common[,1]), elide)
##
in_common[, c(2,6,7)]


