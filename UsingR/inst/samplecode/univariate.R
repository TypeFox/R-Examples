
## @knitr whales
whale <- c(74, 122, 235, 111, 292, 111, 211, 133, 156, 79)


## @knitr 
length(whale)


## @knitr 
sum(whale)


## @knitr 
sum(whale)/length(whale)


## @knitr 
mean(whale)


## @knitr 
whale - mean(whale)                     # mean(whale) recycled
whale^2 / length(whale)
sqrt(whale)


## @knitr hip
hip_cost <- c(10500, 45000, 74100, NA, 83500, 86000, 38200, NA, 
              44300, 12500, 55700, 43900, 71900, NA, 62000)


## @knitr 
sum(hip_cost)


## @knitr 
sum(hip_cost, na.rm = TRUE)


## @knitr 
mean(hip_cost, na.rm=TRUE)


## @knitr 
head(precip)


## @knitr rainiest
head( sort(precip, decreasing=TRUE) )


## @knitr 
head(names(precip))


## @knitr 
test_scores <- c(Alice = 87, Bob = 72, Shirley = 99)


## @knitr 
test_scores <- setNames(c(87, 82, 99), c("Alice", "Bob", "Shirley"))


## @knitr test_scores
test_scores <- c(87, 782, 99)
names(test_scores) <- c("Alice", "Bob", "Shirley")
test_scores


## @knitr 
x <- c(1, "two", "III")
x


## @knitr 
as.numeric("1")
as.character(1)


## @knitr 
whale <- scan("whale.txt")              # or scan(file.choose())


## @knitr colon
1:5                                     # 1 2 3 4 5
1:length(whale)                         # 1 2 3 ... 10
0:(length(whale) - 1)                   # 0 1 2 ... 9


## @knitr 
seq(0, 100, by=10)                      # count by 10s


## @knitr 
seq(0, 100, length.out = 11)            # counts by 10s as well


## @knitr 
rep(5, times=10)                        # 10 5s


## @knitr 
rep(1:3, times=4)


## @knitr 
rep(c(1,2,3), times=c(3,2,1))           # or rep(1:3, 3:1)


## @knitr indexing
whale                                   # all the values
whale[1]                                # first element
whale[2]                                # second element
whale[11]                               # 11th element is not there


## @knitr 
whale[c(1,3,5,7,9)]


## @knitr all_but
whale[-1]


## @knitr 
precip[c("Seattle Tacoma", "New York")] # which is rainier?


## @knitr 
match(c("Seattle Tacoma", "New York"), names(precip))


## @knitr 
precip["Seattle"]                   # needs "Seattle Tacoma" to match


## @knitr 
x <- c(1, 2, 3)
x[1] <- 11                              # 11 is now first value
x


## @knitr 
x[2:3] <- c(12, 13)
x


## @knitr 
x[6] <- 6
x


## @knitr index-table, results="asis", echo=FALSE
tbl <- "
\\code{x[1]}: The first element of \\code{x}.
\\code{x[]}:  All elements of \\code{x}.
\\code{x[length(x)]}: The last element of \\code{x}.
\\code{x[c(2,3)]}: The second and third elements of \\code{x}.
\\code{x[-c(2,3)]}: All but the second and third elements of \\code{x}.
\\code{x[0]}: $0$-length vector of same type as \\code{x}.
\\code{x[1] <- 5}: Assign a value of 5 to first element of \\code{x}.
\\code{x[c(2,3)] <- c(4,5)}: Assign values to second and third elements of \\code{x}. In assignment, recycling of the right hand side may occur. Assignment can grow the length of a data vector.
"
out <- read.dcf(textConnection(tbl))
d <- data.frame(Command=colnames(out), Description=out[1,], stringsAsFactors=FALSE)
cat(booktabs(d, 
         caption="Various uses of indexing", 
         subtop=NULL, 
         rnames=NULL, 
         label="tab:indexing-options",
         colTypes=c("l", "p{3.2in}")
         ))


## @knitr recycling
x[2:3] <- 0
x


## @knitr 
x <- 1:10
x[] <- 1:3                              # 10 on left, 3 on right
x


## @knitr 
n <- 10                                 # get 10 elements of x
x[1 + 0:(n-1) %% length(x)]             # use remainder for indices


## @knitr 
c(class(1), class(pi), class(seq(1, 5, by=1)))


## @knitr 
sqrt(2) * sqrt(2)                       # looks right
sqrt(2) * sqrt(2) - 2                   # a difference


## @knitr 
c("Lincoln", "said", "\"Four score", 'and seven years ago..."')


## @knitr 
sprintf("X%s", 1:10)


## @knitr 
sprintf("%8d", c(1, 12, 123, 1234, 12345))


## @knitr 
paste("X", 1:10, sep="")


## @knitr 
paste("The", "quick", "brown", "fox", "...", sep="_")


## @knitr 
paste(c("Four","The"), c("Score","quick"), c("and","fox"), sep="_")


## @knitr 
x <- c(one=1, two=2, three=3)
out <- paste(names(x), x, sep=", ", collapse=": ")
sprintf("[ %s ]", out)


## @knitr 
x <- paste("X", rep(1:3, 4), sep="")
y <- factor(x)
y


## @knitr 
y[1] <- "X4"
y


## @knitr 
levels(y) <- c(levels(y), "X4")         # add "X4" to existing levels
y[1] <- "X4"
y


## @knitr 
levels(y) <- paste("label", 1:4, sep="")
y


## @knitr 
y <- factor(state.name[1:5])
y
levels(y) <- c("South", "West", "West", "South", "West")
y


## @knitr echo=FALSE
set.seed(1121)
(x <- rnorm(20))


## @knitr 
y <- factor(state.name)[1:5]
y                                       # 50 levels
factor(y, levels=y)                     # levels are actual values


## @knitr 
r <- "red"; b <- "blue"; g <- "green"
factor(c(r,r,r,r,r,g,g,g,g,g,b,b,b,b,b))


## @knitr 
gl(3, 5, labels=c("red", "green", "blue"))


## @knitr 
m <- head(Cars93)
out <- m$Origin : m$AirBags
out
levels(out)


## @knitr 
require(lubridate)
current_time <- now()                   # current time
class(current_time)                     # cf. ?POSIXct 
as.numeric(current_time)                # seconds since Jan 1, 1970
month(current_time, label=TRUE)         # what month?


## @knitr parse_data_time
x <- "15-Feb-2013 07:57:34"
y <- parse_date_time(x, "dbYHMS")
year(y)                                 # can get pieces


## @knitr 
now() - days(1)
now() - hours(24)                       # same thing


## @knitr 
is.na(1)
is.numeric("one")
is.logical("false")                     # "false" is character


## @knitr 
3 < pi
"one" == 1                              # 1 is first coerced to "1" 
sqrt(2) * sqrt(2) == 2                  # floating point gotcha


## @knitr 
isTRUE(all.equal(sqrt(2) * sqrt(2), 2))


## @knitr 
whale <- c(74, 122, 235, 111, 292, 111, 211, 133, 156, 79)
whale > 100
whale == 111


## @knitr 
whale < 100 | whale > 200
whale > 100 & whale < 200


## @knitr 
any(whale > 300)


## @knitr 
all(whale > 50)


## @knitr 
which(whale < 100 | whale > 200)


## @knitr 
292 %in% whale                         # is 292 in whale?


## @knitr 
any(292 == whale)


## @knitr 
match(c(292, 293), whale)


## @knitr 
sum(whale > 200)


## @knitr 
whale <- c(74, 122, 235, 111, 292, 111, 211, 133, 156, 79)
whale[ whale > mean(whale) ]


## @knitr 
whale[whale < mean(whale)-sd(whale) | whale > mean(whale)+sd(whale)]


## @knitr 
hip_cost <- c(10500, 45000, 74100, NA, 83500, 86000, 38200, NA, 
              44300, 12500, 55700, 43900, 71900, NA, 62000)
hip_cost[ !is.na(hip_cost) ]            # store as variable ...


## @knitr 
x <- babies$dwt                    # Dad's weight, 999 is used for NA
x[ x == 999 ] <- NA                       
range(x, na.rm=TRUE)                    # need to be careful of NA


## @knitr 
age <- kid.weights$age
ht <- kid.weights$height
ht[ age >= 48 & age < 60]


## @knitr 
p <- c(2, 3, 5, 7, 11, 13, 17, 19) 


## @knitr 
x <- c(2, 5, 4, 10, 8)
x^2
x - 6
(x - 9)^2      


## @knitr 
rep("a", 10)
seq(1, 99, by=2)
rep(1:3, rep(3,3))
rep(1:3, 3:1)      
c(1:5, 4:1)


## @knitr 
primes_under_20 <- c(1, 2, 3, 5, 8, 13, 21, 34)
ns <- 1:10
recips <- 1/ns
cubes <- (1:6)^3
years <- 1964:2014
subway <- c(14, 18, 23, 28, 34, 42, 50, 59, 66, 72, 79, 86, 96, 103, 110)
by25 <- seq(0,1000, by=25)


## @knitr 
sum(abs(rivers - mean(rivers))) / length(rivers)


## @knitr 
-1:3                                    # like (-1):3


## @knitr 
1:2*3                                   # not like 1:(2*3)


## @knitr 
precip["Juneau"]


## @knitr 
j_cities <- c("Jackson", "Jacksonville", "Juneau")
precip[j_cities]


## @knitr 
precip[grepl("^J", names(precip))]


## @knitr 
paste("Trial", 1:10)


## @knitr 
f <- system.file("DESCRIPTION", package="UsingR")
dname <- dirname(f)
fname <- basename(f)


## @knitr 
paste(dname, fname, sep=.Platform$file.sep)


## @knitr 
require(MASS)
man <- Cars93$Manufacturer
length(man)                             # number of cases
length(levels(man))                     # number of levels


## @knitr 
cyl <- Cars93$Cylinders
levels(cyl)                             # "rotary"
which(cyl == "5")                       # just 5 is also okay
Cars93$Manufacturer[ which(cyl == 5) ]  # which companies


## @knitr 
mtcars$am <- factor(mtcars$am, labels=c("automatic", "manual"))


## @knitr 
require(HistData)
any(Arbuthnot$Female > Arbuthnot$Male)


## @knitr 
A <- c(TRUE, FALSE, TRUE, TRUE)
!A


## @knitr 
A <- c(TRUE, FALSE, TRUE, TRUE)
B <- c(TRUE, FALSE, TRUE, TRUE)
!(A & B)
!A | !B


## @knitr 
names(precip[precip > 50])


## @knitr 
m <- mean(precip)
trimmed_m <- mean(precip, trim=0.25)
any(precip > m + 1.5 * trimmed_m)


## @knitr 
"ABCDE" == "ABCDE"
"ABCDE" < "ABCDEF"
"ABCDE" < "abcde"
"ZZZZZ" < "aaaaa"
"11" < "8"


## @knitr 
commutes <- c(17, 16, 20, 24, 22, 15, 21, 15, 17, 22)
commutes[commutes == 24] <- 18
max(commutes)
min(commutes)
mean(commutes)
sum(commutes >= 20)
sum(commutes < 18)/length(commutes)


## @knitr 
cds <- c(79, 74, 161, 127, 133, 210, 99, 143, 249, 249, 368, 302)
longmos <- c(1, 3, 5, 7, 8, 10, 12)
long <- cds[longmos]
short <- cds[-longmos]
mean(long)
mean(short)


## @knitr 
x <- c(0.57, 0.89, 1.08, 1.12, 1.18, 1.07, 1.17, 1.38, 1.441, 1.72)
names(x) <- 1990:1999


## @knitr 
diff(x)


## @knitr 
which(diff(x) < 0)


## @knitr 
diff(x)/x[-10] * 100


## @knitr eval=FALSE
## function(x) {
##   sum(x)/length(x)
## }


## @knitr eval=FALSE
## function(x) {
##   total <- sum(x)
##   n <- length(x)
##   total / n
## }


## @knitr 
our_mean <- function(x) {
  sum(x)/length(x)
}                                       
our_mean(c(1,2,4,6,8,0))                # try it out


## @knitr 
mean_distance <- function(x) {
  distances <- abs(x - mean(x))
  mean(distances)
}


## @knitr 
f <- function(x) {
  mean(x^2) - mean(x)^2
}
f(1:10)


## @knitr 
iseven <- function(x) x %%2 == 0


## @knitr 
isodd <- function(x) x%%2 == 1


## @knitr 
iseven <- function(x) {
  x <- as.integer(x)
  ans <- x %% 2 == 0
  setNames(ans, x)                      # add names
}
iseven(1:10)


## @knitr 
isprime <- function(x){
  !any(x %% 2:(x-1) == 0)
}


## @knitr wts
wts <- c(38, 43, 48, 61, 47, 24, 29, 48, 59, 24, 40, 27)


## @knitr 
sort(wts)


## @knitr echo=FALSE, cache=FALSE, warning=FALSE, fig.width=6, fig.height=4, out.width="3in", out.height="2in"
xkcd_dotplot(wts, xlab="Weights of four-year olds", pt.size=5)


## @knitr 
mean(wts)


## @knitr 
devs <- wts - mean(wts)
mean(devs)                              # 0 up to rounding errors


## @knitr echo=FALSE, out.width=thirdwidth, cache=FALSE
## illustrate the mean is a balance
## verbose way to make these figures
draw_center <- function() {
  xlim = c(-2.2,3.2)
  x = c(-2,2,3)

  plot(x,1 + 0*x, pch=16, cex=3,col="black",
       bty = "n",
       ylim = c(0,2),
       xlim = xlim,
       xlab = "",
       ylab = "",
       yaxt = "n"
       )

  lines(xlim,c(.5,.5),lwd = 2)
  delta = .2
  xbar = mean(x)
  lines(c(xbar - delta,xbar,xbar+delta),c(0,.5,0),lwd=2)
}

draw_left <- function() {
  xlim = c(-2.2,3.2)
  x = c(-2,2,3)
  ## tilt left (at 2)
  y=0 #8
  plot(NA, NA, pch=16, cex=3,col="black",
       bty = "n",
       ylim = c(0,2),
       xlim = xlim,
       xlab = "",
       ylab = "",
       yaxt = "n"
       )
  points(x,y+1 + 1/8*(x-2), pch=16, cex=3,col="black")
  
  lines(xlim,y+c(0,.5+1/8*(3-2)),lwd = 2)
  xbar = mean(x)
  delta = 0.2
  lines(c(2 - delta,2,2+delta),y+c(0,.5,0),lwd=2)
}

draw_right <- function() {
 xlim = c(-2.2,3.2)
  x = c(-2,2,3)
  ## tilt right
  y = 0 #4
  plot(NA, NA, pch=16, cex=3,col="black",
       bty = "n",
       ylim = c(0,2),
       xlim = xlim,
       xlab = "",
       ylab = "",
       yaxt = "n"
       )
  points(x,y + 1 - 1/6*x,, pch=16, cex=3,col="black")
  delta = 0.2  
  lines(xlim,y+c(1/2 + 1/3,0),lwd = 2)
  xbar = mean(x)
  lines(c(0 - delta,0,0+delta),y+c(0,.5,0),lwd=2)
}
draw_right()
draw_left()
draw_center()


## @knitr 
mean(wts, trim = 0.10)                  # trim 10% of both ends


## @knitr 
mean(exec.pay)
mean(exec.pay, trim = 0.10)


## @knitr Macdonell_height
w <- Macdonell$frequency / sum(Macdonell$frequency) # n_k/n
y <- Macdonell$height
sum(w*y)


## @knitr 
median(wts)


## @knitr 
n <- length(exec.pay); trim = 0.10
lo <- 1 + floor(n * trim)               # floor drops decimal values
hi <- n + 1 - lo
median(sort(exec.pay)[lo:hi])           # compare to median(exec.pay)


## @knitr mean_median_wealth, results="asis", echo=FALSE
yr <- c(1989,1992,1995,1998,2001,2004,2007,2010)
median_hshd_wealth <-c(79100, 75100, 81900, 95600, 106100, 107200, 126400, 77300)
mean_hshd_wealth <-c(313600, 282900, 300400, 377300, 487000, 517100, 584600, 498800)
d <- data.frame(Year=yr, Median=median_hshd_wealth, Mean=mean_hshd_wealth)
d$Ratio <- round(d$Mean/d$Median, 1)

cat(booktabs(d, 
             caption="Median and Mean household net worth in U.S. by year. Source Federal Reserve Board, 2010 SCF Chartbook (\\url{http://www.fas.org/sgp/crs/misc/RL33433.pdf}).", 
             subtop=NULL, 
             rnames=NULL, 
             label="tab:mean-median",
             colTypes=c("l","c","c","c")
             ))



## @knitr 
x <- 0:5                                # 0,1,2,3,4,5
length(x)                               # even.
mean(sort(x)[3:4])                      # median averages 3rd, 4th
median(x)                               # clearly in middle
quantile(x, 0.25)                       # 1 + .25(5) = 2.25.
quantile(x, seq(0, 1, by=0.2))          # quintiles
quantile(x)                             # quartiles are default


## @knitr 
fivenum(x)


## @knitr 2011_income_capital_gains
income <- c("90" = 110651, "95" = 155193, "99" = 366623,  
            "99.5" = 544792,  "99.9" = 1557090, "99.99" = 7969900)
income


## @knitr 
table(wts)
table(wts) == max(table(wts))
which(table(wts) == max(table(wts)))
as.numeric(names( which(table(wts) == max(table(wts))) ))


## @knitr 
range(wts)                              # minimum and maximum values
diff(range(wts))                        # the distance between


## @knitr var
var(wts)


## @knitr 
sum( (wts - mean(wts))^2 ) / (length(wts) - 1)


## @knitr 
hip_cost <- c(10500, 45000, 74100, NA, 83500, 86000, 38200, NA, 
              44300, 12500, 55700, 43900, 71900, NA, 62000)
range(hip_cost, na.rm=TRUE)
sd(hip_cost, na.rm=TRUE)


## @knitr 
z_score <- function(x) (x - mean(x))/sd(x)
z_score(wts)


## @knitr 
scale(wts)[,1]


## @knitr 
x <- c(54, 50, 79, 79, 51, 69, 55, 62, 100, 80)
z <- (x - mean(x))/sd(x)
x[z >= 1.28]


## @knitr 
mean(x) + 1.28 * sd(x)


## @knitr echo=FALSE, warning=FALSE,  out.width=doublewide
z1 <- (wts - mean(wts))/sd(wts)
epay <- exec.pay[50:75]
z2 <- (epay - mean(epay))/sd(epay)
xkcd_dotplot(z1, xlab="Weight of four-year olds", pt.size=4)
xkcd_dotplot(z2, xlab="Executive pay", pt.size=4)


## @knitr 
z <- (exec.pay - mean(exec.pay)) / sd(exec.pay)
out <- abs(z) > 3                       # 199 TRUE or FALSE values
sum(out) / length(z)                    # sum of logical


## @knitr cof
sd(exec.pay)/mean(exec.pay)             # coefficient of variation


## @knitr IQR
median(rivers)                          # center
IQR(rivers)                             # spread


## @knitr 
IQR(rivers)/sd(rivers)


## @knitr mad
mad(rivers)/sd(rivers)                  # 213/496


## @knitr 
ht <- kid.weights$height
mad(ht)/sd(ht)                          # 11.68/10.7


## @knitr echo=FALSE, warning=FALSE,  out.width=doublewide, cache=FALSE
set.seed(4000)
xkcd_dotplot(sample(UsingR::exec.pay, 20), xlab="Executive pay sample", pt.size=4)
xkcd_dotplot(subset(UsingR::kid.weights, age < 60 & age >= 48)$height, xlab="Heights of four-year olds", pt.size=4)


## @knitr 
skew <- function(x) {
  n <- length(x)
  z <- (x - mean(x)) / sd(x)
  sum(z^3) / n
}


## @knitr 
skew(exec.pay)
four_year_hts <- kid.weights[kid.weights$age %in% 48:59, "height"]
skew(four_year_hts)


## @knitr 
kurtosis <- function(x) {
  n <- length(x)
  z <- (x - mean(x)) / sd(x)
  sum(z^4)/n - 3
}
kurtosis(galton$parent)              # height of parents in data set


## @knitr echo=FALSE, warning=FALSE,  out.width=doublewide, cache=FALSE
set.seed(4000)
xz <- rnorm(20); xz <- (xz - mean(xz))/sd(xz)
yz <- rt(20, df=3); yz <- (yz - mean(yz))/sd(yz)
xkcd_stackeddotplot(yz, xz, xlab="Normal and long-tailed data")
##
xkcd_dotplot(galaxies, xlab="Galaxy velocities",  ylimits=c(.84, 1.15))


## @knitr bumpers
stem(bumpers)


## @knitr echo=FALSE, warning=FALSE,  out.width=doublewide, cache=FALSE
hist(faithful$waiting)
hist(faithful$waiting, probability=TRUE)


## @knitr eval=FALSE
## hist(faithful$waiting)


## @knitr 
bins <- seq(40, 100, by=5)


## @knitr 
x <- faithful$waiting
out <- cut(x, breaks=bins)
head(out)


## @knitr 
table(out)


## @knitr echo=FALSE, warning=FALSE,  out.width=doublewide
x0 <- rnorm(10000)                       # 100,000
out <- hist(x0, plot=FALSE,  breaks = "FD")
x <- out$breaks
y <- c(0, out$density)
plot.new()
plot.window(xlim=range(x), ylim=range(y))
for(i in 1:(length(x)-1)) {
  lines(c(x[i+1], x[i+1]), c(0, y[i+1]), col="gray80")
  lines(c(x[i],x[i], x[i+1]), c(y[i], y[i+1], y[i+1]))
}
axis(1)
xkcd_density(x0, title="A density plot of the data")



## @knitr out.width=doublewide, echo=FALSE
## this graphic is based on code that appears following.
b_hist <- hist(bumpers, plot=FALSE)
b_dens <- density(bumpers)
hist(bumpers, probability=TRUE,
     xlim = range(c(b_hist$breaks,  b_dens$x)),
     ylim = range(c(b_hist$density, b_dens$y)))
lines(b_dens, lwd=2)


## @knitr fig.keep="none"
plot( density(bumpers) )


## @knitr plotting_arguments, echo=FALSE, results="asis"
x <- "
\\code{xlim}: Set $x$ coordinate range.
\\code{ylim}: Set $y$ coordinate range.
\\code{xlab}: Set label for $x$ axis.
\\code{ylab}: Set label for $y$ axis.
\\code{main}: Set main title.
\\code{pch}: Adjust plot symbols (\\code{?pch}).
\\code{cex}: Adjust size of text and symbols on a graphic.
\\code{col}: Adjust color of objects drawn. (\\code{?colors}).
\\code{lwd}: Adjust width of lines drawn.
\\code{lty}: Adjust how line is drawn. Can be \\qcode{blank}, \\qcode{solid}, \\qcode{dashed}, \\qcode{dotted}, \\qcode{dotdash}, etc.
\\code{bty}: Adjust box type, if drawn. One of \\qcode{o}, \\qcode{l}, \\qcode{7}, \\qcode{c},  \\qcode{u}, or  \\qcode{]}.
"
out <- read.dcf(textConnection(x))
d <- data.frame(Argument=colnames(out), Description=out[1,], stringsAsFactors=FALSE)
cat(booktabs(d, 
         caption="Standard plotting arguments to modify a graphic", 
         subtop=NULL, 
         rnames=NULL, 
         label="tab:various-plot-options",
         colTypes=c("l", "p{3in}")
         ))


## @knitr 
b_hist <- hist(bumpers, plot=FALSE)
b_dens <- density(bumpers)


## @knitr hist_with_density, fig.keep="none"
hist(bumpers, probability=TRUE,
     xlim = range(c(b_hist$breaks,  b_dens$x)),
     ylim = range(c(b_hist$density, b_dens$y)))
lines(b_dens, lwd=2)


## @knitr overlay_functions, echo=FALSE, results="asis"
x <- "
\\code{points}: Add points to a graphic.
\\code{lines}: Add points connected by lines to a graphic.
\\code{abline}: Add a line of the form $a + bx$, $y=h$, or $x=v$.
\\code{text}: Add text to a graphic.
\\code{mtext}: Add text to margins of a graphic.
"
out <- read.dcf(textConnection(x))
d <- data.frame(Functions=colnames(out), Description=out[1,], stringsAsFactors=FALSE)
cat(booktabs(d, 
         caption="Functions to add layers to a graphic", 
         subtop=NULL, 
         rnames=NULL, 
         label="tab:various-plot-functions",
         colTypes=c("l", "p{3in}")
         ))


## @knitr echo=FALSE, warning=FALSE,  out.width=doublewide
boxplot(bumpers, horizontal=TRUE, main="Bumpers")
boxplot(kid.weights$weight, horizontal=TRUE, main="Weights")


## @knitr eval=FALSE
## boxplot(bumpers, horizontal=TRUE, main="Bumpers")


## @knitr echo=FALSE, warning=FALSE,  out.width=triplewide
qqnorm(rep(Macdonell$finger, Macdonell$frequency))
qqnorm(jitter(HistData::Galton$child, 5))
qqnorm(exec.pay, main="exec.pay")


## @knitr qqnorm_fingers, fig=FALSE, fig.keep="none"
x <- rep(Macdonell$finger, Macdonell$frequency)
qqnorm(x)


## @knitr fig.keep="none"
x <- jitter(HistData::Galton$child, factor=5)     # add noise
qqnorm(x)


## @knitr 
require(LearnEDA)
lengths <- beatles$time / 60
c(mean=mean(lengths), median=median(lengths), 
  longest=max(lengths), shortest=min(lengths))


## @knitr 
nk <- ChestSizes$count
yk <- ChestSizes$chest
n <- sum(nk)
wk <- nk/n
sum(wk * yk)


## @knitr 
x <- c(80,82,88,91,91,95,95,97,98,101,106,106,109,110,111)   
median(x)


## @knitr 
require(LearnEDA)
stem(farms$count)


## @knitr 
2.3 * 10^(-4)


## @knitr fig.keep="none"
hist(firstchi)                # looks like 25 or so
mean(firstchi)                # we were pretty close ...


## @knitr fig.keep="none"
hist(pi2000-.1, prob=TRUE)
lines(density(pi2000))      


## @knitr fig.keep="none"
hist(normtemp$temperature)    # looks like its 98.2 -- not 98.6
mean(normtemp$temperature)    


## @knitr fig.keep="none"
require(MASS)
hist(DDT)
boxplot(DDT)


## @knitr 
c(mean=mean(DDT), sd=sd(DDT))


## @knitr fig.keep="none"
x <- rep(ChestSizes$chest, ChestSizes$count)
hist(x)


## @knitr 
x <- as.numeric(paradise)               # as numbers
x <- x[!is.na(x)]                       # strip NA values


## @knitr 
names(state.area) <- state.abb  
state.area['NJ']
sum(state.area < state.area['NJ'])/50 * 100
sum(state.area < state.area['NY'])/50 * 100


## @knitr 
hist(state.area)                        # 50,000 cuts of last case
state.area[state.area > 5e5]


## @knitr 
sum(pi2000 <= 3)/length(pi2000) * 100
sum(pi2000 >= 5)/length(pi2000) * 100


## @knitr 
sum(rivers < 500) / length(rivers)
sum(rivers < mean(rivers)) / length(rivers)
quantile(rivers, 0.75)


## @knitr 
times <- nym.2002$time                  # easier to use
range(times)                            # looks like minutes
sum(times < 3*60)/length(times) * 100   # 2.6% beat 3 hours
quantile(times,c(.10, .25))             # 3:28 to 3:53
quantile(times,c(.90))                  # 5:31


## @knitr 
mean(rivers)
median(rivers)
mean(rivers, trim=.25)


## @knitr 
stem(islands)                           # quite skewed
c(mean=mean(islands),
  median=median(islands),
  trimmed=mean(islands,trim=0.25))


## @knitr 
(OBP['bondsba01']- mean(OBP)) / sd(OBP)


## @knitr 
z <- scale(x)[,1]                       # matrix notation
mean(z)                                 # basically 0
sd(z)


## @knitr 
c(mad=mad(exec.pay), IQR=IQR(exec.pay), sd=sd(exec.pay))


## @knitr 
amt <- npdb$amount
summary(amt)
sum(amt < mean(amt))/length(amt) * 100


## @knitr 
sd(rivers) / mean(rivers)


## @knitr 
ia_times <- diff(babyboom$running.time)
sd(ia_times) / mean(ia_times)


## @knitr 
skew(babyboom$wt)
ia_times <- diff(babyboom$running.time)
skew(ia_times)


## @knitr fig.keep="none"
hist(hall.fame$HR)


## @knitr fig.keep="none"
require(HistData)
chest <- rep(ChestSizes$chest, ChestSizes$count)
qqnorm(chest)


## @knitr fig.keep="none"
hist(cfb$AGE)


## @knitr 
skew(Cars93$Price) > skew(Cars93$MPG.highway)


## @knitr 
Mode <- function(x) {
  tbl <- table(x)
  ind <- which(tbl == max(tbl))
  vals <- names(ind)
  as(vals, class(x)[1])                    # unnecessary!
}


## @knitr 
x <- babies$smoke
x <- factor(x, labels=c("never", "now", "until current",
                        "once, quit", "unknown"))
table(x)


## @knitr 
out <- table(x)
prop <- 100 * out / sum(out)               # or prop.table(out)
round(prop, digits = 2)                    # format output


## @knitr fig.keep="none"
barplot(table(x), horiz=TRUE, main="Smoking data")


## @knitr barcharts, echo=FALSE, warning=FALSE,  out.width=doublewide
barplot(table(x), horiz=TRUE, main="Smoking data")
dotchart(table(x), main="Smoking data")


## @knitr fig.keep="none", eval=FALSE
## dotchart(table(x))


## @knitr 
bumps <- cut(bumpers, c(0, 1000, 2000, 3000, 4000))
table(bumps)


## @knitr 
summary(Cars93$Cylinders)


## @knitr 
chars <- unlist(strsplit(lorem, split=""))
table(chars)


## @knitr 
sort(table(chars))


## @knitr fig.keep="none"
require(MASS)
dotchart(table(Cars93$Cylinders))


