## ----setup, echo = FALSE, message=FALSE, warning=FALSE, results='hide'----
knitr::opts_chunk$set(
  cache=FALSE,
  comment="#>",
  collapse=TRUE, 
  echo=TRUE
)
library(knitr); library(eeptools)

## ------------------------------------------------------------------------
age_calc(dob = as.Date('1995-01-15'), enddate = as.Date('2003-02-16'), 
         units = "years")
age_calc(dob = as.Date('1995-01-15'), enddate = as.Date('2003-02-16'), 
         units = "months")
age_calc(dob = as.Date('1995-01-15'), enddate = as.Date('2003-02-16'), 
         units = "days")

## ------------------------------------------------------------------------
x <- data.frame(sid = c(101, 101, 102, 103, 103, 103, 104, 105, 105, 106, 106),
                 grade = c(9, 10, 9, 9, 9, 10, 10, 8, 9, 7, 7))
retained_calc(df = x, sid = "sid", grade = "grade", grade_val = 9)

## ------------------------------------------------------------------------
df <- data.frame(sid = c(rep(1,3), rep(2,4), 3, rep(4,2)),
                   schid = c(1, 2, 2, 2, 3, 1, 1, 1, 3, 1),
                   enroll_date = as.Date(c('2004-08-26',
                   '2004-10-01', '2005-05-01', '2004-09-01',
                   '2004-11-03', '2005-01-11', '2005-04-02',
                   '2004-09-26', '2004-09-01','2005-02-02'), format='%Y-%m-%d'),
                   exit_date = as.Date(c('2004-08-26', '2005-04-10',
                    '2005-06-15', '2004-11-02', '2005-01-10',
                    '2005-03-01', '2005-06-15', '2005-05-30',
                    NA, '2005-06-15'), format='%Y-%m-%d'))

moves <- moves_calc(df, sid = "sid", schid = "schid", enroll_date = "enroll_date", 
                    exit_date = "exit_date")
moves


## ----statamode-----------------------------------------------------------
vecA <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
statamode(vecA, method = "stata")
vecB <- c(1, 1, 1, 3:10)
statamode(vecB, method = "last")
vecC <- c(1, 1, 1, NA, NA, 5:10)
statamode(vecC, method = "last")
vecA <- c(LETTERS[1:10]); vecA <- factor(vecA)
statamode(vecA, method = "last")
vecB <- c("A", "A", "A", LETTERS[3:10]); vecB <- factor(vecB)
statamode(vecB, method = "last")
vecA <- c(LETTERS[1:10])
statamode(vecA, method = "sample")
vecB <- c("A", "A", "A", LETTERS[3:10])
statamode(vecB, method = "stata")
vecC <- c("A", "A", "A", NA, NA, LETTERS[5:10])
statamode(vecC, method = "stata")

## ----sim-----------------------------------------------------------------
library(MASS)
#Examples of "sim" 
set.seed (1)
J <- 15
n <- J*(J+1)/2
group <- rep (1:J, 1:J)
mu.a <- 5
sigma.a <- 2
a <- rnorm (J, mu.a, sigma.a)
b <- -3
x <- rnorm (n, 2, 1)
sigma.y <- 6
y <- rnorm (n, a[group] + b*x, sigma.y)
u <- runif (J, 0, 3)
dat <- cbind (y, x, group)
# Linear regression 
dat <- as.data.frame(dat)
dat$group <- factor(dat$group)
M3 <- glm (y ~ x + group, data=dat)
cases <- expand.grid(x = seq(-2, 2, by=0.1), 
                     group=seq(1, 14, by=2))
cases$group <- factor(cases$group)
sim.results <- gelmansim(mod = M3, newdata = cases, n.sims=200, na.omit=TRUE)
head(sim.results)

## ----map-----------------------------------------------------------------
crimes <- data.frame(state = tolower(rownames(USArrests)), USArrests)
require(reshape) # for melt
crimesm <- melt(crimes, id = 1)
states_map <- map_data("state")
p1 <- ggplot(crimes, aes(map_id = state)) + geom_map(aes(fill = Murder), 
                                                     linetype = 1, map = states_map) + 
       expand_limits(x = states_map$long, y = states_map$lat) + labs(title="USA Crime")
p1 <- p1 + coord_map()
p1 + theme_dpi_map()

## ----lmautoplot----------------------------------------------------------
data(mpg)
mymod <- lm(cty~displ + cyl + drv, data=mpg)
autoplot(mymod)

## ----crossplot-----------------------------------------------------------
sampDat <- data.frame(cbind(x=seq(1,3,by=1), y=sample(LETTERS[6:8], 60, 
                                                        replace=TRUE)),
                        fac=sample(LETTERS[1:4], 60, replace=TRUE))
varnames<-c('Quality','Grade')
crosstabplot(sampDat, "y", "fac", varnames = varnames,  label = TRUE, 
             title = "Crosstab Plot", shade = FALSE)
crosstabplot(sampDat, "y", "fac", varnames = varnames,  label = FALSE, 
             title = "Crosstab Plot", shade = TRUE)

## ------------------------------------------------------------------------
library(eeptools)
data("stuatt")
head(stuatt)

## ------------------------------------------------------------------------
data(stulevel)
head(stulevel)


## ------------------------------------------------------------------------
data("midsch")
head(midsch)


