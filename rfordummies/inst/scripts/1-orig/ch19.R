# Chapter 19 - Ten Things You Can Do in R That You Would’ve Done in Microsoft Excel

# Adding Row and Column Totals

iris.num <- iris[, -5]

colSums(iris.num)
colMeans(iris.num)

apply(iris.num, 2, min)
apply(iris.num, 2, max)

sapply(iris.num, min)
sapply(iris.num, max)


# Formatting Numbers

format(12345.6789, digits=9, decimal.mark=",",
   big.mark=" ",small.mark=".", , small.interval=3)

x <- colMeans(mtcars[, 1:4])
format(x, digits=2, nsmall=2)

x <- seq(0.5, 0.55, 0.01)
sprintf("%.1f %%", 100*x)

set.seed(1)
x <- 1000*runif(5)
sprintf("$ %3.2f", x)

stuff <- c("bread", "cookies")
price <- c(2.1, 4)
sprintf("%s costed $ %3.2f ", stuff, price)


# Sorting Data

with(mtcars, mtcars[order(hp), ])
with(mtcars, mtcars[order(hp, decreasing=TRUE), ])

# Making Choices with If

mtcars <- within(mtcars,
   mpgClass <- ifelse(mpg < mean(mpg), "Low", "High"))

mtcars[mtcars$mpgClass == "High", ]


# Calculating Conditional Totals

with(mtcars, mean(mpg))
with(mtcars, mean(mpg[hp < 150]))
with(mtcars, mean(mpg[hp >= 150]))
with(mtcars, length(mpg[hp > 150]))


# Transposing Columns or Rows

x <- matrix(1:12, ncol=3)
x
t(x)

t(mtcars[1:4, ])


# Finding Unique or Duplicated Values

unique(mtcars$cyl)
dupes <- duplicated(iris)
head(dupes)
which(dupes)
iris[dupes, ]
iris[!dupes, ]
nrow(iris[!dupes, ])


# Working with Lookup Tables

index <- match("Toyota Corolla", rownames(mtcars))
index
mtcars[index, 1:4]


# Working with Pivot Tables

with(mtcars, tapply(hp, list(cyl, gear), mean))
aggregate(hp~cyl+gear+am, mtcars, mean)


# Using the Goal Seek and Solver

sales <- function(price) { 100 - 0.5 * price }
revenue <- function(price) { price * sales(price) }


par(mfrow=c(1, 2))
curve(sales, from=50, to=150, xname="price", ylab="Sales", main="Sales")
curve(revenue, from=50, to=150, xname="price", ylab="Revenue", main="Revenue")
par(mfrow=c(1, 1))

optimize(revenue, interval=c(50, 150), maximum=TRUE)
