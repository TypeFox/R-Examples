## stackListItems.R
## Paul Johnson <pauljohn at ku.edu>
## 2010-09-07

## I asked about this in r-help last week and promised
## a summary of answers.

## We face this problem all the time. A procedure
## generates a list of data frames. How to stack them
## together?

## The short answer is that the plyr package's rbind.fill
## method is probably the fastest method that is not
## prone to trouble and does not require much user caution.

## result <- rbind.fill(mylist)

## A slower alternative that also works is

## result <- do.call("rbind", mylist)

## That is always available in R and it works well enough, even
## though it is not quite as fast. Both of these are much faster than
## a loop that repeatedly applies "rbind".

## Truly blazing speed can be found if we convert this into
## matrices, but that is not possible if the list actually
## contains data frames.



## Here is a test case

df1 <- data.frame(x = rnorm(100), y = rnorm(100))
df2 <- data.frame(x = rnorm(100), y = rnorm(100))
df3 <- data.frame(x = rnorm(100), y = rnorm(100))
df4 <- data.frame(x = rnorm(100), y = rnorm(100))

mylist <-  list(df1, df2, df3, df4)

## Here's the way we have done it. We understand this,
## we believe the result, it is easy to remember. It is
## also horribly slow for a long list.

resultDF <- mylist[[1]]
for (i in 2:4) resultDF <- rbind(resultDF, mylist[[i]])


## It works better to just call rbind once, as in:

resultDF2 <- rbind( mylist[[1]], mylist[[2]], mylist[[3]], mylist[[4]])


## That is faster because it calls rbind only once.

## But who wants to do all of that typing? How tiresome.
## Thanks to Erik Iverson in r-help, I understand that

resultDF3 <- do.call("rbind", mylist)

## is doing the EXACT same thing.
## Erik explained that "do.call( "rbind", mylist)"
## is *constructing* a function call from the list of arguments.
## It is shorthand for "rbind(mylist[[1]], mylist[[2]], mylist[[3]])"
## assuming mylist has 3 elements.

## Check the result:
all.equal( resultDF2, resultDF3)

## You often see people claim it is fast to allocate all
## of the required space in one shot and then fill it in.
## I got this algorithm from code in the
## "complete" function in the "mice" package.
## It allocates a big matrix of 0's and
## then it places the individual data frames into that matrix.

m <- 4
nr <- nrow(df1)
nc <- ncol(df1)
resultDF4 <- as.data.frame(matrix(0, nrow = nr*m, ncol = nc))
for (j in  1:m) resultDF4[(((j-1)*nr) + 1):(j*nr), ] <- mylist[[j]]

## This is a bit error prone for my taste. If the data frames have
## different numbers of rows, some major code surgery will be needed.

##
## Dennis Murphy pointed out the plyr package, by Hadley Wickham.
## Dennis said " ldply() in the plyr package. The following is the same
## idea as do.call(rbind, l), only faster."

library("plyr")
resultDF5  <- ldply(mylist, rbind)
all.equal(resultDF, resultDF5)



## Plyr author Hadley Wickham followed up with "I think all you want here is rbind.fill:"

resultDF6 <- rbind.fill(mylist)
all.equal(resultDF, resultDF6)


## Gabor Grothendieck noted that if the elements in mylist were matrices, this would all work faster.

mylist2 <- lapply(mylist, as.matrix)

matrixDoCall <-  do.call("rbind", mylist2)

all.equal(as.data.frame(matrixDoCall), resultDF)


## Gabor also showed a better way than 'system.time' to find out how
## long this takes on average using the rbenchmark package. Awesome!

#> library(rbenchmark)
#> benchmark(
#+ df = do.call("rbind", mylist),
#+ mat = do.call("rbind", L),
#+ order = "relative", replications = 250
#+ )



## To see the potentially HUGE impact of these changes, we need to
## make a bigger test case. I just used system.time to evaluate, but
## if this involved a close call, I'd use rbenchmark.

phony <- function(i){
  data.frame(w = rnorm(1000), x = rnorm(1000), y = rnorm(1000), z = rnorm(1000))
}
mylist <- lapply(1:1000, phony)




### First, try my usual way
resultDF <- mylist[[1]]
system.time(
for (i in 2:1000) resultDF <- rbind(resultDF, mylist[[i]])
           )
## wow, that's slow:
## user  system elapsed
## 168.040   4.770 173.028


### Now do.call method:
system.time( resultDF3 <- do.call("rbind", mylist) )
all.equal(resultDF, resultDF3)

## Faster! Takes one-twelfth as long
##   user  system elapsed
##  14.64    0.85   15.49


### Third, my adaptation of the complete function in the mice
### package:
m <- length(mylist)
nr <- nrow(mylist[[1]])
nc <- ncol(mylist[[1]])

system.time(
   resultDF4 <- as.data.frame(matrix(0, nrow = nr*m, ncol = nc))
 )

colnames(resultDF4) <- colnames(mylist[[1]])
system.time(
   for (j in  1:m) resultDF4[(((j-1)*nr) + 1):(j*nr), ] <- mylist[[j]]
)

all.equal(resultDF, resultDF4)
##Disappointingly slow on the big case:
#   user  system elapsed
# 80.400   3.970  84.573


### That took much longer than I expected, Gabor's
### hint about the difference between matrix and data.frame
### turns out to be important. Do it again, but don't
### make the intermediate storage thing a data.frame:
mylist2 <- lapply(mylist, as.matrix)

m <- length(mylist2)
nr <- nrow(mylist2[[1]])
nc <- ncol(mylist2[[1]])

system.time(
   resultDF4B <- matrix(0, nrow = nr*m, ncol = nc)
 )

colnames(resultDF4B) <- colnames(mylist[[1]])
system.time(
   for (j in  1:m) resultDF4B[(((j-1)*nr) + 1):(j*nr), ] <- mylist2[[j]]
)

### That's FAST!
###    user  system elapsed
###   0.07    0.00    0.07

all.equal(resultDF, as.data.frame(resultDF4B))



### Now the two moethods from plyr.


system.time( resultDF5  <- ldply(mylist, rbind))

## Just about as fast, much less error prone
##  user  system elapsed
##  1.290   0.000   1.306

all.equal(resultDF, resultDF5)


system.time(resultDF6 <- rbind.fill(mylist))
##   user  system elapsed
##  0.450   0.000   0.459

all.equal(resultDF, resultDF6)



## Gabor was right. If we have matrices, do.call is
## just about as good as anything.

system.time(matrixDoCall <-  do.call("rbind", mylist2) )
##   user  system elapsed
##  0.030   0.000   0.032


all.equal(as.data.frame(matrixDoCall), resultDF)



