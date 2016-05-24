### R code from vignette source 'RUnit.Rnw'

###################################################
### code chunk number 1: RUnit.Rnw:251-292 (eval = FALSE)
###################################################
## library(RUnit)
## 
## ##  define sample functions to be tested
## foo <- function(x) {
##   x <- x*x
##   x <- 2*x
##   return(x)
## }
## test.foo <- function() {
## 
##   checkTrue(is.numeric(foo(1:10)))
##   checkEquals(length(foo(1:10)), 10)
##   checkEqualsNumeric(foo(1), 2)
## }
## 
## bar <- function(x, y=NULL) {
## 
##   if (is.null(y)) {
##     y <- x
##   }
## 
##   if (all(y > 100)) {
##     ##  subtract 100
##     y <- y - 100
##   }
## 
##   res <- x^y
##   return(res)
## }
## 
## track <- tracker(); ## initialize a tracking "object"
## track$init(); ## initialize the tracker
## a <- 1:10
## d <- seq(0,1,0.1)
## 
## resFoo <- inspect(foo(a), track=track); ## execute the test function and track
## resBar <- inspect(bar(d), track=track); ## execute the test function and track
## 
## resTrack <- track$getTrackInfo(); ## get the result of Code Inspector (a list)
## 
## printHTML.trackInfo(resTrack) ; ## create HTML sites


###################################################
### code chunk number 2: RUnit.Rnw:312-321 (eval = FALSE)
###################################################
## foo <- function(x)
## {
##   y <- 0
##   for(i in 1:x)
##   {
##       y <- y + x
##   }
##   return(y)
## }


###################################################
### code chunk number 3: RUnit.Rnw:324-337 (eval = FALSE)
###################################################
## foo.mod <- function(x)
## {
##    track$bp(1) ;
##    y <- 0
##    track$bp(2);
##    for(i in 1:x)
##   {
##       track$bp(4) ;
##      y <- y +x
##   }
##    track$bp(6);
##    return(y)
## }


###################################################
### code chunk number 4: RUnit.Rnw:342-345 (eval = FALSE)
###################################################
## if(any(a==1)) {
##   print("do TRUE")
## } else print ("do FALSE");


###################################################
### code chunk number 5: RUnit.Rnw:348-355 (eval = FALSE)
###################################################
## if(any(a==1)) {
##    track$bp(2);
##    print("do TRUE")
## }else{
##  track$bp(3);
##  print("do FALSE");
## }


