### R code from vignette source 'datatable-intro.Rnw'

###################################################
### code chunk number 1: datatable-intro.Rnw:15-19
###################################################
if (!exists("data.table",.GlobalEnv)) library(data.table)  
# In devel won't call library, but R CMD build/check will.
rm(list=as.character(tables()$NAME),envir=.GlobalEnv)
# for development when we Sweave this file repeatedly. Otherwise first tables() shows tables from last run


###################################################
### code chunk number 2: datatable-intro.Rnw:36-38
###################################################
DF = data.frame(x=c("b","b","b","a","a"),v=rnorm(5))
DF


###################################################
### code chunk number 3: datatable-intro.Rnw:41-43
###################################################
DT = data.table(x=c("b","b","b","a","a"),v=rnorm(5))
DT


###################################################
### code chunk number 4: datatable-intro.Rnw:47-49
###################################################
CARS = data.table(cars)
head(CARS)


###################################################
### code chunk number 5: datatable-intro.Rnw:53-54
###################################################
tables()


###################################################
### code chunk number 6: datatable-intro.Rnw:66-67
###################################################
sapply(DT,class)


###################################################
### code chunk number 7: datatable-intro.Rnw:83-85
###################################################
tables()
DT


###################################################
### code chunk number 8: datatable-intro.Rnw:90-92
###################################################
DT[2,]         # select row 2
DT[x=="b",]    # select rows where column x == "b"


###################################################
### code chunk number 9: datatable-intro.Rnw:98-99
###################################################
cat(try(DT["b",],silent=TRUE))


###################################################
### code chunk number 10: datatable-intro.Rnw:103-105
###################################################
setkey(DT,x)
DT


###################################################
### code chunk number 11: datatable-intro.Rnw:113-114
###################################################
tables()


###################################################
### code chunk number 12: datatable-intro.Rnw:119-120
###################################################
DT["b",]


###################################################
### code chunk number 13: datatable-intro.Rnw:125-127
###################################################
DT["b",mult="first"]
DT["b",mult="last"]


###################################################
### code chunk number 14: datatable-intro.Rnw:132-133
###################################################
DT["b"]


###################################################
### code chunk number 15: datatable-intro.Rnw:138-148
###################################################
grpsize = ceiling(1e7/26^2)   # 10 million rows, 676 groups
tt=system.time( DF <- data.frame(
  x=rep(LETTERS,each=26*grpsize),
  y=rep(letters,each=grpsize),
  v=runif(grpsize*26^2),
  stringsAsFactors=FALSE)
)
head(DF,3)
tail(DF,3)
dim(DF)


###################################################
### code chunk number 16: datatable-intro.Rnw:157-160
###################################################
tt=system.time(ans1 <- DF[DF$x=="R" & DF$y=="h",])   # 'vector scan'
head(ans1,3)
dim(ans1)


###################################################
### code chunk number 17: datatable-intro.Rnw:165-167
###################################################
DT = as.data.table(DF)       # but normally use fread() or data.table() directly, originally 
system.time(setkey(DT,x,y))  # one-off cost, usually


###################################################
### code chunk number 18: datatable-intro.Rnw:169-173
###################################################
ss=system.time(ans2 <- DT[list("R","h")])   # binary search
head(ans2,3)
dim(ans2)
identical(ans1$v, ans2$v)


###################################################
### code chunk number 19: datatable-intro.Rnw:175-176
###################################################
if(!identical(ans1$v, ans2$v)) stop("vector scan vs binary search not equal")


###################################################
### code chunk number 20: datatable-intro.Rnw:186-189
###################################################
system.time(ans1 <- DT[x=="R" & y=="h",])   # works but is using data.table badly
system.time(ans2 <- DF[DF$x=="R" & DF$y=="h",])   # the data.frame way
mapply(identical,ans1,ans2)


###################################################
### code chunk number 21: datatable-intro.Rnw:204-206
###################################################
identical( DT[list("R","h"),],
           DT[.("R","h"),])


###################################################
### code chunk number 22: datatable-intro.Rnw:208-209
###################################################
if(!identical(DT[list("R","h"),],DT[.("R","h"),])) stop("list != . check")


###################################################
### code chunk number 23: datatable-intro.Rnw:229-230
###################################################
DT[,sum(v)]


###################################################
### code chunk number 24: datatable-intro.Rnw:235-236
###################################################
DT[,sum(v),by=x]


###################################################
### code chunk number 25: datatable-intro.Rnw:241-246
###################################################
ttt=system.time(tt <- tapply(DT$v,DT$x,sum)); ttt
sss=system.time(ss <- DT[,sum(v),by=x]); sss
head(tt)
head(ss)
identical(as.vector(tt), ss$V1)


###################################################
### code chunk number 26: datatable-intro.Rnw:248-249
###################################################
if(!identical(as.vector(tt), ss$V1)) stop("by check failed")


###################################################
### code chunk number 27: datatable-intro.Rnw:257-262
###################################################
ttt=system.time(tt <- tapply(DT$v,list(DT$x,DT$y),sum)); ttt
sss=system.time(ss <- DT[,sum(v),by="x,y"]); sss
tt[1:5,1:5]
head(ss)
identical(as.vector(t(tt)), ss$V1)


###################################################
### code chunk number 28: datatable-intro.Rnw:264-265
###################################################
if(!identical(as.vector(t(tt)), ss$V1)) stop("group check failed")


