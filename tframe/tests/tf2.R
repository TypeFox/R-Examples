 require("stats")
 require("tframe")

 Sys.info()

t1<-ts(c(1,2,3,4,5),start=c(1991,1))
t2<-ts(c(2,3,4,5,6,7,8),start=c(1992,1))
t3<-ts(c(NA,2,3,4,5),start=c(1991,1))

latestStart(t1,t2,t3) # 1992 1 corresponding to the starting date of 
                       # the object which starts latest (t2)
ok <- all(c(1992,1) == latestStart(t1,t2,t3))
ok

latestStart(t1,t3)     # both start in 1991 1 (NAs count as data)
ok <- ok & all(c(1991,1) == latestStart(t1,t3))
ok

latestStart(tbind(t1,t2,t3)) # tbind gives a single object starting in 1991 1
ok <- ok & all(c(1991,1) == latestStart(tbind(t1,t2,t3)))
ok

latestStart(t2, tbind(t1,t2,t3))
ok <- ok & all(c(1992,1) == latestStart(t2, tbind(t1,t2,t3)))
ok

latestStartIndex(t1,t2,t3)  # position of t2 in the argument list
ok <- ok & 2 == latestStartIndex(t1,t2,t3)
ok

latestEndIndex(t1,t2,t3)  # position of t2 in the argument list
ok <- ok & 2 == latestEndIndex(t1,t2,t3)
ok

earliestEndIndex(t1,t2,t3)  # position of t2 in the argument list
ok <- ok & 1 == earliestEndIndex(t1,t2,t3)
ok

earliestStart(t1,t2) # 1991 1
ok <- ok & all(c(1991,1) == earliestStart(t1,t2))
ok

earliestEnd(t1,t2)   # 1995 1
ok <- ok & all(c(1995,1) == earliestEnd(t1,t2))
ok

latestEnd(t1,t2)     # 1998 1
ok <- ok & all(c(1998,1) == latestEnd(t1,t2))
ok

if (! ok) stop("some tests FAILED")
