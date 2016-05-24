require(its)

addDimnames <- function(mat) {
  if(is.null(dimnames(mat))) {dimnames(mat) <- list(NULL,NULL)}
  if(is.null(dimnames(mat)[[1]])&(nrow(mat)>0)) {dimnames(mat)[[1]] <- 1:nrow(mat)}
  if(is.null(dimnames(mat)[[2]])&(ncol(mat)>0)) {dimnames(mat)[[2]] <- 1:ncol(mat)}
  return(mat)
}

test <- function(x) {
  if(!x) stop()
}

today <- as.POSIXct(format(Sys.time(),"%Y-%m-%d"))
mytimes <- seq.POSIXt(from=today,by="DSTday",length.out=10)
attr(mytimes,"tzone") <- ""
mat <- addDimnames(matrix(1:30,10,3))
dimnames(mat)[[2]] <- c("A","B","C")
its.format("%Y-%m-%d")
x <- its(mat,mytimes)
moretimes <- seq.POSIXt(from=today+1,by="DSTday",length.out=11)
more <- addDimnames(matrix(1:33,11,3))
dimnames(more)[[2]] <- c("A","B","C")
x2 <- its(more,moretimes)


##its-arith**********************************************************
##arith-methods------------------------------------------------------
test.arith <- function() {
  y1 <- x[,1]+x[,2]
  y2 <- its(x@.Data[,1,drop=FALSE]+x@.Data[,2,drop=FALSE],mytimes)
  test(all.equal(y1@.Data,y2@.Data))
  test(all.equal(y1@dates,y2@dates))
  y1 <- x[,1]+x[,2]*pi
  y2 <- its(x@.Data[,1,drop=FALSE]+x@.Data[,2,drop=FALSE]*pi,mytimes)
  test(all.equal(y1@.Data,y2@.Data))
  test(all.equal(y1@dates,y2@dates))
  y1 <- x[,1]+pi
  y2 <- its(x@.Data[,1,drop=FALSE]+pi,mytimes)
  test(all.equal(y1@.Data,y2@.Data))
  test(all.equal(y1@dates,y2@dates))
}

##extractor**********************************************************
test.extractor <- function() {
  y1 <- x[,1,dates=dates(x)[1:5]]
  y2 <- x[1:5,1]
  y3 <- y1+y2
  test(all(dates(y1)==dates(x)[1:5]))
  test(all(core(y1)==core(x)[1:5,1]))
}

##names**************************************************************
test.names <- function() {
  test(all(names(x)==dimnames(core(x))[[2]]))
  y1 <- x
  names(y1) <- letters[1:ncol(y1)]
  test(all(names(y1)==letters[1:ncol(y1)]))
}

##dates**************************************************************
test.dates <- function() {
  test(all(dates(x)==x@dates))
  y1 <- x
  dates(y1) <- moretimes[1:nrow(y1)]
  test(all(dates(y1)==moretimes[1:nrow(y1)]))
}

##core***************************************************************
test.core <- function() {
  test(all(core(x)==x@.Data))
  y1 <- x
  core(x) <- addDimnames(matrix(101:130,10,3))
  test(all(core(x)==addDimnames(matrix(101:130,10,3))))
}

##its-cumdif*********************************************************
##cumsum-method------------------------------------------------------
test.cumsum <- function() {
  foo <- cumsum(x)
  test(all.equal(foo@.Data[,1],cumsum(x@.Data[,1])))
  test(all.equal(foo@dates,mytimes))
}

##diff-method--------------------------------------------------------
test.diff <- function() {
  foo <- diff(cumsum(x))
  bar <- alignedIts(foo,x,print=FALSE)
  test(all.equal(bar[[1]],bar[[2]]))
}

##its-def************************************************************
##-Functions-
##is.its-function----------------------------------------------------
test.is.its <- function() {
  test(is.its(x))
  test(!is.its(x@.Data))
  test(!is.its(x@dates))
}

##as.its-function----------------------------------------------------
test.as.its <- function() {
  foo <- as.numeric(mat[,1,drop=F])
  class(foo) <- c("POSIXt","POSIXct")
  bar <- its(mat[,-1],foo)
  waz <- as.its(mat)
  test(all.equal(bar,waz))
}

##its-function-------------------------------------------------------
test.its.creation <- function() {
  test(all.equal(x@dates,mytimes))
  test(all.equal(x@.Data/mat,x@.Data/x@.Data))

  ##parameters
  ## years <- 100:105
  ## hoursecs <- 60*60
  ## regdaysecs <- 24*hoursecs
  ## monthdays <- c(28,29,30,31)
  ## monthsecs <- c(monthdays*regdaysecs,monthdays*regdaysecs-hoursecs,monthdays*regdaysecs+hoursecs)
  ## weeksecs <- 7*regdaysecs ##+hoursecs*c(-1,0,1)
  ## daysecs <-  regdaysecs+hoursecs*c(-1,0,1)
  ## its.format("%Y-%m-%d")
}

##***newIts
##newIts-from
test.new.its <- function() {
  its.format("%Y-%m-%d")
  mystarts <- c("2003-01-01","2002-12-31","2003-11-17","2004-10-27")
  myends <- c("2003-02-01","2003-12-31","2004-12-17","2004-11-01")

  for(ddd in mystarts) {
    TEST <- newIts(start=ddd)
    test(start(TEST)==ddd)
  }

  ##newIts-to
  for(ddd in myends) {
    TEST <- newIts(start="2002-11-17",end=ddd)
    test(end(TEST)==ddd)
  }

  for(i in 1:3) {
    TEST <- newIts(start=mystarts[i],end=myends[i])
    test((start(TEST)==mystarts[i])&(end(TEST)==myends[i]))
  }

  ##newIts-by
  now <- as.POSIXct(format(Sys.time(),"%Y-%m-%d"))
  its.end.date <- now+100*24*60*60

  day.range <- seq.POSIXt(from=now,to=its.end.date,by="DSTday")
  TEST <- newIts(end=format(its.end.date,"%Y-%m-%d"))
  test(all.equal(day.range,TEST@dates))

  month.range <- seq.POSIXt(from=as.POSIXct("2003-10-01"),to=as.POSIXct("2010-12-01"),by="months")
  TEST <- newIts(start="2003-10-01",end="2010-12-01",by="month")
  test(all.equal(month.range,TEST@dates))

  week.range <- seq.POSIXt(from=now,to=its.end.date,by="weeks")
  TEST <- newIts(end=format(its.end.date,"%Y-%m-%d"),by="week")
  test(all.equal(week.range,TEST@dates))

  ##newIts-ncol
  ncol(newIts(end=format(its.end.date,"%Y-%m-%d"),ncol=5))==5

  ##***extractIts permutations
  ##-weekday
  weekDaySelection <- c(0,6)
  nowt <- newIts(extract=TRUE,weekday=TRUE,select=weekDaySelection)
  test(length(nowt)==0)
  weekDaySelection <- 1:5
  TEST1 <- newIts(extract=TRUE,weekday=TRUE)
  TEST2 <- newIts(extract=TRUE,weekday=TRUE,select=weekDaySelection)
  TEST3 <- newIts()
  test(length(TEST1)>0)
  test(length(TEST2)==length(TEST1))
  test(length(TEST3)>length(TEST2))
  test(all(as.POSIXlt(TEST1@dates)$wday%in%weekDaySelection))
  test(all(weekDaySelection%in%as.POSIXlt(TEST1@dates)$wday))

  ##-find ("all","last","first")
  test(identical(newIts(extract=TRUE,select=0:6,period="week"),newIts()))
  test(identical(newIts(extract=TRUE,select=0:6,period="week",find="all"),newIts()))
  TEST1 <- newIts(extract=TRUE,period="week",find="first",partial=FALSE)
  TEST2 <- newIts(extract=TRUE,select=0,period="week")
  test(all(TEST1%in%TEST2))
  TEST1 <- newIts(extract=TRUE,period="week",find="last",partial=FALSE)
  TEST2 <- newIts(extract=TRUE,select=6,period="week")
  test(all(TEST1%in%TEST2))
  TEST1 <- newIts(weekday=TRUE,extract=TRUE,period="week",find="first",partial=FALSE)
  TEST2 <- newIts(weekday=TRUE,extract=TRUE,select=1,period="week")
  test(all(TEST1%in%TEST2))
  TEST1 <- newIts(weekday=TRUE,extract=TRUE,period="week",find="last",partial=FALSE)
  TEST2 <- newIts(weekday=TRUE,extract=TRUE,select=5,period="week")
  test(all(TEST1%in%TEST2))

  ##-period
  test(all(as.POSIXlt(newIts(extract=TRUE,period="week",find="first",partial=FALSE)@dates)$wday==0))
  test(all(as.POSIXlt(newIts(extract=TRUE,period="week",find="last",partial=FALSE)@dates)$wday==6))
  test(all(as.POSIXlt(newIts(extract=TRUE,period="month",find="first",partial=FALSE)@dates)$mday==1))
  test(all(as.POSIXlt(newIts(extract=TRUE,period="month",find="last",partial=FALSE)@dates)$mday%in%28:31))
  test(identical(newIts(extract=TRUE,period="week",select=0:6),newIts()))
  test(identical(newIts(extract=TRUE,period="month",select=1:31),newIts()))

  ##-partials
  TEST <- newIts(start="2003-11-18")
  TEST1 <- newIts(start="2003-11-18",period="week",find="first",extract=TRUE,partial=TRUE)
  TEST2 <- newIts(start="2003-11-18",period="week",find="first",extract=TRUE,partial=FALSE)
  TEST@dates[1]==TEST1@dates[1]
  test((nrow(TEST1)-1)==nrow(TEST2))

  ##-select
  for(i in 0:6) {
    TEST <- newIts(extract=TRUE,period="week",select=i)
    test(all(as.POSIXlt(TEST@dates)$wday==i))
  }

  for(i in 30:31) {
    TEST <- newIts(extract=TRUE,period="month",select=i)
    test(all(as.POSIXlt(TEST@dates)$mday==i))
  }
}

##its-disp***********************************************************
##plot-method--------------------------------------------------------
##create 5 sinusoids differing in phase by pi/6
test.plot.its <- function() {
  its.format("%Y-%m-%d %X")
  sintimes <- seq.POSIXt(from=Sys.time(),by=24*60*60,length.out=100)
  sintimes.num <- as.numeric(sintimes)
  sinmat <- addDimnames(matrix(NA,100,5))
  for(j in 1:5){sinmat[,j] <- sin((6*pi*sintimes.num/(sintimes.num[100]-sintimes.num[1]))-j*pi/6)}
  dimnames(sinmat)[[2]] <- LETTERS[1:5]
  sinx <- its(sinmat,sintimes)
  par(mfrow=c(3,3))
  ##line,point
  plot(sinx,type="p",main="Point")
  plot(sinx,type="l",main="Line")
  plot(sinx,type="b",main="Both")
  ##colour, width, type cycling
  plot(sinx,lwdvec=1:3,main="Width")
  plot(sinx,ltypvec=1:3,main="Type")
  plot(sinx,colvec=c(1,2,7),main="Colour")
  ##axis
  plot(sinx,format="%B",main="Label")
  plot(sinx,at=sintimes[c(1,100)],main="Position")
  ##NA handling
  sinx[10:20,] <- sinx[10:20,]*NA
  sinx[,2] <- sinx[,2]*NA
  plot(sinx,interp="n",main="NAs")
}

##print-method-------------------------------------------------------
print(x)
##its-file***********************************************************
##writecsvIts-function-----------------------------------------------
test.read.write.its <- function() {
  file <- tempfile()
  writecsvIts(x,file,col.names=FALSE)
  writecsvIts(x,file,row.names=FALSE,col.names=FALSE)
  writecsvIts(x,file)
  ##readcsvIts-function------------------------------------------------
  foo <- its(readcsvIts(file))
  y <- its(x)
  test(identical(as.numeric(foo@.Data),as.numeric(y@.Data)))
  test(all.equal(foo@dates,y@dates))
  test(identical(dimnames(foo),dimnames(y)))
}

##its-fin************************************************************
##accrueIts-function-------------------------------------------------
##test(all.equal(accrueIts(x)-lagIts(x)[-1,]/(365),accrueIts(x)*0))
##its-info***********************************************************
##summary-method-----------------------------------------------------
test.its.summary <- function() {
  foo <- summary(x)
  test(all.equal(as.numeric(foo[1,]),seq(1,21,10)))
  test(all.equal(as.numeric(foo[6,]),seq(10,30,10)))
  test(all.equal(as.numeric(foo[8,]),rep(10,3)))
}

##start-method-------------------------------------------------------
test.its.start <- function() {
  test(identical(start(x,format="%Y-%m-%d-%X"),format.POSIXct(mytimes[1],format="%Y-%m-%d-%X")))
  test(identical(start(x[2:nrow(x),],format="%Y-%m-%d-%X"),format.POSIXct(mytimes[2],format="%Y-%m-%d-%X")))
}

##end -method--------------------------------------------------------
test.its.end <- function() {
  test(identical(end(x,format="%Y-%m-%d %X"),format.POSIXct(mytimes[10],format="%Y-%m-%d %X")))
  test(identical(end(x[1:(nrow(x)-1),],format="%Y-%m-%d %X"),format.POSIXct(mytimes[nrow(x)-1],format="%Y-%m-%d %X")))
}
##its-join***********************************************************
##alignedIts-function---------------------------------------------
test.aligned.its <- function() {
  its.format("%Y-%m-%d")
  isub <- seq(1,9,2)
  xsub <- x[isub,]
  foo <- alignedIts(x,xsub,print=F)

  ##identical is not working here
  test(identical(foo[[1]],xsub))
  test(identical(foo[[2]],xsub))

  test(identical(foo[[1]]@dates,xsub@dates))
  test(identical(core(foo[[1]]),core(xsub)))
}

##appendIts-function-------------------------------------------------
test.append.its <- function() {
  its.format("%Y-%m-%d %X")
  xx <- its(mat,mytimes)
  ## these operations change the order of the attributes of the date
  ## after this, identical can't be used to compare series
  ## because the attributes order does not match
  later <- mytimes+366*24*60*60
  over <- mytimes+5*24*60*60

  xlate <- its(mat,later)
  xover <- its(mat,over)

  foo <- appendIts(xx,xlate)
  bar <- appendIts(xlate,xx)
  test(all.equal(foo,bar))
  test(all.equal(foo[1:10],xx))
  ## test(identical(foo[11:20],xlate))
  test(all.equal(foo[11:20],xlate))
  foo <- try(appendIts(x,xx[(2:(nrow(x)-1)),],but=FALSE),silent=TRUE)
  test(identical(grep("appendor data must extend appendee data",foo)>0,TRUE))
  foo <- try(appendIts(x,xover,but=FALSE),silent=TRUE)
  test(identical(grep("overlap data does not match",foo)>0,TRUE))
  foo <- try(appendIts(x,xover),silent=TRUE)
  test(foo=="Error in appendIts(x, xover) : overlap not allowed\n")
  dimnames(xlate)[[2]][1] <- "Z"
  foo <- try(appendIts(x,xlate),silent=TRUE)
  test(foo=="Error in appendIts(x, xlate) : names of the two inputs must match\n")
  ##10 cases
  ##    S1  E1  S2  E2
  ## 1   1   2   3   4
  ## 2   1   3   2   4
  ## 3   1   4   2   3
  ## 4   2   3   1   4
  ## 5   2   4   1   3
  ## 6   3   4   1   2
  ## 7   1   2   2   2
  ## 8   1   2   3   3
  ## 9   2   3   1   1
  ##10   1   2   1   1

  ## 1   1   2   3   4
  x1 <- xx[1:4,]
  x2 <- xx[5:8,]
  foo <- appendIts(x1,x2,but=FALSE,matchnames=FALSE)
  test(identical(xx[1:8,],foo))
  foo <- appendIts(x1,x2,but=FALSE,matchnames=TRUE)
  test(identical(xx[1:8,],foo))
  foo <- appendIts(x1,x2,but=TRUE,matchnames=FALSE)
  test(identical(xx[1:8,],foo))
  foo <- appendIts(x1,x2,but=TRUE,matchnames=TRUE)
  test(identical(xx[1:8,],foo))
  ## 2   1   3   2   4
  x1 <- xx[1:4,]
  x2 <- xx[3:6,]
  foo <- appendIts(x1,x2,but=FALSE,matchnames=FALSE)
  test(identical(xx[1:6,],foo))
  foo <- appendIts(x1,x2,but=FALSE,matchnames=TRUE)
  test(identical(xx[1:6,],foo))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=FALSE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=TRUE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
  ## 3   1   4   2   3
  x1 <- xx[1:4,]
  x2 <- xx[2:3,]
  foo <- try(appendIts(x1,x2,but=FALSE,matchnames=FALSE),silent=TRUE)
  test(identical(grep("appendor data must extend appendee data",foo)>0,TRUE))
  foo <- try(appendIts(x1,x2,but=FALSE,matchnames=TRUE),silent=TRUE)
  test(identical(grep("appendor data must extend appendee data",foo)>0,TRUE))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=FALSE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=TRUE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
  ## 4   2   3   1   4
  x1 <- xx[2:3,]
  x2 <- xx[1:4,]
  foo <- try(appendIts(x1,x2,but=FALSE,matchnames=FALSE),silent=TRUE)
  test(identical(grep("appendor data must extend appendee data",foo)>0,TRUE))
  foo <- try(appendIts(x1,x2,but=FALSE,matchnames=TRUE),silent=TRUE)
  test(identical(grep("appendor data must extend appendee data",foo)>0,TRUE))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=FALSE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=TRUE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
  ## 5   2   4   1   3
  x1 <- xx[2:4,]
  x2 <- xx[1:3,]
  foo <- appendIts(x1,x2,but=FALSE,matchnames=FALSE)
  test(identical(xx[1:4,],foo))
  foo <- appendIts(x1,x2,but=FALSE,matchnames=TRUE)
  test(identical(xx[1:4,],foo))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=FALSE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=TRUE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
  ## 6   3   4   1   2
  x1 <- xx[3:4,]
  x2 <- xx[1:2,]
  foo <- appendIts(x1,x2,but=FALSE,matchnames=FALSE)
  test(identical(xx[1:4,],foo))
  foo <- appendIts(x1,x2,but=FALSE,matchnames=TRUE)
  test(identical(xx[1:4,],foo))
  foo <- appendIts(x1,x2,but=TRUE,matchnames=FALSE)
  test(identical(xx[1:4,],foo))
  foo <- appendIts(x1,x2,but=TRUE,matchnames=TRUE)
  test(identical(xx[1:4,],foo))
  ## 7   1   2   2   2
  x1 <- xx[1:4,]
  x2 <- xx[4,]
  foo <- appendIts(x1,x2,but=FALSE,matchnames=FALSE)
  test(identical(xx[1:4,],foo))
  foo <- appendIts(x1,x2,but=FALSE,matchnames=TRUE)
  test(identical(xx[1:4,],foo))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=FALSE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=TRUE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
  ## 8   1   2   3   3
  x1 <- xx[1:4,]
  x2 <- xx[5,]
  foo <- appendIts(x1,x2,but=FALSE,matchnames=FALSE)
  test(identical(xx[1:5,],foo))
  foo <- appendIts(x1,x2,but=FALSE,matchnames=TRUE)
  test(identical(xx[1:5,],foo))
  foo <- appendIts(x1,x2,but=TRUE,matchnames=FALSE)
  test(identical(xx[1:5,],foo))
  foo <- appendIts(x1,x2,but=TRUE,matchnames=TRUE)
  test(identical(xx[1:5,],foo))
  ## 9   2   3   1   1
  x1 <- xx[2:4,]
  x2 <- xx[1,]
  foo <- appendIts(x1,x2,but=FALSE,matchnames=FALSE)
  test(identical(xx[1:4,],foo))
  foo <- appendIts(x1,x2,but=FALSE,matchnames=TRUE)
  test(identical(xx[1:4,],foo))
  foo <- appendIts(x1,x2,but=TRUE,matchnames=FALSE)
  test(identical(xx[1:4,],foo))
  foo <- appendIts(x1,x2,but=TRUE,matchnames=TRUE)
  test(identical(xx[1:4,],foo))
  ##10   1   2   1   1
  x1 <- xx[2:4,]
  x2 <- xx[2,]
  foo <- appendIts(x1,x2,but=FALSE,matchnames=FALSE)
  test(identical(xx[2:4,],foo))
  foo <- appendIts(x1,x2,but=FALSE,matchnames=TRUE)
  test(identical(xx[2:4,],foo))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=FALSE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
  foo <- try(appendIts(x1,x2,but=TRUE,matchnames=TRUE),silent=TRUE)
  test(identical(grep("overlap not allowed",foo)>0,TRUE))
}

##union-method-------------------------------------------------------
test.union.its <- function() {
  its.format("%Y-%m-%d")
  isub <- seq(1,9,2)
  ioth <- setdiff(1:10,isub)
  xsub <- x[isub,]
  xun <- union(x,xsub)
  test(identical(xun[,1:3],x))
  test(identical(xun[isub,4:6],xsub))
  test(all(is.na(xun[ioth,4:6])))
}

##intersect-method---------------------------------------------------
test.intersect.its <- function() {
  its.format("%Y-%m-%d")
  isub <- seq(1,9,2)
  xsub <- x[isub,]
  xin <- intersect(x,xsub)
  test(identical(xin[,1:3],xsub))
  test(identical(xin[,4:6],xsub))
}

##its-lags***********************************************************
##lagIts-function----------------------------------------------------
test.lag.its <- function() {
  foo <- lagIts(x)
  test(all(foo[-1,]==x[-nrow(x),]))
  test(all(foo@dates==x@dates))
}

##lagdistIts-function------------------------------------------------
test.lagdist.its <- function() {
  foo <- lagdistIts(x[,1],1,3)
  test(all.equal(foo[,1],lagIts(x[,1],1)))
  test(all.equal(foo[,2],lagIts(x[,1],2)))
  test(all.equal(foo[,3],lagIts(x[,1],3)))
  test(all.equal(x@dates,foo@dates))
}

##its-subset*********************************************************
##rangeIts-function--------------------------------------------------
test.its.range <- function() {
  its.format("%Y-%m-%d")
  now <- format.POSIXct(Sys.time(),format=its.format())
  tomorrow <- format.POSIXct(Sys.time()+24*60*60,format=its.format())
  test(identical(appendIts(rangeIts(x,start=tomorrow),rangeIts(x,end=now)),x))
  test(nrow(rangeIts(x,start=now,end=now))==1)
  test(nrow(rangeIts(x,start=now,end=tomorrow))==2)
}

test.subset.its <- function() {
  ##[-method-----------------------------------------------------------
  i1 <- rep(c(TRUE,FALSE),5)
  j1 <- c(TRUE,FALSE,TRUE)
  test(all(x[2:8,2:3]==mat[2:8,2:3]))
  test(all(x[3,3]==mat[3,3]))
  test(all(x[0,0]==mat[0,0]))
  test(all(x[-2,]==mat[-2,]))
  test(all(x[,-2]==mat[,-2]))
  test(all(x[,j1]==mat[,j1]))
  test(all(x[,j1]==mat[,j1]))
  test(all(x[i1,j1]==mat[i1,j1]))
  test(all(mat[c(TRUE,FALSE,TRUE),]==x[c(TRUE,FALSE,TRUE),]))
  test(all(mat[,c(TRUE,FALSE)]==x[,c(TRUE,FALSE)]))
  mat2 <- mat
  x2 <- x
  mat2[c(TRUE,FALSE),c(TRUE,FALSE)] <- c(1000,2000)
  x2[c(TRUE,FALSE),c(TRUE,FALSE)] <- c(1000,2000)
  test(all(mat2==x2))
}

##its-times**********************************************************
##daysecondIts-function----------------------------------------------
##weekdayIts-function------------------------------------------------
test.times.its <- function() {
  foo <- as.POSIXlt(mytimes)$wday
  test(identical((foo>0 & foo<6),weekdayIts(x)))
}

##collapseIts--------------------------------------------------------
test.collapse.its <- function() {
  foo <- newIts(ncol=4,period="week",find="last",extract=T)[1:5,]
  foo[,1] <- c(NA,NA,NA,1,1)
  foo[,2] <- 1
  foo[,4] <- c(NA, 1, 1,1,1)
  dimnames(foo)[[2]] <- c("A","A","B","A")
  test(all.equal(collapseIts(foo),foo[,2:3]))
  foo[5,1] <- 1.01
  bar <- try(collapseIts(foo),silent=TRUE)
  test(bar=="Error in collapseIts(foo) : column data must match in collapse function\n")
}

##its-utilities******************************************************
######################################################################################################################################################################################
##fromirtsIts-function------------------------------------------------
##identical(fromirtsIts(irts(x@dates,x)),x)

##locf---------------------------------------------------------------
test.locf.its <- function() {
  foo <- x
  foo[2:4,] <- NA
  test(identical(dates(x),dates(foo)))
}

##-Utility Methods-
##validity check-----------------------------------------------------
##dates<--method-----------------------------------------------------
##[-method-----------------------------------------------------------

##-Utility Functions-
##addDimnames-function-----------------------------------------------
##gap.its-function---------------------------------------------------
##overlaps.its-function----------------------------------------------
##overlapmatches.its-function----------------------------------------
##namesmatch.its-function--------------------------------------------
##its.format-function------------------------------------------------
##simdates.its-function----------------------------------------------
##***extractDates
test.extract.its <- function() {
  its.format("%Y-%m-%d")
  years <- 100:105
  hoursecs <- 60*60
  regdaysecs <- 24*hoursecs
  monthdays <- c(28,29,30,31)
  monthsecs <- c(monthdays*regdaysecs,monthdays*regdaysecs-hoursecs,monthdays*regdaysecs+hoursecs)
  weeksecs <- 7*regdaysecs ##+hoursecs*c(-1,0,1)
  daysecs <-  regdaysecs+hoursecs*c(-1,0,1)
  weekDaySelection <- 1:5

  TEST <- newIts(start="2003-11-17",end="2005-12-25")
  ##-select
  test(all((as.numeric(extractIts(TEST,period="week",select=2)@dates)-as.numeric(extractIts(TEST,period="week",select=1)@dates))
           %in% daysecs))
  ##-weekday
  test(all(as.POSIXlt(extractIts(TEST,weekday=TRUE)@dates)$wday
           %in%1:5))
  test(all(as.POSIXlt(extractIts(TEST,weekday=TRUE,find="last")@dates)$wday
           %in%1:5))
  test(all(as.POSIXlt(extractIts(TEST,weekday=TRUE,select=weekDaySelection)@dates)$wday
           %in%weekDaySelection))
  test(all(as.POSIXlt(extractIts(TEST,weekday=TRUE,select=weekDaySelection,period="week",find="last")@dates)$wday==5))
  test(all(as.POSIXlt(extractIts(TEST,weekday=TRUE,select=weekDaySelection,period="week",find="first")@dates)$wday==1))
  ##-find
  test(all(as.POSIXlt(extractIts(TEST,weekday=TRUE,period="week",find="first")@dates[-1])$wday==1))
  TESTX <-extractIts(TEST[1:(length(TEST@dates)-2)],weekday=TRUE,period="week",find="last")@dates
  test(all(as.POSIXlt(TESTX)$wday[-length(TESTX)]==5))
  ##-period
  test(all(as.POSIXlt(extractIts(TEST,weekday=FALSE,period="year",find="first",partial=FALSE)@dates)$yday==0))
  test(all(as.POSIXlt(extractIts(TEST,weekday=FALSE,period="month",find="last",partial=FALSE)@dates)$mday%in%monthdays))
  test(all(as.POSIXlt(extractIts(TEST,weekday=FALSE,period="month",find="first",partial=FALSE)@dates)$mday==1))
  test(all(as.POSIXlt(extractIts(TEST,weekday=FALSE,period="week",find="first",partial=FALSE)@dates)$wday==0))
  test(all(as.POSIXlt(dates(extractIts(newIts(start="2001-12-21",end="2002-01-10")[-10:-11,],find="last",period="week",partial=FALSE)))$wday ==6))
  ##firstlast
  test(start(extractIts(TEST,per="y",find="f",firstlast=TRUE))==start(TEST)&&end(extractIts(TEST,,per="y",find="f",firstlast=TRUE))==end(TEST))
  ##-select
  test(all((as.numeric(extractIts(TEST,period="week",select=2)@dates)-as.numeric(extractIts(TEST,period="week",select=1)@dates))%in%daysecs))
  test(all(as.POSIXlt(extractIts(TEST,period="week",select=2)@dates)$wday==2))
  test(all(as.POSIXlt(extractIts(TEST,period="week",select=2,weekday=TRUE)@dates)$wday==2))
}

testIts <- function() {
  test.arith()
  test.extractor()
  test.names()
  test.dates()
  test.core()
  test.cumsum()
  test.diff()
  test.is.its()
  test.as.its()
  test.its.creation()
  test.new.its()
  ##test.plot.its()
  test.read.write.its()
  test.its.summary()
  test.its.start()
  test.its.end()
  test.aligned.its()
  test.append.its()
  test.union.its()
  test.intersect.its()
  test.lag.its()
  test.lagdist.its()
  test.its.range()
  test.subset.its()
  test.times.its()
  test.collapse.its()
  test.locf.its()
  test.extract.its()
}

cat(paste("* its test suite successful  *\n",R.version.string,"\n"))

##debug(testIts)
print(system.time(testIts()))
