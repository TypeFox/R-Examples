AkuData = function(n)
{
  #data("testDataBasis",envir=environment())
  testDataBasis = getTestDataBasis()
  ix=sample(1:dim(testDataBasis)[1],size=100+round(n/0.4),replace=T,prob=testDataBasis[,9])
  x=testDataBasis[ix,1:8]
  sstat = (x %% 1000) - 100
  sstat[sstat>500] = NaN
  sstat[apply(!is.na(sstat),2,cumsum)>n] = NaN
  rows = rowSums(!is.na(sstat))>0
  x = x[rows,]
  sstat = as.vector(sstat[rows,])
  id=as.vector(row(x))
  age = as.vector(x %/% 1000)
  R = as.numeric(sstat>0)
  q     = as.vector(1+ (col(x)-1) %%4)
  year  = as.vector(2014+ (col(x)-1) %/%4)
  month  = as.vector(((row(x)-1) %%3)) + 1 + (q-1)*3
  edu   = rep(sample(1:4,size=dim(x)[1],replace=T,prob=c(3,6,4,2)),dim(x)[2])
  sex   = rep(sample(c(0,1),size=dim(x)[1],replace=T),dim(x)[2])
  id = 1:dim(x)[1]
  famid  = rep(sample(1:60,size=dim(x)[1],replace=T),dim(x)[2]) + 100*(id %/%100)
  z=data.frame(id,year,q,month,R,age,sex,edu,famid)
  z =z[is.finite(sstat),]
  sstat = sstat[is.finite(sstat)]
  rownames(z) = NULL
  z$year      = factor(z$year)
  z$q         = factor(z$q)
  z$month     = factor(z$month)
  z$age      = factor(z$age)
  z$sex      = factor(z$sex)
  z$edu      = factor(z$edu)
  z$unemployed  = as.numeric(sstat==200)
  z$workforce   = as.numeric(sstat==100 | sstat==200)
  z
}

# stackoverflow questions 30357330
pkgEnvAkuData <- new.env(parent=emptyenv())
  if(!exists("testDataBasis", pkgEnvAkuData)) {
  data("testDataBasis", package="CalibrateSSB", envir=pkgEnvAkuData)
}

getTestDataBasis <- function() {
  pkgEnvAkuData[["testDataBasis"]]
}
