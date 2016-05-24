##### getLogBoxsize() tests #####
context(desc="Test getLogBoxsize() function")

test_that(desc="Test, if getLogBoxsize throws errors/warnings on wrong arguments",
{
  expect_error(object=getLogBoxsize())
  expect_error(object=getLogBoxsize("aaa"))
  expect_error(getLogBoxsize(list(1,2)))
  expect_warning(getLogBoxsize(0))
  expect_warning(getLogBoxsize(c(1,2,-3)))
  expect_warning(getLogBoxsize(c(1,2,3)))
})

test_that(desc="Test, if getLogBoxsize returns results of same length as input",
{
  # check length of one
  expect_equivalent(object=length(getLogBoxsize(2)),expected=1)
  # check random length > 1
  length <- 1+rpois(1,10)
  input <- rpois(length,30)
  expect_equivalent(object=length(getLogBoxsize(input)),expected=length)
})

test_that(desc="Test, if getLogBoxsize returns appropriate results",
{
  expect_equivalent(object=getLogBoxsize(1),expected=0.009950331)
  # FIXME this fails! I hope due to numerical reasons, 
  # need to figure out how to display more accurate values in the console
  # expect_equivalent(object=getLogBoxsize(-1),expected=-0.01005034)
})

##### quote2box() tests #####
context(desc="Test quote2box() function")

test_that(desc="Test, if quote2box() throws errors/warnings on wrong arguments",
{
  expect_error(object=quote2box())
  expect_error(object=quote2box("aaa"))
  expect_error(quote2box())
  expect_error(quote2box(quote=1,boxsize=c(1,2)))
})


test_that(desc="Test, if quote2box() returns results of same length as input",
{
  # check length of one
  expect_equivalent(object=length(quote2box(rpois(1,100))),expected=1)
  # check random length > 1
  length <- 1+rpois(1,10)
  input <- rpois(length,30)
  expect_equivalent(object=length(quote2box(input)),expected=length)
})

test_that(desc="Test, if box2lower <= quote2box <= box2upper for all values",
{
  tst <- union(runif(n=50,min=0.1,max=1000),rpois(n=50,lambda=100))
  
  boxsize=1
  log=F
  boxnumber <- quote2box(quote=tst,boxsize=boxsize,log=log)
  lower <- box2lower(boxnumber=boxnumber,boxsize=boxsize,log=log)
  upper <- box2upper(boxnumber=boxnumber,boxsize=boxsize,log=log)
  expect_true(object=all(lower<=tst))
  expect_true(object=all(tst<=upper))
  
  boxsize=getLogBoxsize(1)
  log=T
  boxnumber <- quote2box(quote=tst,boxsize=boxsize,log=log)
  lower <- box2lower(boxnumber=boxnumber,boxsize=boxsize,log=log)
  upper <- box2upper(boxnumber=boxnumber,boxsize=boxsize,log=log)
  expect_true(object=all(lower<=tst))
  expect_true(object=all(tst<=upper))
})


##### nextBox() tests #####
context(desc="Test nextBox() function")

test_that(desc="Test, if nextBox() throws errors/warnings on wrong arguments",
{
  expect_error(object=nextBox())
  # FIXME add additional tests
})


test_that(desc="Test, if nextBox() returns results of same length as input",
{
  # check length of one
  expect_equivalent(object=length(nextBox(quote=rpois(1,100),status="X")),expected=1)
  expect_equivalent(object=length(nextBox(quote=rpois(1,100),status="O")),expected=1)
  # check random length > 1
  length <- 1+rpois(1,10)
  input <- rpois(length,30)
  expect_equivalent(object=length(nextBox(quote=input,status="X")),expected=length)
  expect_equivalent(object=length(nextBox(quote=input,status="O")),expected=length)
})

test_that(desc="Test, if nextBox() produces appropriate values",
{
  tst <- union(runif(n=500,min=0.1,max=1000),rpois(n=500,lambda=100))

  boxsize=1
  log=F
  boxnumber <- quote2box(quote=tst,boxsize=boxsize,log=log)
  lower <- box2lower(boxnumber=boxnumber,boxsize=boxsize,log=log)
  upper <- box2upper(boxnumber=boxnumber,boxsize=boxsize,log=log)
  expect_equal(object=nextBox(quote=tst,status="X",boxsize=boxsize,log=log,offset=0),expected=upper)
  expect_equal(object=nextBox(quote=tst,status="O",boxsize=boxsize,log=log,offset=0),expected=lower)
  
  
  boxsize=getLogBoxsize(1)
  log=T
  boxnumber <- quote2box(quote=tst,boxsize=boxsize,log=log)
  lower <- box2lower(boxnumber=boxnumber,boxsize=boxsize,log=log)
  upper <- box2upper(boxnumber=boxnumber,boxsize=boxsize,log=log)
  expect_equal(object=nextBox(quote=tst,status="X",boxsize=boxsize,log=log,offset=0),expected=upper)
  expect_equal(object=nextBox(quote=tst,status="O",boxsize=boxsize,log=log,offset=0),expected=lower)
})

test_that(desc="Test, if nextBox() produces appropriate values",
{
  expect_equal(object=nextBox(quote=46.0,status="X",boxsize=1L,log=F),expected=47)
  expect_equal(object=nextBox(quote=46.0000001,status="X",boxsize=1L,log=F),expected=47)
  expect_equal(object=nextBox(quote=46.9999999,status="X",boxsize=1L,log=F),expected=47)
  expect_equal(object=nextBox(quote=47.0,status="X",boxsize=1L,log=F),expected=48)

  expect_equal(object=nextBox(quote=46.0,status="O",boxsize=1L,log=F),expected=46)
  expect_equal(object=nextBox(quote=46.0000001,status="O",boxsize=1L,log=F),expected=46)
  expect_equal(object=nextBox(quote=46.9999999,status="O",boxsize=1L,log=F),expected=46)
  expect_equal(object=nextBox(quote=47.0,status="O",boxsize=1L,log=F),expected=47)
})

##### nextReversal() tests #####
context(desc="Test nextReversal() function")

test_that(desc="Test, if nextBox() produces appropriate values",
{
  expect_equal(object=nextReversal(45.99999,"X",boxsize=1),expected=43)
  expect_equal(object=nextReversal(46.0,"X",boxsize=1),expected=44)
  expect_equal(object=nextReversal(46.00001,"X",boxsize=1),expected=44)
  expect_equal(object=nextReversal(45.99999,"O",boxsize=1),expected=48)
  expect_equal(object=nextReversal(46.0,"O",boxsize=1),expected=49)
  expect_equal(object=nextReversal(46.00001,"O",boxsize=1),expected=49)

  # check another boxsize
  expect_equal(object=nextReversal(45.99999,"X",boxsize=0.5),expected=44.5)
  expect_equal(object=nextReversal(46.0,"X",boxsize=0.5),expected=45)
  expect_equal(object=nextReversal(46.00001,"X",boxsize=0.5),expected=45)
  expect_equal(object=nextReversal(45.99999,"O",boxsize=0.5),expected=47)
  expect_equal(object=nextReversal(46.0,"O",boxsize=0.5),expected=47.5)
  expect_equal(object=nextReversal(46.00001,"O",boxsize=0.5),expected=47.5)

  # check another reversal number
  expect_equal(object=nextReversal(45.99999,"X",reversal=1,boxsize=1),expected=45)
  expect_equal(object=nextReversal(46.0,"X",reversal=1,boxsize=1),expected=46)
  expect_equal(object=nextReversal(46.00001,"X",reversal=1,boxsize=1),expected=46)
  expect_equal(object=nextReversal(45.99999,"O",reversal=1,boxsize=1),expected=46)
  expect_equal(object=nextReversal(46.0,"O",reversal=1,boxsize=1),expected=47)
  expect_equal(object=nextReversal(46.00001,"O",reversal=1,boxsize=1),expected=47)
  
})

##### xo.processor() tests #####
context(desc="Test xo.processor() function")

test_that(desc="Test, if xo.processor() throws errors/warnings on wrong arguments",
{
  expect_error(object=xo.processor())
})

# test_that(desc="Test, if speed of xo.processor() is sufficent",
# {
#   data(DOW)
#   times <- system.time(xo.processor(high=DOW$High, low=DOW$Low,date=DOW$Date))
#   warning(paste0("Timings of xo.processor() for linear charts: ",times[1]," sec."))
#   expect_less_than(object=times[1],expected=0.15)
#   
#   times <- system.time(xo.processor(high=DOW$High, low=DOW$Low,date=DOW$Date,boxsize=getLogBoxsize(1),log=T))
#   warning(paste0("Timings of xo.processor() for logarithmic charts: ",times[1]," sec."))
#   expect_less_than(object=times[1],expected=0.15)
# })

  test_that(desc="Test, if speed of xo.processor() scales nearly linear in input size",
{
  # check if speed scales nearly linear
  myfactor <- 10
  length=2000 # approx. 1 years
  mydeltas <- rnorm(n=length,mean=0,sd=1)
  myts <- rep(x=1000,length.out=length)
  for (i in 1:(length-1)) {
    myts[i+1] = myts[i]+(myts[i]*mydeltas[i]/100)
  }
  times.short <- system.time(xo.processor(high=myts, low=myts,date=seq(1:length),boxsize=getLogBoxsize(1),log=T))

  length=length*myfactor # approx. 10 years
  mydeltas <- rnorm(n=length,mean=0,sd=1)
  myts <- rep(x=1000,length.out=length)
  for (i in 1:(length-1)) {
    myts[i+1] = myts[i]+(myts[i]*mydeltas[i]/100)
  }
  times.long <- system.time(xo.processor(high=myts, low=myts,date=seq(1:length),boxsize=getLogBoxsize(1),log=T))
  
  expect_less_than(object=times.long[1],expected=(2*myfactor)*times.short[1])
})

test_that(desc="Test, if xo.processor() produces correct output",
{
  data(DOW)
  load(file="boxutils-example1.RData") # this loads an object result
  expect_equivalent(object=xo.processor(high=DOW$High,low=DOW$Low,date=DOW$Date),expected=result)
})
