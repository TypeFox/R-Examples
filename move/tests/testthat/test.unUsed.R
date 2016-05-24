context("unUsedRecords")
test_that('timestamps ordering in stack',{
  expect_is(u<-new('.unUsedRecordsStack', 
                 timestampsUnUsedRecords=Sys.time()+1:3,   
                 sensorUnUsedRecords=factor(rep(1,3)),
                 trackIdUnUsedRecords=factor(c('a','a','b')),
                 dataUnUsedRecords=data.frame(1:3)
  ),'.unUsedRecordsStack')
  t<-as.POSIXct(rep(1,3), origin='1970-1-1', tz='UTC')
  u@timestampsUnUsedRecords<-t
  expect_error(validObject(u),"The data set includes double timestamps per ID in the unused records \\(first one:a 1 1970-01-01 00:00:01\\)")
  t<-as.POSIXct(c(2,1,1), origin='1970-1-1')
  u@timestampsUnUsedRecords<-t
  expect_error(validObject(u),"The data set includes un ordered timestamps in the unUsedRecordsStack")
  
  })
context("examples")
test_examples()