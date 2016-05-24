ConvertTimeFormat_KML <-
function(sTime,iAddHours,iAddDays)
{
  sOrigTime <- as.POSIXct(sTime)
  sNewTime <- sOrigTime + (iAddHours * 60 * 60) + (iAddDays * 60 * 60 * 24)
  # adding time gives a format like this;
  # "2007-09-18 00:15:34 EST"
  #  12345678901234567890123
  # we need a format like this;
  # "2005-08-28T05:21:00+10:00"
  #  1234567890123456789012345
  
  sTimeString <- substring(sNewTime,12,19)
  if (sTimeString == "")
    sTimeString <- "00:00:00"
  
  paste(substring(sNewTime,1,10),"T",sTimeString,sep="")
}
