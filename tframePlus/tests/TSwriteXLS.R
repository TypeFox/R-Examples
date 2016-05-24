require("tframePlus")

 z <- ts(1:10, start=c(1999,2), freq=12)
 seriesNames(z) <- "ser 1"
 ok <- TSwriteXLS(z, FileName="test.xls") 

 zz <- tbind(z, diff(z))
 seriesNames(zz) <- c("ser 1", "diff")
 ok <- TSwriteXLS(zz, FileName="test.xls",  SheetNames="2 series")

 zz <- ts(1:10, start=c(1999,1), freq=1)
 seriesNames(zz) <- "annual ser"
 ok <- TSwriteXLS(z, zz, FileName="test.xls",  SheetNames=c("monthly", "annual"))

 unlink("test.xls") 
