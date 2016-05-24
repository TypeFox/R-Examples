dayOfYear <-
  function(dat) {       
    i <- as.numeric(format(as.Date(dat,origin="1970-01-01"),'%j')) 
    i
  }

