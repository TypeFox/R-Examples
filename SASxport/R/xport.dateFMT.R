`xport.dateFMT` <-
  function(when, fill=16)
  {
    if(missing(when)) when <- Sys.time()
    
    date.format <- sprintf("%%-%d.%ds", fill, fill)
    toupper(sprintf(date.format, format(when, "%d%b%y:%H:%M:%S")))
  }

