timestamp =  function(dataset, Date, Time) {
  
  if(missing(Date))     {Date = "Date"}
  if(missing(Time))     {Time = "Time"}
  
  date = as.vector(dataset[,Date])
  TIME = as.vector(dataset[,Time])
  
  dateC = as.character(date)
  nc = nchar(dateC)
  dc = nc - 4
  sc = nc - 1
  dateA = substr(date, 1, dc)
  dateB = substr(date, sc, nc)
  date1 = paste(dateA, dateB, sep = "")
  date2 = as.Date(date1, "%m/%d/%y")
  TS = paste(date2, TIME)
  mydata_ts = data.frame(TS, dataset)
  
  return = mydata_ts
  }