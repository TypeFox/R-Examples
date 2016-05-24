access <-
function(tickers, start, end, data)
{
  datelist = seq(as.Date(start), as.Date(end), by="day")
  dat = data[as.character(datelist), tickers]
  return(dat)
}
