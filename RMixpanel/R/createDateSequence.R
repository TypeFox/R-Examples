createDateSequence = function(
  from,       # Start date in either format <"yyyy-mm-dd"> or <yyyymmdd>. Inclusive.
  to=from     # End date in either format <"yyyy-mm-dd"> or <yyyymmdd>. Inclusive.
  )
### MP, 2015-09
### Generate Regular Sequences of daily Dates.
{
  from = as.character(from)
  to = as.character(to)
  
  if (nchar(from) == 8)
    from = paste(substring(from, c(1, 5, 7), c(4, 6, 8)), collapse="-")

  if (nchar(to) == 8)
    to = paste(substring(to, c(1, 5, 7), c(4, 6, 8)), collapse="-")
  
  ## Validity test.
  ## assertCondition(nchar(from) == 10 && nchar(to) == 10) #, "Date format is not valid. Pls Supply yyyymmdd or 'yyyy-mm-dd'.")
  
  ## Convert to series of class <date>
  dates = as.character(seq(as.Date(from), as.Date(to), "day"))
  dates
}
