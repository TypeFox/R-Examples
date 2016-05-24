AddMetFiles <- 
  function(month, Year, met, script.file, control.file) 
{
  ## if month is one, need previous year and month = 12
  if (month == 0) {
    month <- 12
    Year <- as.numeric(Year) - 1
  }
  
  if (month < 10) {
    month <- paste("0", month, sep = "")
  }
  
  ## add first line
  line <- paste("echo", met, ">>", control.file, sep=" ")
  cat(line, file = script.file, sep = "\n")
  
  line <- paste("echo RP", Year, month, ".gbl >> ", control.file, sep = "")
  cat(line, file = script.file, sep = "\n")
}