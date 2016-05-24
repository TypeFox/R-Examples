# -----------------------------------------------
# Function: Write data to clipboard.
# -----------------------------------------------
toclipboard <- function(data,
                  sep = "\t",
                  quote = FALSE,
                  eol = "\n",
                  na = "NA",
                  dec = ".",
                  row.names = FALSE,
                  col.names = TRUE) {
  
  if (tolower(.Platform$OS.type)=="windows") {
  write.table(data, file = "clipboard", sep = sep,
              quote = quote, eol = eol, na = na, dec = dec,
              row.names = row.names, col.names = col.names)
  } else {
  print("The operating system you are using is currently not supported by 'toclipboard'. No output has been written.")
  }
  
}