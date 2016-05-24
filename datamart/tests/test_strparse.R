library(datamart)

test_strparse <- function() {
  pat <- "^[\\s]+(?<vval>.*)"
  print(strparse(pat, ""))
  print(strparse(pat, "Depends: xdata"))
  print(strparse(pat, "  xdata"))
  print(strparse(pat, "\txdata"))
  print(strparse(pat, "\t\txdata"))
  pat <- "^(?<vname>[^:\\s]+)[\\s]*:[\\s]*(?<vval>.*)$"
  print(strparse(pat, ""))
  print(strparse(pat, "Depends: xdata"))
  print(strparse(pat, "Depends  : xdata"))
  print(strparse(pat, "Depends\t:\txdata"))
  print(strparse(pat, "  xdata"))
  print(strparse(pat, "\txdata"))
  print(strparse(pat, "\t\txdata"))
  mlines <- c(
        'VARIABLE LABELS weight "weight".',
        'VARIABLE LABELS altq "Year of birth".',
        'VARIABLE LABELS hhg "Household size".',
        'missing values all (-1).',
        'EXECUTE.'
  )
  pat <- 'VARIABLE LABELS (?<name>[^\\s]+) \\"(?<lbl>.*)\\".$'
  print(strparse(pat, mlines))
}

test_strparse()

