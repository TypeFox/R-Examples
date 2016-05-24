## Read Life Expectancy table
## Return imputed life expectancies

getLifeExp <-
function(){
  if (any(is.na(DALYget("LE"))))
    stop("Life Expectancy table contains empty cell(s).",
         call. = FALSE)
  if (!all(grepl("^[[:digit:]]*\\.?[[:digit:]]+", DALYget("LE"))))
    stop("Life Expectancy table contains non-numeric value(s).",
         call. = FALSE)

  lineM <- numeric(0)
  lineF <- numeric(0)
  ages <- DALYget("ages")

  for (i in seq(20)){
    lineM <- c(lineM, seq(from = as.double(DALYget("LE", i, 1)),
                          to = as.double(DALYget("LE", i + 1, 1)),
                          length.out = 10 * (ages[i+1] - ages[i]) + 1))
    lineF <- c(lineF, seq(from = as.double(DALYget("LE", i, 2)),
                          to = as.double(DALYget("LE", i + 1, 2)),
                          length.out = 10 * (ages[i+1] - ages[i]) + 1))
    lineM <- lineM[-length(lineM)]
    lineF <- lineF[-length(lineF)]
  }

  lineM <- c(lineM, as.double(DALYget("LE", 21, 1)))
  lineF <- c(lineF, as.double(DALYget("LE", 21, 2)))

  return(c(lineM, lineF))
}