## helper function for 'sensitivity()'
## extract unique input parameters from 'DALY' object

as_column <-
function(x, outcome, name){
  ## check if M==F
  ## check if ages equal
  sex_strat <- !identical(x[, , 1], x[, , 2])
  age_strat <- !identical(x[, 1, ], x[, 2, ], x[, 3, ], x[, 4, ], x[, 5, ])

  ## empty data.frame
  data <- matrix(x, ncol = 10)
  colnames(data) <-
    paste(paste(name, outcome, sep = ""),
          rep(c("M", "F"), each = 5), rep(1:5, 2), sep = ".")
  data <- as.data.frame(data)

  ## fill data.frame
  if (!sex_strat){
    data <- data[, 1:5]
    colnames(data) <- paste(paste(name, outcome, sep = ""), 1:5, sep = ".")

    if (!age_strat){
      data <- as.data.frame(data[, 1])
      colnames(data) <- paste(paste(name, outcome, sep = ""), sep = ".")
    }

  } else {
    if (!age_strat){
      data <- data[, c(1, 6)]
      colnames(data) <-
        paste(paste(name, outcome, sep = ""), c("M", "F"), sep = ".")
    }
  }

  ## return data.frame
  return(data)
}