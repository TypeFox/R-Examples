# based on "singlechoice 15_09_29.R"

singlechoice <- function(val1, val2, k) {

  # normal approach
  # p_win <- pnorm((val1 - val2)/(200 * sqrt(2)))

  # logistic approach
  p_win <- 1 - 1/(1 + 10^((val1 - val2)/400))

  kp <- k * p_win

  val1new <- val1 - kp + k
  val2new <- val2 + kp - k

  return(round(c(val1new, val2new)))

}

