msg_format <-
  "<Status code: %s>
  Message: %s
  Explanation: %s"

msg_list <- function(a, b, c){
  list(status_code = a,
       message = b,
       explanation = c
  )
}

is_three_digit_code <- function(x) grepl("\\d{3}", x)

nozero <- function(x) names(x[sapply(x, length) > 0])

stopcode <- function(x, y) stop(sprintf('%s: %s\n', x, y), call. = FALSE)
