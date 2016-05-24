boa.getinput <- function(message, n = 1, evaluate = TRUE)
{
   repeat {
      cat(message)
      input <- scan(what = "", n = n, sep = "\n")
      if(evaluate) value <- try(eval(parse(text = input)), TRUE)
      else value <- try(parse(text = input), TRUE)
      if(inherits(value, "try-error")) cat("Warning: syntax error\n")
      else break
   }

   value
}
