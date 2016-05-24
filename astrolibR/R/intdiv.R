intdiv = function(dividend, divisor) {
   a = as.integer(dividend)
   b = as.integer(divisor)

   sign(a) * (abs(a) %/% b)
}

