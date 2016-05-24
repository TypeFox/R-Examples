G1 <-
function(s, a = 1)
###Cosh negentropy function + derivatives
list(Gs = logb(cosh(a * s))/a, gs = tanh(a * s), gps = a * (1 - tanh(a * s)^
	2))

