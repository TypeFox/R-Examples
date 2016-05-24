"dtri" <-

function(x, a = 0., b = a + 2.)



{



	x <- (x - a)/(b - a)



	# zero out density beyond range of support.



	temp <- (2. * ifelse(x < 0.5, x, 1. - x))/(b - a)



	ifelse(temp < 0., 0., temp)



}

