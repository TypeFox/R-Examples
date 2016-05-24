"ptri" <-

function(x, a = 0., b = a + 2.)



{



	x <- (x - a)/(b - a)



	temp <- ifelse(x < 0.5, 2. * x^2., 1. - 2. * (1. - x)^2.)



	temp <- ifelse(x < 0., 0., temp)



	ifelse(x > 1., 1., temp)



}

