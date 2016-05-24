"qtri" <-

function(p, a = 0., b = a + 2.)



{



	ifelse(p < 0.5, sqrt(p/2.), 1. - sqrt((1. - p)/2.)) * (b - a) +



		a



}

