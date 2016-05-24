"random.walk" <-

function(n, p = 0.5)



{



	cumsum(rbern(n, p) * 2. - 1.)



}

