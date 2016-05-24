`qgpd2` <-
function(N , sigma = 1 , xi = 1 , u = 0, la = 1 ) # returns N observation return level
	u + ( sigma * ( ( N * la )^( xi ) - 1)) / xi

