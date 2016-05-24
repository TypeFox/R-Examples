hillshading = function(cgrad,sv){
	if (nargs() < 2) {
		cat("USAGE: hillshading(cgrad,sunvector) \n")
		return()
	}

	hsh = cgrad[,,1]*sv[1]+cgrad[,,2]*sv[2]+cgrad[,,3]*sv[3]
	hsh = (hsh + abs(hsh))/2. 
}
