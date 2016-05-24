aa.compare <-
function(x, seqlength){
	specimen = nrow(x) #Records number of specimen with the same species name
	aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
	VLF.count = mat.or.vec(nr = 20, nc = seqlength) # creates a 4x216 matrix of 0's
	single <- rep(0,216) # creates a vector of 216 0's
	shared <- rep(0,216) # creates a vector of 216 0's

	for(i in 1:seqlength){
		for(n in 1:nrow(x)){
			if(is.na(x[n,i+2]) == FALSE){
				z = 1
				match = FALSE
				while(match == FALSE && z <= length(aa)){
					if(x[n, i+2] == aa[z]){
						VLF.count[z,i] = VLF.count[z,i] + 1
						match = TRUE
					}
					else{
						z = z + 1
					}
				}
			}
		}
	}
	
	for(i in 1:seqlength){
		for(z in 1:20){
			if(VLF.count[z,i] == 1)
			{
				single[i] = single[i] + 1
			}
			else{
				if(VLF.count[z,i] > 1){
					shared[i] = shared[i] + VLF.count[z,i]
				}
			}
		}
	}
	return(rbind(single,shared))
}
