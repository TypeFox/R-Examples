compare <-
function(x,seqlength){
	specimen = nrow(x) #Records number of specimen with the same species name
	A.VLF.count = rep(0,seqlength) # creates a vector of seqlength 0's
	G.VLF.count = rep(0,seqlength) # creates a vector of seqlength 0's
	C.VLF.count = rep(0,seqlength) # creates a vector of seqlength 0's
	T.VLF.count = rep(0,seqlength) # creates a vector of seqlength 0's
	single <- rep(0,seqlength) # creates a vector of seqlength 0's
	shared <- rep(0, seqlength) # creates a vector of seqlength 0's
	slength<-seqlength+2
	for(i in 3:slength)
	{
		for(n in 1:(specimen)){
			if(is.na(x[n,i]) == FALSE){
				if(x[n,i] == "A"){
					A.VLF.count[i-2] = A.VLF.count[i-2] + 1
				}
				else{
					if(x[n,i] == "C"){
						C.VLF.count[i-2] = C.VLF.count[i-2] + 1
					}
					else{
						if(x[n,i] == "G"){
							G.VLF.count[i-2] = G.VLF.count[i-2] + 1
						}
						else{
							if(x[n,i] == "T"){
							T.VLF.count[i-2] = T.VLF.count[i-2] + 1
							}
						}
					}
				}
			}
		}
		if(A.VLF.count[i-2] == 1)
		{
			single[i-2] = single[i-2] + 1
		}
		else{
			if(A.VLF.count[i-2] > 1){
				shared[i-2] = shared[i-2] + A.VLF.count[i-2]
			}
		}
		if(C.VLF.count[i-2] == 1)
		{
			single[i-2] = single[i-2] + 1
		}
		else{
			if(C.VLF.count[i-2] > 1){
				shared[i-2] = shared[i-2] + C.VLF.count[i-2]
			}
		}
		if(G.VLF.count[i-2] == 1)
		{
			single[i-2] = single[i-2] + 1
		}
		else{
			if(G.VLF.count[i-2] > 1){
				shared[i-2] = shared[i-2] + G.VLF.count[i-2]
			}
		}
		if(T.VLF.count[i-2] == 1)
		{
			single[i-2] = single[i-2] + 1
		}
		else{
			if(T.VLF.count[i-2] > 1){
				shared[i-2] = shared[i-2] + T.VLF.count[i-2]
			}
		}
	}
	return(rbind(single,shared))
}
