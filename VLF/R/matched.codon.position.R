matched.codon.position <-
function(matched){
	ones <- 0
	twos <- 0
	threes <- 0
	for(i in 1:length(matched)){
		if(matched[[i]][5] == 1){
			ones <- ones + 1
		}
		else{
			if(matched[[i]][5] == 2){
				twos <- twos + 1
			}
			else{
				if(matched[[i]][5] == 3){
					threes <- threes + 1
				}
			}
		}
	}
	return(rbind(ones,twos,threes))
}
