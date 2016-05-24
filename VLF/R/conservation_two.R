conservation_two <-
function(modal1, modal2, p,seqlength){
	both <- modal1 + modal2
	first <- 0
	second <- 0
	third <- 0
	i <- 1
	while(i <= seqlength){
		if(both[i] >= p){
			first = first + 1
		}
		if(both[i+1] >= p){
			second = second + 1
		}
		if(both[i+2] >= p){
			third = third + 1
		}
		i = i+3
	}
	combined <- rbind(first,second,third)
	return(combined)
}
