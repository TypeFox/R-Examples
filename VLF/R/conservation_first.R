conservation_first <-
function(modal, p,seqlength){
	first <- 0
	second <- 0
	third <- 0
	i <- 1
	while(i <= seqlength){
		if(modal[i] >= p){
			first = first + 1
		}
		if(modal[i+1] >= p){
			second = second + 1
		}
		if(modal[i+2] >= p){
			third = third + 1
		}
		i = i+3
	}
	combined <- rbind(first,second,third)
	return(combined)
}
