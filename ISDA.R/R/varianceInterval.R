varianceInterval <-
function(interval){
	xmin= interval$minValue
	xmax= interval$maxValue
	m =length(xmin)
	variance = 0
	if(m != length(xmax)){
		print("the lenght of xmin and xmax must be the same")
	}else{
		term1temp =0
		term2temp =0
		for(i in 1:m){
			term1temp = term1temp + ((xmax[i]^2) +(xmin[i]*xmax[i])+(xmin[i]^2))
			term2temp = term2temp+ (xmin[i]+xmax[i])
		}
		term1temp = (1/(3*m))*term1temp
		term2temp = (1/(4*(m^2)))*(term2temp^2)
		variance = term1temp-term2temp
	}
	variance
}

