modeInterval <-
function(intervals,type='T') {
	vetorMin =intervals$minValue
	vetorMax =intervals$maxValue
	numberLines= length(vetorMin)
	getBreaks <- function(nclasses=1,...){
		breaks= NULL
		if(type=="T"){
			numbers= c(vetorMax,vetorMin)
			temp= table(numbers)
			name = names(temp)
			breaks = sort(as.numeric(name))
		}else{
			minValue = min(vetorMin)
			maxValue = max(vetorMax)
			
			
			range = maxValue-minValue
			classRange = range/nclasses
			breaks= minValue
			
			for (i in 1:nclasses) {
				if(i != nclasses){
					breaks = c(breaks,breaks[i]+classRange)
				}else{
					
					breaks = c(breaks,maxValue)
				}
			}
			
		}
	
		breaks
	}
#	probablity calculation for any case
	calcProbInt <- function(bks,...){
		vectResBks = NULL
		value = 0
		lenBks= length(bks)
		
		j = 1
		for(i in 1:(lenBks-1)){
			while(j<numberLines+1)	{
				
				
				if((vetorMin[j]<=bks[i] && vetorMax[j]>=bks[i+1])){
					denominator = vetorMax[j]- vetorMin[j]
					numerator = bks[i+1]-bks[i]
					valueTemp = numerator/denominator
					value= value + valueTemp/numberLines
					
					j=j+1
				}else if(vetorMin[j]>=bks[i] && vetorMin[j]<bks[i+1] && vetorMax[j]>=bks[i+1]){
					tamInter = vetorMax[j]- vetorMin[j]
					tamBks = bks[i+1] - vetorMin[j]
					valueTemp = tamBks/tamInter
					value= value + valueTemp/numberLines
					
					j=j+1
				}else if(vetorMin[j]<=bks[i] && vetorMax[j]>bks[i] && vetorMax[j]<=bks[i+1]){
					tamInter = vetorMax[j]- vetorMin[j]
					tamBks = vetorMax[j] - bks[i]
					valueTemp = tamBks/tamInter
					value= value + valueTemp/numberLines
					
					j=j+1
				}else if((vetorMin[j]>=bks[i] &&vetorMax[j]<=bks[i+1])){
					denominator = 1
					numerator = 1
					valueTemp = numerator/denominator
					value= value + valueTemp/numberLines
					
					j=j+1
				}else{
					
					j= j+1
				}		
			}
			vectResBks = c(vectResBks,value)
			
			value = 0
			j = 1
		}
	
		vectResBks
	}
	bks = NULL
	x=NULL
	if(length(vetorMin) != length(vetorMax) ){
		print("the length of the two vectors must be the same")
		stopifnot(vetorMin) == length(vetorMax)
	}else if(type == "T"){
		bks=getBreaks()
		x = calcProbInt(bks)	
	}else if(type == "ST"){
		nclasses = ceiling(log2(numberLines)+1)
		bks=getBreaks(nclasses = nclasses)
		x = calcProbInt(bks)
	}else if(type == "SC"){
		intervalTemp = interval(vetorMin,vetorMax)
		nclasses = (3.5* sdInterval(intervalTemp))/(numberLines^(1/3))
		bks=getBreaks(nclasses = ceiling(nclasses))
		x = calcProbInt(bks)
	}
	temp= max(x)
	return = NULL
	for (i in 1:length(x)) {
		if(x[i]== temp){
			return = interval(bks[i],bks[i+1])
		
		}
	}
	return
}

