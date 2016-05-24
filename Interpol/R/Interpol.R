Interpol <-
function(data,dims,method="linear") {

	## catch missing or bad input
	if(missing(data))
          stop("data argument missing")
        else
          #if(!is.vector(data))
          #  stop("invalid argument: data argument is required to be a vector of D features")
        if(!is.numeric(dims) | (dims < 2))
          stop("invalid argument: target dimension is required to be a positive integer value > 1")
	if(!is.character(method) | !(method=="linear" | method=="spline" | method=="natural" | method =="fmm" | method =="periodic" | method =="average"))
	  stop("invalid argument: method is required to be linear, spline, natural, fmm or periodic")
	if(missing(method))
		method = "linear"

	max.length = 0
	for(m in 1:length(data)){
		y = unlist(data[[m]])
		if(max.length < length(y)){
			max.length = length(y)	
		}
	}
	################################averaging
	#if(method == "average" | dims < (max.length / 2)) {
	if(method == "average") {
		output <- c()
		length = dims
		
		for(n in 1:length(data)) {
			y = unlist(data[[n]])
			x = seq(1:length(y))
			factor = (length(y))/(length)
			steps = (1:length)*factor
			#print(steps)
			values = mean(y[1:floor(steps[1])])
			for(i in 1:(length(steps)-1)){
				tmp = mean(y[floor(steps[i]):floor(steps[i+1])])
				#print(paste("unten: ", floor(steps[i])))
				#print(paste("oben: ", floor(steps[i+1])))
				#print(paste("werte:", y[floor(steps[i]):floor(steps[i+1])]))
				#print(tmp)
				#print("-------------------------------")
				values = c(values, tmp)
			}			
			output <- rbind(output,values)
		}
		colnames(output) <- 1:length(output[1,])
		for(i in 1:length(output[1,])){
			colnames(output)[i] <- paste("V", i, sep="")
		}
		rownames(output) <- 1:length(output[,1])

		NAS = which(is.na(output))
		output[NAS] = output[NAS-1]
		return(output)
	}
	
	################################linear interpolation
	if(method == "linear") {

		output <- c()
		length = dims

 		if(length == 2) {
			for(n in 1:length(data)) { 
				y = unlist(data[[n]])
				cut1 = 1
				cut2 = floor(length(y) / 2)
				cut3 = length(y)
				value1 = mean(y[cut1:cut2])
				value2 = mean(y[cut2:cut3])
				output <- rbind(output, c(value1, value2)) 
			}
		}
		else {
			output <- c()
			for(n in 1:length(data)) {
				list <- unlist(data[[n]])
				if(!(length(list)==length)){

					factor = (length(list) - 1) / (length - 1)
					output.tmp = list[1]
					for(i in 1:(length-2)){
						normIndex = factor * i
						floorIndex = floor(normIndex)+1
						tmp = list[floorIndex]*(floorIndex+1-normIndex) + list[floorIndex+1] * (normIndex-floorIndex)
						output.tmp = c(output.tmp, tmp)			
					}
					output.tmp = c(output.tmp, list[length(list)])
				}else{output.tmp = list}
				output <- rbind(output, output.tmp)
			}
		}

		rownames(output) <- names(data)
		return(output)

	}
	################################various spline interpolations
	if(method == "spline") { method = "monoH.FC" } #replaces spline with monoH.FC, other methods are used as is.

	output <- c()
	length = dims
		
	for(n in 1:length(data)) {
		y = unlist(data[[n]])
		x = seq(1:length(y))
		f_fmm = splinefun(x,y, method=method)
		factor = length(y)/length
		output <- rbind(output,f_fmm((1:length)*factor))
	}
	names(output) <- names(data)
	return(output)

}

