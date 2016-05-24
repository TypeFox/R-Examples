# flat pattern filters

# filter.function: "iqr", "rsd", "sd", "mad", "mean", "median"
# by.percent - if T the number of variables to filter will be a percentage of the number of 
# variables in the dataset; percentage is given by the "red.value" parameter
# if red.value = "auto", will calculate number of variables to remove automatically 
# (as implemented in metaboanalyst)
# by.threshold - if T, defines filtering will select variables where values of filtering
# function are below a given threshold
# threshold is defined by red.value that defines the minimum value of the function 
# needed to keep the variable

"flat_pattern_filter" = function(dataset, filter.function = "iqr", by.percent = T, 
                                 by.threshold = F, red.value = 0){
	
  if (by.percent & by.threshold) 
    warning("Both by.percent and by.threshold are T; filtering by percentage")
  if (!by.percent & !by.threshold) 
    stop ("Either by.percent or by.threshold need to be T")
    
  # apply filter function to all variables
  filter.values = apply_filter_function(dataset$data, filter.function)
  
  if (by.percent)
    dataset$data = flat_pattern_filter_percentage(dataset$data, filter.values, red.value)
  else if(by.threshold)
    dataset$data = flat_pattern_filter_threshold(dataset$data, filter.values, red.value)
  
	add.desc = paste("Flat pattern filtering with function", filter.function, sep=" ")
	dataset$description = paste(dataset$description, add.desc, sep="; ")
	dataset
}


# FLAT PATTERN FILTERING

# method: iqr, rsd, rnsd, sd, mad, mean, median
"apply_filter_function" = function(datamat, filter.fn = "iqr"){
	if (filter.fn == "iqr"){
		filter.values = apply(datamat, 1, IQR, na.rm=T);
	}
  else if (filter.fn == "rsd"){
		sds = apply(datamat, 1, sd, na.rm=T)
		mns = apply(datamat, 1, mean, na.rm=T)
		filter.values = abs(sds/mns)
	}
  else if (filter.fn == "rnsd"){
		mads = apply(datamat, 1, mad, na.rm=T)
		meds = apply(datamat, 1, median, na.rm=T)
		filter.values = abs(mads/meds)
	}
  else if (filter.fn == "sd"){
		filter.values = apply(datamat, 1, sd, na.rm=T)
	}
  else if (filter.fn == "mad"){
		filter.values = apply(datamat, 1, mad, na.rm=T)
	}
  else if (filter.fn == "mean"){
		filter.values = apply(datamat, 1, mean, na.rm=T)
	}
  else if (filter.fn == "median"){
		filter.values = apply(datamat, 1, median, na.rm=T)
	}
	filter.values
}

"flat_pattern_filter_percentage" = function(datamat, filter.values, percentage = "auto") {
  
  rk = rank(-filter.values, ties.method="first")
	var.num = nrow(datamat)
	if (percentage == "auto") {
		if (var.num < 250) datamat = datamat[rk < var.num*0.95,]
    else if (var.num < 500) datamat = datamat[rk < var.num*0.90,]
		else if (var.num < 1000) datamat = datamat[rk < var.num*0.75,]
		else datamat = datamat[rk < var.num*0.6,]
	} 
  else if (is.numeric(percentage) ) {
		var.remain = rk < var.num*((100-percentage)/100)
		datamat = datamat[var.remain,]
	}
  else stop("Invalid value for percentage parameter")
	datamat
}


"flat_pattern_filter_threshold" = function(datamat, filter.values, threshold = 0) {
  var.remain = filter.values >= threshold
  datamat[var.remain,]
}
