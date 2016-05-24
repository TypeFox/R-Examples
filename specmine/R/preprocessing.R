#convert transmittance to absorbance
transmittance_to_absorbance = function(dataset, percent = T){
  datamat = dataset$data
  if (!percent){
    absorbance.datamat = 2-log10(datamat*100)
  } else {
    absorbance.datamat = 2-log10(datamat)
  }
  dataset$data = absorbance.datamat
  dataset$labels$val = "A"
  dataset
}

#convert absorbance to transmittance
absorbance_to_transmittance = function(dataset){
  datamat = dataset$data
  transmittance.datamat = 100*(10^(-datamat))
  dataset$data = transmittance.datamat
  dataset$labels$val = "%T"
  dataset
}

# Smoothing - hyperSpec (spc.bin and spc.loess)
smoothing_interpolation = function(dataset, method = "bin", reducing.factor = 2, x.axis = NULL, p.order = 3, window = 11, deriv = 0){
  if (method == "bin") {
		dataset = smoothing_spcbin_hyperspec(dataset, reducing.factor)
	} 
  else if (method == "loess") {
		dataset = smoothing_spcloess_hyperspec(dataset, x.axis)
	} else if (method == "savitzky.golay"){
		dataset = savitzky_golay(dataset, p.order, window, deriv)
	}
	dataset
}

# spc.bin hyperSpec smoothing interpolation
smoothing_spcbin_hyperspec = function(dataset, reducing.factor = 2) {
	hyper.object = convert_to_hyperspec(dataset)
	smooth.result = hyperSpec::spc.bin(hyper.object, reducing.factor, na.rm = TRUE)
	res.dataset = convert_from_hyperspec(smooth.result)
  res.dataset$description = paste(dataset$description, "smoothed with hyperSpec spc.bin", sep="-")
  res.dataset$type = res.dataset$type
  res.dataset
}

# spc.loess hyperSpec smoothing interpolation 
smoothing_spcloess_hyperspec = function(dataset, x.axis = NULL){
  hyper.object = convert_to_hyperspec(dataset)
	if (is.null(x.axis)){
		smooth.result = hyperSpec::spc.loess(hyper.object, hyperSpec::wl(hyper.object), na.rm = TRUE)
	} else {
		smooth.result = hyperSpec::spc.loess(hyper.object, x.axis, na.rm = TRUE)
	}
  res.dataset = convert_from_hyperspec(smooth.result)
  res.dataset$description = paste(dataset$description, "smoothed with hyperSpec spc.loess", sep="-")
  res.dataset$type = res.dataset$type
  res.dataset
}

savitzky_golay = function(dataset, p.order, window, deriv = 0){
    if (window %%2 != 1 || window < 0) 
        stop("window size (window) must be a positive odd number")
    if (p.order >= window) 
        stop("window size (window) is too small for the polynomial order (p.order)")
    if (p.order < deriv) 
        stop("polynomial order p (p.order) should be geater or equal to differentiation order (deriv)")
    X = t(dataset$data)
    half_window = (window -1)/2
    b = outer(-half_window:half_window, 0:p.order, "^")
    A = MASS::ginv(b)
    result = matrix(data = 0, ncol=ncol(X),nrow=nrow(X))
    for (i in 1:nrow(X)){
        first.values = X[i,1] - abs( X[i,1:(half_window)] - X[i,1] )
        last.values = tail(X[i,], n = 1) + abs(X[i,(ncol(X)-half_window+1):ncol(X)] - tail(X[i,],n=1))
        all = c(first.values, X[i,], last.values)
        result[i,] = factorial(deriv) * convolve(all, A[deriv+1,], type="f")
    }
    colnames(result) = colnames(X)
    dataset$data = t(result)
    rownames(dataset$data) = colnames(X)
    colnames(dataset$data) = rownames(X)
    dataset$description = paste(dataset$description, "smoothed with savitzky-golay filter", sep = "-")
    dataset
} 


# DATA CORRECTION - functions to do spectra correction

"data_correction" = function(dataset, type = "background", method = "modpolyfit", ...){
	if (type == "background"){
		dataset = background_correction(dataset)
	} 
  else if (type == "offset"){
		dataset = offset_correction(dataset)
	} 
  else if (type == "baseline"){
		dataset = baseline_correction(dataset, method, ...)
	} 
	dataset
}

background_correction = function(dataset) {
  hyper.object = convert_to_hyperspec(dataset)
	background = hyperSpec::apply(hyper.object, 2, quantile, probs = 0.05)
	correction.result = hyperSpec::sweep(hyper.object, 2, background, "-")
  res.dataset = convert_from_hyperspec(correction.result)
  res.dataset$description = paste(dataset$description, "background correction", sep="; ")
  res.dataset$type = dataset$type
  res.dataset
}

offset_correction = function(dataset){
  hyper.object = convert_to_hyperspec(dataset)
	offsets = hyperSpec::apply(hyper.object, 1, min)
	correction.result = hyperSpec::sweep(hyper.object, 1, offsets, "-")
  res.dataset = convert_from_hyperspec(correction.result)
  res.dataset$description = paste(dataset$description, "offset correction", sep="; ")
  res.dataset$type = dataset$type
  res.dataset
}

# ... - extra parameters to baseline function
baseline_correction = function(dataset, method = "modpolyfit", ...){
	rnames = rownames(dataset$data)
	cnames = colnames(dataset$data)
	samples.df = t(dataset$data)
	bl = baseline::baseline(samples.df, method = method, ...)
	samples.df = baseline::getCorrected(bl)
  dataset$data = t(samples.df)
	rownames(dataset$data) = rnames
	colnames(dataset$data) = cnames
  dataset$description = paste(dataset$description, "baseline correction", sep="; ")
	dataset
}

# shifting spectra
# method - "constant" - uses a constant shift that is added to the x.values
# method - "interpolation" - uses interpolation - linear or spline according to "interp.function"
# shift.val - value of the shift (for constant and interpolation methods); can be a single value for all spectra
#			  "auto" - shifts are automatically determined
# or a vector of length = number of samples; can also be the string "auto" for automatic calculation of shifts
"shift_correction" = function(dataset, method = "constant", shift.val = 0, interp.function = "linear",
                              ref.limits = NULL) {
  
  x.vals = get_x_values_as_num(dataset)
  
  if(! length(shift.val) %in% c(1,num_samples(dataset)) ) {
    stop("Shift.val parameter has incorrect size: should be 1 or number of samples in the dataset")
  }
  else if (length(shift.val) == 1) {
	if (shift.val == "auto") {
        if (is.null(ref.limits) | length(ref.limits) != 2) {
			stop("Parameter ref.limits incorrect for automatic determination of shifts")
		}
        else { 
			shift.val = calculate_shifts(dataset, ref.limits)
		}
	}
  }
  if (method == "constant") {
    new.x.values = x.vals + shift.val
    dataset = set_x_values(dataset, new.x.values)
  }
  else if (method == "interpolation") {
    if (interp.function == "spline") {
      interp_fn = function(data, shift, x.values) {
        spline (x.values + shift, data, xout = x.values, method = "natural")$y
      }     
    }
    else if(interp.function == "linear") {
      interp_fn = function(data, shift, x.values) {
        approx(x.values + shift, data, xout = x.values, method = "linear")$y
      }
    }
    else stop("Interpolation function not defined")
    
    if (length(shift.val) == 1)
      newdata = apply(dataset$data, 2, interp_fn, shift = shift.val, x.values = x.vals)
    else {
      newdata = matrix(NA, nrow(dataset$data), ncol(dataset$data))
      for (i in 1:length(shift.val))
        newdata[,i] = interp_fn(dataset$data[,i], shift.val[i], x.vals)
    }
    rownames(newdata) = rownames(dataset$data)
    colnames(newdata) = colnames(dataset$data)
    dataset$data = newdata
  }
  else stop("Method is not defined")
  
  dataset
}

# calculate shifts based on a band of the spectra (see hyperSpec vignette sect. 12.2.1)
"calculate_shifts" = function(dataset, ref.limits = NULL)
{ 
  #x.vals = get.x.values.as.num(dataset)
  dataM = subset_x_values_by_interval (dataset, ref.limits[1], ref.limits[2])
  #xvalsM = x.vals[x.vals >= ref.limits[1] & x.vals <= ref.limits[2]]
  bandpos = apply (t(dataM$data), 1, find_max, get_x_values_as_num(dataM))
  refpos = find_max (colMeans(t(dataM$data)), get_x_values_as_num(dataM))
  refpos - bandpos
}

"find_max" = function (y, x){
  pos = which.max (y) + (-1:1)
  X = x [pos] - x [pos [2]]
  Y = y [pos] - y [pos [2]]
  X = cbind (1, X, X^2)
  coef = qr.solve (X, Y)
  - coef [2] / coef [3] / 2 + x [pos [2]]
}

# multiplicative scatter correction

"msc_correction" = function(dataset) {
  temp = t(dataset$data)
  newdata = pls::msc(temp)
  dataset$data = t(newdata)
  dataset
}

# first derivative

first_derivative = function(dataset) {
  new.data = apply(dataset$data, 2, diff)
  dataset$data = new.data
  dataset
}
