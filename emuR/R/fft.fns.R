# 
# 
# 
# "fcalc"<- function(fftdata, byrow = T, samfreq = 20000, nyq = samfreq/2,
#                    low = 0, high = nyq, fun = sum, ...)
# {
#   ## fftdata: a matrix of fft values as returned by muspec() and muslice()
#   ## low: starting frequency over which the function in fun is to be applied
#   ## high: ending frequency in Hz over which the function is to eb applied
#   ## fun: apply which function? Can be any S-PLUS arithmetic function
#   ## e.g. fft.calc(fftdata, low=1000, high=4000, fun=var) calculates
#   ## the variance of each spectrum in the 1-4 kHz range
#   ## written by Jonathan Harrington, 1992	
#   mat <- NULL
# 
#   for(j in 1:length(low)) {
#     values <- fft.extract(fftdata, samfreq, nyq, low[j], high[j])
#     newdata <- 10^(values/10)
# 
#     if(!is.matrix(newdata))
#       resdata <- fun(newdata, ...)
#     else if(byrow)
#       resdata <- apply(newdata, 1, fun, ...)
#     else
#       resdata <- apply(newdata, 2, fun, ...)
#     
#     mat$vals <- cbind(mat$vals, resdata)
#     olab <- paste(paste(low[j], high[j], sep = "-"), "Hz", sep = "")
#     mat$dimlabs <- c(mat$dimlabs, olab)
#   }
#   
#   mat$vals <- 10 * log(mat$vals, base = 10)
#   
#   if(ncol(mat$vals) > 1)
#     dimnames(mat$vals) <- list(NULL, mat$dimlabs)
#   else
#     mat$vals <- c(mat$vals)
#   
#   mat$vals
# }
# 
# 
# "fft.complex"<- function(data, normlen = T)
# {
#   ## data: usally a  matrix of sampled data values with successive columns
#   ## corresponding to successive segments
#   ## data can also be a vector corresponding to a single segment
#   ## calculates the real and imaginary parts of the fft for each column 
#   ## if normlen is True, the fft is normalised by dividing each
#   ## value by the square root of its length
#   if(is.matrix(data) == F) {
#     n <- length(data)
#     fftres <- fft(data)
#   }
#   else {
#     fftres <- apply(data, 2, fft)
#     n <- nrow(data)
#   }
#   if(normlen)
#     fftres <- fftres/sqrt(n)
#   fftres
# }
# 
# 
# "fft.extract"<- function(fftdata, samfreq = 20000, nyq = samfreq/2,
#                          low = 0, high = nyq)
# {
#   ## assume a sampling frequency of n, i.e. 20 kHz
#   ## extract fft components in the frequency range low to high, where low
#   ## and high are in Hertz
#   
#   if(!is.matrix(fftdata)) fftdata <- rbind(fftdata)
#   fftlen <- ncol(fftdata)
# 
#   if(low != 0)
#     left <- round(((low/nyq) * (fftlen - 1)) + 1)
#   else 
#     left <- 1
# 
#   if(high != nyq)
#     right <- round(((high/nyq) * (fftlen - 1)) + 1)
#   else 
#     right <- fftlen
# 
#   fftdata[, left:right]
# }
# 
# 
# "fft.mod"<- function(fftvals, unreflected = F, transpose = T)
# {
#   ## calculate the Mod values from fftvals. fftvals consists
#   ## of real and imag. parts and is usually the output of
#   ## Fft.complex. Successive columns of fftvals relate
#   ## to successive segments. fftvals may also be a vector
#   ## that relates to a single segment.
#   ## If unreflected is true, only the first n/2
#   ## points are returned per column. If transpose is true
#   ## the resulting matrix is transposed such that the nth row
#   ## that is derived corresponds to the nth row of the segment
#   ## list that was used to produce fftvals
#   if(is.matrix(fftvals) == F) fftvals <- cbind(fftvals)
#   n <- nrow(fftvals)
#   fft.mod <- apply(fftvals, 2, Mod)
# 
#   if(unreflected) fft.mod <- fft.mod[1:(n/2),  ]
# 
#   if(transpose) fft.mod <- t(fft.mod)
#   
#   fft.mod
# }
# 
# 
# "fft.ratio"<- function(fftdata, samfreq = 20000, nyq = samfreq/2,
#                        low.a = 0, high.a = nyq, 
#                        low.b = 0, high.b = nyq)
# {
#   ## finds energy ratios
#   ## fftdata: a matrix of spectral values, as returned by muspec, muslice
#   ## low.ahigh.a, low.b, high.b: values in Hertz
#   ## calculates the ratio of energy in the range low.a to high.a : low.b to high.b
#   ## e.g. to calculate the ratio of energy from 0-300 Hz to 2000-4000 Hz, 
#   ## low.a=0, low.b=300, high.a=2000, high.b=4000. The default for low.a, low.b
#   ## is 0 Hz, for high.a, high.b 10000 Hz. Therefore to calculat the ratio
#   ## of energy from 2-4 kHz to the rest of the spectrum, enter low.a=2000, 
#   ## high.a=4000 leaving low.b, high.b unspecified
#   ## written by Jonathan Harrington, 1992
#   values.a <- fft.extract(fftdata, samfreq, nyq, low.a, high.a)
#   newdata <- 10^(values.a/20)
#   if(!is.matrix(newdata))
#     resdata.a <- sum(newdata)
#   else resdata.a <- apply(newdata, 1, sum)
#   values.b <- fft.extract(fftdata, samfreq, nyq, low.b, high.b)
#   newdata.b <- 10^(values.b/20)
#   if(!is.matrix(newdata.b))
#     resdata.b <- sum(newdata.b)
#   else resdata.b <- apply(newdata.b, 1, sum)
#   resdata.a/resdata.b
# }
# 
# "fplot"<- function(fftdata, labs = NULL, which = NULL, colour = T,
#                    linetype = F, samfreq = 20000, nyq = samfreq/2,
#                    xlab = "Frequency (Hz)", ylab = "Intensity (dB)",
#                    low = 0, high = nyq, dbrange = NULL, axes = T,
#                    main = "", average = F, smoothing = F, points = 20,
#                    coeff = F, type = "l", super = F, legn="tl", cex=1)
# {
#   ## a matrix of fft values, as returned from muspec() or muslice()
#   ## low: plot from this frequency in Hz
#   ## high: plot up to this frequency in Hz (default range is 0-10000)
#   ## dbrange: specify a range for the y-axis in db
#   ## axes: if F, no axes will be plotted
#   ## main: provide a main axis title
#   ## super: superimpose FFTs that occur in successive rows of fftdata
#   ## on the same plot
#   ## labs: provide a label file which is parallel to fftdata;
#   ## the resulting plot will be color-coded (use super=T to superimpose the ffts
#   ## colmain: specify a color for the axes
#   ## returns: spectral plots, db against Hz
#   ## examples: mu.sub37
#   ## written by Jonathan Harrington, 1992
#   ## assume a sampling frequency of n, i.e. 20 kHz
#   flag <- F
#   colour.flag <- colour
#   linetype.flag <- linetype
# 
#   if(!is.logical(colour))   stop("colour must be T or F")
#   if(!is.logical(linetype)) stop("linetype must be T or F")
# 
#   mat <- NULL
#   if(!is.matrix(fftdata))
#     fftdata <- rbind(fftdata)
#   if(is.null(labs))
#     {legn <- F
#      flag <- T
#      labs <- rep(".", nrow(fftdata))}
# 
#   if(!is.null(which)) {
#     temp <- muclass(labs, which)
#     fftdata <- fftdata[temp,  ]
#     labs <- labs[temp]
#   }
#   if(average)
#     {
#       av.data <- fplot.mean(fftdata, labs)
#       fftdata <- av.data$mspec
#       labs <- av.data$lab
#     }
#   if(smoothing) {
#     if(super) {
#       dbrange <- range(fftdata)
#       zdat <- fftdata
#     }
#     cdat <- cepstrum(fftdata, points = points, spectrum = T)
#     fftdata <- -cdat$cep
#     if(!is.matrix(fftdata))
#       fftdata <- rbind(fftdata)
#   }
# 
#   col.lty <- mu.colour(labs, colour, linetype)
#   colour <- col.lty$colour
#   lty <- col.lty$linetype
# 
#   fftlen <- ncol(fftdata)
#   if(is.null(dbrange)) dbrange <- range(fftdata)
#   if(low != 0)
#     left <- round(((low/nyq) * (fftlen - 1)) + 1)
#   else 
#     left <- 1
# 
#   if(high != nyq)
#     right <- round(((high/nyq) * (fftlen - 1)) + 1)
#   else 
#     right <- fftlen
# 
#   nums <- seq(low, high, length = (right - left) + 1)
# 
#   for(i in 1:nrow(fftdata)) {
#     plot(nums, fftdata[i, left:right], type = type, xlab = "",
# 	 ylab = "", col = colour[i], ylim = dbrange, 
# 	 , lty = lty[i])
#     par(new = T)
#   }
# 
#   if(axes) 
#     title(main = main, xlab = xlab, ylab = ylab, col = 1, cex=cex)
# 
#   if(legn != F){
#     legn <- mu.legend(legn, c(low, high), dbrange)
#     legend(legn$x, legn$y, col.lty$legend$lab, col = col.lty$legend$col, 
# 	   lty = col.lty$legend$lty, cex = cex)
#   }
#   if(smoothing) {
#     if(super) {
#       par(new = T)
#       if(flag) labs <- NULL
#       
#       fplot(zdat, labs = labs, samfreq = samfreq, nyq = nyq, 
# 	    xlab = "", ylab = "", low = low, high = high, 
# 	    dbrange = dbrange, main = "", cex=cex, 
# 	    colour=colour.flag, linetype=linetype.flag)
#       par(new = F)
#     }
#     if(coeff)
#       mat$coeff <- cdat$coeff
#   }
#   mat
# }
# 
# 
# "fplot.mean"<- function(fftdata, labs)
# {
#   ## fftdata: spectral values
#   ## labs: a parallel label file
#   ## returns: $mspec an averaged spectrum, per label-type in labs
#   ## returns: $lab: the label corresponding to the row of mspec
#   mat <- NULL
#   fftdata <- 10^(fftdata/20)
# 
#   for(j in unique(labs)) {
#     temp <- labs == j
#     vals <- fftdata[temp,  ]
#     mvals <- apply(vals, 2, mean)
#     mat$mspec <- rbind(mat$mspec, mvals)
#     mat$lab <- c(mat$lab, j)
#   }
# 
#   mat$mspec <- 20 * log(mat$mspec, base=10)
#   mat
# }
# 
# 
# "moment"<- function(specvals, least = T, nyq = 10000, low = 0, high = nyq)
# {
#   ## specvals: the output of muspec; dB-FFT values
#   ## least: this normalises each spectrum so that its minimum
#   ## db value is set to 0
#   ## low, high: (in Hz). Over which spectral range do you
#   ## want to calculate the first two spectral moments
#   ## returns: $first: the first spectral moment (spectral centre
#   ## of gravity;
#   ## $second: the second spectral moment (variance, or moment
#   ## of inertia
#   matout <- NULL
#   fftlen <- ncol(specvals)
#   if(!is.matrix(specvals))
#     specvals <- rbind(specvals)
#   if(least) {
#     minv <- apply(specvals, 1, min)
#     minv <- rep(minv, times = rep(fftlen, length(minv)))
#     minv <- matrix(t(minv), ncol = fftlen, byrow = T)
#     specvals <- specvals - minv
#   }
#   x <- seq(0, nyq, length = fftlen)
#   x <- rep(x, nrow(specvals))
#   x <- matrix(t(x), ncol = fftlen, byrow = T)
#   if((low != 0) | (high != nyq)) {
#     x <- fft.extract(x, low = low, high = high)
#     specvals <- fft.extract(specvals, low = low, high = high)
#   }
#   prodvals <- x * specvals
#   sumvals <- apply(prodvals, 1, sum)
#   sumspecvals <- apply(specvals, 1, sum)
#   matout$first <- sumvals/sumspecvals
#   sqprodvals <- specvals * (x^2)
#   sumsqprodvals <- apply(sqprodvals, 1, sum)
#   matout$second <- sqrt((sumsqprodvals/sumspecvals) - (matout$first^2))
#   matout
# }
# 
# 
# 
# "moments" <-
# function(count, x, minval=F)
# 
# {
# # compute moments. x is a numeric class
# # count is the frequency with which that 
# # particular class occurs
# # This function gives exactly the same
# # results as those for the mean, variance
# # skewness and kurtosis in example Table 3.13.1
# # p. 87, Snedecor & Cochran, 'Statistical Methods'
# # 6th Edition, 1975. Let the arguments count and x
# # equal f and U respectively in their example
# # the centre of gravity with minval = F.
# # the first two moments in this function
# # also give the same results as in Harrington & Cassidy.
# if(minval)
# count <- count - min(count)
# if(missing(x))
# {
# if(is.spectral(count))
# x <- trackfreq(count)
# else
# x <- 0:(length(count)-1)
# }
# k <- 1
# mom1 <- sum((x - 0)^k * count) / sum(count)
# # the variance
# k <- 2
# mom2 <- sum((x - mom1)^k * count) / sum(count)
# 
# # third moment
# k <- 3
# mom3 <- (sum((x - mom1)^k * count) / sum(count)) / (mom2 * sqrt(mom2))
# 
# # fourth moment
# k <- 4
# # peaked distributions show positive kurtosis
# # flat-topped distributions show negative kurtosis
# mom4 <- (sum((x - mom1)^k * count) / sum(count)) / mom2^2 - 3
# c(mom1, mom2, mom3, mom4)
# }
# 
# 
# 
# 
# # Local Variables:
# # mode:S
# # S-temp-buffer-p:t
# # End:
