#calculates MODIFIED Short Time Fourier Transforms.
#if center = T, remove means
#if calc.null, calculate 'null hypothesis' by randomising data (sampling w/o replacement) and calculating FFT on that.
#if reassign, calculate reassigned stft


#MV method wrapper
#maybe we should implement sums?

stft <- function(X, start=0, end=1, length=NULL,  time.format = c("auto"), type = c("mv", "svm", "sum"), mv.indices, date.col, reassign = TRUE,plot.it = FALSE,...){
type = match.arg(type)
call <- match.call()
if (is.list(X)){
#if (!is.null(X$freq)) freq = X$freq

X = get.intervals(X, start, end, length, time.format, incl.date=TRUE, simplify  = TRUE)
if (missing( date.col)) date.col = TRUE
}

if (length(dim(X)) < 2) X = matrix(X, ncol = 1)

#is first col date-like?
if (missing(date.col)){

if ((ncol(X)> 1) && ( X[1,1] >= 365 * 60*60*24)){
print("Assuming first column is time.")
date.col = T
} else {
date.col = F
}
}
if (missing(mv.indices)) mv.indices <- 1:(min(ncol(X) - date.col, 3))

if (type == "svm"){
if (date.col){
obj1 = stftcalc(cbind(X[,1], sqrt(rowSums(X[,mv.indices + 1, drop = F]^2))), reassign = reassign,...)
} else {
obj1 = stftcalc( sqrt(rowSums(X[,mv.indices, drop = F]^2)),reassign = reassign, ...)
}
obj1$type = "svm"
} else {
ind = mv.indices[1]
if (date.col) ind = c(1,ind +1)

obj1 = stftcalc(X[,ind], reassign = reassign,...) 

if (length(mv.indices) > 1){

if (type == "mv"){
for (ind in mv.indices[-1] ) {

if (date.col) ind = c(1,ind +1)
obj = stftcalc(X[,ind], reassign = F, ...)
obj1$values = pmax(obj1$values, obj$values)
}
obj1$type = c("mv", paste(mv.indices, collapse=""))
} else {
for (ind in mv.indices[-1] ) {

if (date.col) ind = c(1,ind +1)
obj = stftcalc(X[,ind], reassign = F, ...)
obj1$values = sqrt(obj1$values^2 + obj$values^2)
}
obj1$type = c("sum", paste(mv.indices, collapse=""))

}
}


}
obj1$call = call
if (plot.it) plot(obj1)
obj1
}


stftcalc <- function(X, win=10, 
                 inc=  win/2, coef=Inf, 
		 wtype="hanning.window", freq , center = T, calc.null = F , pvalues = F, time = NULL, reassign = T , quiet = F)
  {
call = match.call()
if (length(dim(X)) ==2) {
if (is.null(time)) time = X[,1]

X = X[,2]

}

if ((length(time) >1) && (missing(freq))) freq =(length(time) -1) /( max(time) - min(time) )
if (missing(freq) ) freq = 1
inc0 = inc
win0 = win
coef0 = coef
win = round(win * freq)
inc = max(round(inc * freq),1)
coef = min(coef, floor(win/2))

Xdel = shift(X, c(1,0), fill = "edge")
    numcoef <- 2*coef
    if (win < numcoef)
      {
	win <- numcoef
if (!quiet)	cat ("stft: window size adjusted to", win, ".\n")
      }
    numwin <- trunc ((length(X) - win) / inc)

    ## compute the windows coefficients
    wincoef <- eval(parse(text=wtype))(win)

    ## create a matrix Z whose columns contain the windowed time-slices
pval = rep(0, numwin+1)
    z <- matrix (0, numwin + 1, win)
    y <- matrix(0, numwin+1, win)
    ydel <- matrix(0, numwin+1, win)
    st <- 1
if (!quiet) pb <- txtProgressBar(min = 0, max = 100,style=1)
    for (i in 0:numwin)
      {
	z[i+1, 1:win] <- (X[st:(st+win-1)] - mean(X[st:(st+win - 1)])* center) * wincoef
	y[i+1,] <- fft(z[i+1,] )
if (reassign){	
	z[i+1, 1:win] <- (Xdel[st:(st+win-1)] - mean(Xdel[st:(st+win - 1)])* center) * wincoef
	ydel[i+1,] <- fft(z[i+1,] )
}
if (pvalues){
temp = sample(X[st:(st + win - 1)])
temp = (temp - mean(temp) * center)*wincoef
temp = Mod(fft(temp))[1:coef]^2
pval[i+1] = wilcox.test( (Mod(y[ i+1,])^2 - mean(Mod(y[ i+1,])^2))^2, (temp - mean(temp))^2)$p.value
}
	st <- st + inc

if (!quiet) setTxtProgressBar(pb, 90* i/(numwin ))
      }
if (reassign){
yfreqdel = cbind(y[, win],y[, 2: win - 1])# t(apply(y, 1, function(t) shift(t, 1, fill = "loop")))

ydel = Arg(y*Conj(ydel)) *(freq/(2*pi))
yfreqdel = -(win/(2*pi*freq)) * Arg(y * Conj(yfreqdel)) + win/(2*freq)
} else {
yfreqdel = NULL
}



null.logmean = NULL
null.logsd = NULL
if (calc.null){
tmpdat = stft(sample(X), win = win0, 
                 inc= inc0, coef=coef0, 
		 wtype=wtype, freq = freq, center = T,  calc.null = F , quiet= T)
null.logmean = log(sqrt(mean((tmpdat$values)^2)))
null.logsd = log(sd(tmpdat$values))
}
if (!quiet)setTxtProgressBar(pb, 100)
if (is.null(time)){
    Y<- list (values = cbind(Mod(y[,1]) ,2*Mod(y[,(2):coef])), windowsize=win, increment=inc,
		  windowtype=wtype, center = center, sampling.frequency = freq, null.logmean = null.logmean, null.logsd = null.logsd, principals = (freq * (1:coef  - 1 ) / win)[apply( Mod(y[,(1):coef]),1, which.max)], frequency = (freq * (1:coef  - 1 ) / win), times =  (win/2 +  inc * 0:(nrow(y) - 1))/(freq), p.values = pval, LGD = yfreqdel[,1:coef], CIF = ydel[,1:coef] )
} else {
if (length(time) == 1){
times = (time  + (win/2 +  inc * 0:(nrow(y) - 1))/freq)
} else {
times = time[win/2 + inc* 0:(nrow(y) - 1) ] 
}

   Y<- list (values = cbind(Mod(y[,1]) ,2*Mod(y[,(2):coef])), windowsize=win, increment=inc,
		  windowtype=wtype, center = center, sampling.frequency = freq, null.logmean = null.logmean, null.logsd = null.logsd, principals = (freq * (1:coef  - 1 ) / win)[apply( Mod(y[,(1):coef]),1, which.max)], frequency = (freq * (1:coef  - 1 ) / win), times = times, p.values = pval, LGD = yfreqdel[,1:coef], CIF = ydel[,1:coef]  )
}
if (!quiet)close(pb)
Y$call = call
    class(Y) <- "stft"
    return(Y)
  }

hanning.window = function (n) {
    if (n == 1) 
        c <- 1
    else {
        n <- n - 1
        c <- 0.5 - 0.5 * cos(2 * pi * (0:n)/n)
    }
    return(c)
}
uniform.window = function(n){
rep(1, n)
}


#topthresh - threshold frequency at which to put higher frequency bins into a top plot # proportional for pval plot, else absolute?
#reassign - use reassigned stft?
plot.stft <- function (x,  mode = c("decibels", "modulus", "pval"), log = "", showmax = TRUE, median = FALSE, xaxis = TRUE, topthresh, reassign = (!(is.null(x$LGD)) && !("mv" %in% x$type)), ylim, xlim,new = TRUE, zlim.raw,zlim.quantile, cex,col = gray (63:0/63),...)
  {
    xv <- x$values
if (missing(cex)){
cex = 5
if (nrow(xv) > 500) cex = 2
}
if (missing(topthresh)){
topthresh = Inf
if (x$sampling.frequency > 30) topthresh = 15
}

if (median) xv = apply(xv,2, function(t) (runmed(t, k = 1 + 2 * min((length(t)-1)%/% 2, ceiling(0.005*length(t))) )))

mode = match.arg(mode)
if (mode == "decibels"){
xv = log(xv)
if (missing(zlim.raw)){

if (missing(zlim.quantile)){
zlim.raw = c(median(xv), Inf)
if (!is.null(x$null.logmean)) zlim.raw = c(x$null.logmean, Inf)

} else {
zlim.raw = c(quantile(xv, zlim.quantile[1]), quantile(xv, zlim.quantile[2]))
}
}
if (length(zlim.raw) == 1) zlim.raw = c(zlim.raw, zlim.raw + abs(zlim.raw)*0.0001)
xv = constrain(xv, zlim.raw[1], zlim.raw[2])

} else if (mode == "pval"){
xv = t(apply(xv, 1, function(t)  -log10(1-pexp(t^2, 1/mean(t^2)) )))
if (missing(zlim.raw)){
if (missing(zlim.quantile)){
zlim.raw = c(0, 15)
} else {
zlim.raw = c(quantile(xv, zlim.quantile[1]), quantile(xv, zlim.quantile[2]))
}
}
 
if (length(zlim.raw) == 1) zlim.raw = c(zlim.raw, zlim.raw + abs(zlim.raw)*0.0001)
xv = constrain(xv, zlim.raw[1], zlim.raw[2])

} else {
if (missing(zlim.raw)){
if (missing(zlim.quantile)){
zlim.raw = c(0, Inf)
} else {
zlim.raw = c(quantile(xv, zlim.quantile[1]), quantile(xv, zlim.quantile[2]))
}
}

if (length(zlim.raw) == 1) zlim.raw = c(zlim.raw, zlim.raw + abs(zlim.raw)*0.0001)
xv = constrain(xv, zlim.raw[1], zlim.raw[2])
}

timegrid = x$times
if (missing(ylim)) ylim = range( x$frequency)
if (missing(xlim)) xlim = range(timegrid)
frequency= x$frequency
if (topthresh < Inf){

if (new){
 layout(matrix(c(1,2,2,2), ncol = 1))
 par(mar = c(0,1,0,0))
 par(oma = c(5, 4, 4, 2) + 0.1)

binwidth = ceiling(length(timegrid) /100)
topind = 1:(floor(length(timegrid) / binwidth) * binwidth)

#res = epoch.apply(x$values, epoch.size = binwidth, function(t) apply(t,2,median)) #apply(x$values, 2, function(t)  apply(matrix(t, binwidth), 2, median))
res = apply(x$values[topind,], 2, function(t)  apply(matrix(t, binwidth), 2, median))
timegridtop = matrix(timegrid[topind], binwidth)[1,]
ind = ceiling(ncol(xv) * 1:20/20)
ylim =  c(0, quantile(sqrt(rowSums(res^2)), 0.95)*1.1)
if (xaxis){

plot(  convert.time(timegridtop), sqrt(rowSums(res^2)) , type="l", xaxt = "n", ylim =ylim, ...)
axis( 1, convert.time(timegridtop), labels = F)
timegridtop = convert.time(timegridtop)
} else {

plot(  (timegridtop), sqrt(rowSums(res^2))  , type="l", xaxt="n", ylim = ylim,...)
axis(1, pretty(timegridtop), labels = F)

}

for (k in ind){
if (frequency[k] <= (topthresh *2/3)){
colour = "black"
} else if (frequency[k] %bt% c(topthresh * 2/3, topthresh)){
colour = "red"
} else {
colour = "blue"
}

lines(timegridtop,sqrt(apply((res)[, k:ncol(res), drop=F]^2, 1, sum)), col=  colour)
abline(h = 0)
}



}
ylim[2] = min(ylim[2], topthresh)
#xv[, which(frequency > topthresh)]
xv = xv[,which(frequency <= topthresh)]
#do top thresholding
}

if(log == "y"){
frequency[1] = frequency[2]^2/frequency[3]
frequency = c(frequency, tail(frequency,1)^2/tail(frequency,2)[1])
ylim[1] = max(min(frequency), ylim[1])
}
if (reassign){
		frequency =  as.vector(x$CIF[,1:ncol(xv)])
if (is.numeric(col)){
	colours = col2rgb(palette()[col] )
	colours = rgb(colours[1], colours[2], colours[3], alpha = 255*(conv01(as.vector(xv)))  , maxColorValue = 255)
} else {
	colours = col[1+ (length(col) - 1) * ( conv01(as.vector(xv)))]
}
}
if (xaxis){
	if (reassign){
		time = convert.time(rep(x$times, ncol(xv) )+ as.vector(x$LGD[,1:ncol(xv)] ))
		if (new){
			plot( time, frequency,pch= ".", cex = cex , col = colours , log = log,ylim = ylim , xlim = convert.time(xlim), ...)
		}else{
			points ( time, frequency,pch= ".", cex = cex , col = colours , ylim = ylim ,  ...)
		}
#####
	}else {
	time = timegrid
	plot(convert.time(  seq(min(time), max(time), len = 20) ), rep(1,20), col=0, xlab = "time", ylab = "frequency", ylim = ylim, xlim = convert.time(xlim), xpd = NA, log =log)
#	par(new = T)
		image( x = convert.time(time) , y = frequency[1:ncol(xv) ],   z=xv, col=col, log = log, xaxt = "n", ylim = ylim,xlim = convert.time(xlim), add= T,...)
	}
} else {
	if (reassign){

		time = (rep(x$times, ncol(xv) )+ as.vector(x$LGD[,1:ncol(xv)] ))
		if (new){
			plot( time, frequency,pch= ".", cex = cex , col = colours , log = log,ylim = ylim , xlim = xlim, ...)
		} else {
			points( time, frequency,pch= ".", cex = cex , col = colours , ylim = ylim , ...)
		}
	#######
	} else {
	time = timegrid
	    image( x = time , y = frequency[1:ncol(xv)],   z=xv, col=col, log = log, ylim = ylim,xlim=xlim,...)
	}
}
axis(2, pretty(constrain(frequency, ylim[1], ylim[2]), min(floor(topthresh), floor(max(frequency))) ), labels = F, tcl = -0.2)
if (as.numeric(showmax) > 0){
#points ( time, x$principals, col=2 * (rowMeans(xv) > 1 * x$null.logmean)  , pch=".", cex = 3)

pseudolines(timegrid, x$principals, col = 2, pch = ".", lwd = 2, cex = 2)

}
if (as.numeric(showmax) > 1){
pseudolines(timegrid, frequency[ apply(x$values, 1, function(t) which.max(replace(t, which.max(t), -Inf)))], col=3, pch = ".", lwd = 2, cex = 2)
}

}

#example code
#plot(stft(subs(mag, 0.94,0.96), win = 1024, plot = F, coef = 512), zlog = T, log="y")
#plot(stft(subs(mag, 0.7,8), win = 1024, plot = F, coef = 512), zlog = T, log="y")
#plot(stft(subs(mag, 0.0001,0.005), win = 1024, plot = F, coef = 512), zlog = T)



#plot(stft(subs(mag, 0.7,0.8), win = 1024, plot = F), zlog = T, log = "y")
#plot( stft(rep(1, 1000) + c(sin(1:500/ 10 * 2*pi), rep(0, 500)) + c(rep(0, 300),sin(1:500/ 20 * 2*pi), rep(0, 200)) , freq = 1, plot.it = F), log="x")
#stft(sin(1:1000 / (1 +sqrt(1000:1)) * 2 * pi), freq = 1)
# stft(rep(1, 1000) + sin(1:1000/ 10 * 2*pi), freq = 1)

#subsets a proportion of the dataset, or a certain length of the dataset starting at a specific proportion position
subs <- function(x, a,b){
len = length(x)
if (a > 1){
return(x[a : (a+b - 1)])
} else if (b > 1) {
return( x [ floor(a * len): (floor(a*len) + b - 1)])
} else {
return (x[floor(a*len) : floor(b*len)])
}
}

#plot(stft(subs(mag, 0.7,0.8), win = 1024, plot = F, coef = 512), zlog = T, log="y")

#gets fft components corresponding to frequencies
getfreqs = function(x, frequencies){
n = length(x)
fftobj = fft(x)
frequencies = c(frequencies, n - frequencies) + 1
frequencies = frequencies[which(frequencies != n+1)]
fftobj = replace(fftobj, (1:n)[ - frequencies], 0)

return(Re(fft(fftobj, inverse=T))/n)
}


print.stft = function(x,...){
cat("STFT object:\n")
cat(format.GRtime(x$times[1],  format = "%y-%m-%d %H:%M:%OS3 (%a)")," to ", format.GRtime((tail(x$times,1)), format = "%y-%m-%d %H:%M:%OS3 (%a)"), "\n")
cat(nrow(x$values), "increments of" , round(x$increment/x$sampling.freq, 3), "s \n")
cat("Window size: " , x$windowsize, "(", round(x$windowsize/x$sampling.frequency, 3), "s ) -> f resolution: ", round(x$frequency[2],3), "Hz\n")
if ("svm" %in% x$type) cat("[SVM]")
if ("mv" %in% x$type) cat("[MV-", x$type[2], "]")
if ("sum" %in% x$type) cat("[SUM-", x$type[2], "]")
if (!is.null(x$LGD)) cat ("[Reassign]")
cat("\n------ \n")
cat("{" ,format(x$call), "}\n")
}

