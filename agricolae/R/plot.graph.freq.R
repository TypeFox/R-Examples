`plot.graph.freq` <-
function (x, breaks=NULL,counts=NULL,frequency=1, plot=TRUE, nclass=NULL,xlab="",ylab="",axes = "",las=1,...)
{
if(axes=="") ejes=TRUE
else {
ejes<-FALSE
if (axes) ejes<-TRUE
}
if (xlab=="") xlab= deparse(substitute(x))
if (is.numeric(x) & is.null(counts)) {
x<-na.omit(x)
	# histogram
if (is.null(nclass)) {
if (is.null(breaks)) {
breaks <- sturges.freq(x)$breaks
}
}
else {
breaks <- sturges.freq(x,k=nclass)$breaks
}

k<-length(breaks)
n<- length(x)
counts <- rep(0,k-1)
for (i in 1:n) {
for (j in 1:(k-2)) {
if( (x[i] >= breaks[j]) && (x[i] < breaks[j + 1])) counts[j]<-counts[j]+1
}
}
for (i in 1:n) {
	if( (x[i] >= breaks[k-1]) && (x[i] <= breaks[k])) counts[k-1]<-counts[k-1]+1
}
    k <- length(counts)
    mids <- rep(0, k)
    ancho <- rep(0, k)
    for (i in 1:k) {
        mids[i] <- (breaks[i] + breaks[i + 1])/2
        ancho[i] <- (breaks[i + 1] - breaks[i])
    }
    altura <- round(1.1 * max(counts), 0)
}
#############
else  {
if( is.list(x)) {
breaks<- x$breaks
counts <- x$counts
}
else breaks <- x
k<-length(counts)
mids<-rep(0,k)
ancho<-rep(0,k)
for (i in 1:k) {
mids[i]<-(breaks[i]+breaks[i+1])/2
ancho[i]<-(breaks[i+1]-breaks[i])
}
}
################
a<-breaks[1]-ancho[1]/2
b<-breaks[k+1]+ancho[k]/2
relative<-round(counts/sum(counts),4)
density <- relative/ancho
histogram<-structure(list(breaks=breaks,counts=counts,mids=mids,relative=relative,density=density),class="graph.freq")

if(plot) {
x <- c(a, b)
if(frequency==1)height<-round(1.1*max(counts),1)
if(frequency==2)height<-round(1.1*max(relative),4)
if(frequency==3)height<-round(1.1*max(density),4)
y <- c(0, height)
#suppressWarnings(warning(plot(x,y, type = "n", xlab=xlab,ylab=ylab,...)))
if(ejes){
suppressWarnings(warning(plot(x,y, type = "n", xlab=xlab,ylab=ylab,axes=FALSE,...)))
axis(1,breaks,las=las)->ax; axis(2,las=las)->ay
}
else suppressWarnings(warning(plot(x,y, type = "n", xlab=xlab,ylab=ylab,axes=axes,...)))
if (frequency==1) {
for (j in 1:k) {
suppressWarnings(warning(rect(breaks[j], 0, breaks[j + 1], counts[j], ...)))
}}
if (frequency==2) {
for (j in 1:k) {
suppressWarnings(warning(rect(breaks[j], 0, breaks[j + 1], relative[j], ...)))
}}
if (frequency==3) {
for (j in 1:k) {
suppressWarnings(warning(rect(breaks[j], 0, breaks[j + 1], density[j], ...)))
}}
}
    invisible(histogram)
}
