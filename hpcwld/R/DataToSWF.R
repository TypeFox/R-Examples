DataToSWF <-
function(Frame,filename="output.swf") {
if(is.null(Frame) || is.null(filename))
	stop("'Frame' must be defined!")
if (is.null(Frame$T) || is.null(Frame$S) || is.null(Frame$N) || is.null(Frame$D)) 
	stop("'T','S','N','D' in data frame must be defined!")
ToSWF(T=Frame$T,S=Frame$S,N=Frame$N,D=Frame$D,filename)
}
