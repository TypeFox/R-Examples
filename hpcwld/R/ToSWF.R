ToSWF <-
function(T, S, N, D, filename="output.swf") {
if(is.null(T) || is.null(S) || is.null(N) || is.null(D) || is.null(filename))
	stop("'T','S','N', and 'D' must be defined!")
Lst=list(T,S,N,D)
for(i in 1:4)
	if (!is.numeric(Lst[[i]]) || length(Lst[[i]]) == 0) 
		stop("Only numeric vectors allowed!")
n=length(T)
if(length(S)<n || length(N)<n || length(D)<n)
	stop("Length of 'T' vector should be not more than others!")

F1=1:n
F2=vector("numeric",length=n)
for(i in 2:n) F2[i]=F2[i-1]+T[i-1]
F3=D[1:n]
F4=F9=S[1:n]
F5=F8=N[1:n]
F6=F7=F10=F11=F14=F15=F16=F17=F18=rep(-1,times=n)
F12=F13=rep(1,times=n)
write.table(format(data.frame(F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15,F16,F17,F18)),file=filename,quote=FALSE,col.names=FALSE,row.names=FALSE)
}
