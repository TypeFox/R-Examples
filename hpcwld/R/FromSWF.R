FromSWF <-
function(filename) {
if(is.null(filename))
	stop("The variable 'filename' must be defined!")
if(!is.character(filename) || length(filename)==0)
	stop("Incorrect 'filename'!")
	
inp=scan(filename,list(0, t=0, D=0, S=0, N=0, 0,0,0,0,0,0,0,0,0,0,0,0,0))
n=length(inp$t)
T=inp$t[2:n]-inp$t[1:(n-1)]
S=inp$S[1:(n-1)]
D=inp$D[1:(n-1)]
N=inp$N[1:(n-1)]
res=data.frame(T,S,N,D)
res
}
