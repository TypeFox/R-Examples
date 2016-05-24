`twoDirections` <-
function(x,w) {
zinc<-pava(x,w); zdec<--pava(-x,w)
sinc<-sum(w*(zinc-x)^2); sdec<-sum(w*(zdec-x)^2)
if (sinc < sdec) return(zinc) else return(zdec)
}

