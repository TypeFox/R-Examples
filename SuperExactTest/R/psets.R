#Compute multi set intersection test
#Author: Minghui Wang, minghui.wang@mssm.edu
#Date: 20 July, 2014

#Input:
#x           number of overlap elements between the subsets
#L           a vector of subset sizes
#n           background size
#log.p       logical; if TRUE, probability p is given as log(p).
#lower.tail  logical; if TRUE (default), probability is P[overlap <= x], otherwise, P[overlap > x].

#distribution function
cpsets <- function(x,L,n,lower.tail=TRUE,log.p=FALSE,simulation.p.value=FALSE,number.simulations=1000000){
	if(length(L)<2) stop('L should have at least two entries\n')
	if(n<1 | any(L>n) | any(L<x)){
		warning('Background population size is too small.')
		return(0)
	}
	if(simulation.p.value) return(cpsets.simulation(x,L,n,lower.tail,log.p,number.simulations))
	L=sort(L)
	.C("C_pmvhyper",as.integer(x),length(L),as.integer(L),as.integer(n),as.numeric(0.0),as.integer(lower.tail),as.integer(log.p))[[5]]
}
cpsets.simulation <- function(x,L,n,lower.tail=TRUE,log.p=FALSE,number.simulations=1000000){
	nL=length(L)
	cc=sapply(1:number.simulations,function(i){
		Ls=sapply(L,function(m) sample.int(n=n,size=m,replace=FALSE))
		sum(table(unlist(Ls))==nL) <= x
	})
	p=sum(cc)/number.simulations
	if(lower.tail==F) p=1-p
	ifelse(log.p,log(p),p)
}
#density function
dpsets <- function(x,L,n,log.p=FALSE){
	if(any(L>n) | any(L<x)) return(0)
	nL=length(L)
	if(nL<2) stop('L should have at least two entries\n')
	L=sort(L)
	.C("C_dmvhyper",as.integer(x),length(L),as.integer(L),as.integer(n),as.numeric(0.0),as.integer(log.p))[[5]]
}
#Sample usage:
#fake data
#n=500; A=260; B=320; C=430; D=300; x=170; L=c(A,B,C,D)
#(d=dpsets(x,L,n))
#(p=cpsets(x,L,n,lower.tail=F))
