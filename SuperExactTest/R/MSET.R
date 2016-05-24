#Statistical Test of Multi-Set Intersections
#Author: Minghui Wang, minghui.wang@mssm.edu
#Date: 20 July, 2014

#a wrapper to call cpsets
MSET <- function(x,n,lower.tail=TRUE,log.p=FALSE){
#Input:
#x           a list of sets.
#n           background size.
#lower.tail  logical; if TRUE (default), probability is P[overlap < m], otherwise, P[overlap >= m], where m is the number of elements shared by all sets.
#log.p       logical; if TRUE, probabilities p are given as log(p).
	L=sapply(x,length)
	if(n<1 | any(L>n) | any(L==0)) stop('Invalid input\n')
	intersects=intersect.list(x)
	Obs=length(intersects)
	Exp=n*do.call(prod,as.list(L/n))
	p=cpsets(max(Obs-1,0,na.rm=TRUE),L,n,lower.tail,log.p)
	list(intersects=intersects, FE=Obs/Exp, p.value=p)
}
#Sample usage:
##not run###
#sample data
#x=list(S1=c(letters[1:20]),S2=c(letters[10:26]),S3=c(sample(letters[1:26],10,FALSE)),S4=c(sample(letters[1:26],10,FALSE)))
#MSET(x,26,FALSE)
