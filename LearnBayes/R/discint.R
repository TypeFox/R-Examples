"discint" <-
function(dist,prob)
#
# DISC_INT Computes a highest probability interval for a discrete distribution.  
#	LIST=DISCINT(DIST,PROB) gives a list, where LIST.set is the set of values and
# 	LIST.prob is the exact probability context EPROB, where DIST=[VALUE,PROBABILITY]
#	is the matrix which contains the discrete distribution and PROB
#	is the probability content desired.
#------------------------
# Written by Jim Albert
# albert@bgnet.bgsu.edu
# November 2004
#------------------------
{	
x=dist[,1]; p=dist[,2]; n=length(x)

sp=sort(p,index.return=TRUE)
ps=sp$x
i=sp$ix[seq(n,1,-1)]; ps=p[i]; xs=x[i]
cp=cumsum(ps)
ii=1:n
j=ii[cp>=prob]; j=j[1]
eprob=cp[j]; set=sort(xs[1:j])
v=list(prob=eprob,set=set)
return(v)
}

