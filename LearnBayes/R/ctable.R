ctable=function(y,a)
#
# C_TABLE Bayes factor for testing independence in a contingency table.
#	BF=C_TABLE(Y,A) returns the Bayes factor BF against independence in a 
#	2-way contingency table using uniform priors, where Y is a matrix
#	containing the 2-way table, and A is a matrix of prior parameters
#------------------------
# Written by Jim Albert
# albert@bgnet.bgsu.edu
# November 2004
#------------------------
{
ldirich=function(a)
{
val=sum(lgamma(a))-lgamma(sum(a))
return(val)
}
ac=colSums(a); ar=rowSums(a)
yc=colSums(y); yr=rowSums(y)

d=dim(y); oc=1+0*yc; or=1+0*yr; I=d[1];J=d[2]

lbf=ldirich(c(y)+c(a))+ldirich(ar-(J-1)*or)+ldirich(ac-(I-1)*oc)-
    ldirich(c(a))-ldirich(yr+ar-(J-1)*or)-ldirich(yc+ac-(I-1)*oc)

bf=exp(lbf)
return(bf)

}
