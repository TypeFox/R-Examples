bootstrapping <-
function(pos, neg)
{
rval.pos=NULL
rval.neg=NULL
for (i in 1:length(pos)){
pp=pos[i]
nn=neg[i]
n.pos=sum(sample(c(rep(T, pp), rep(F, nn)), pp+nn, replace=T))
n.neg=pp+nn-n.pos
rval.pos=c(rval.pos, n.pos)
rval.neg=c(rval.neg, n.neg)
}
return (list(rval.pos, rval.neg))
}
