# P(B=b | W=w) = dhyper(b,3,5,3-w)
require(MASS)              # for fractions()
outer(0:2,0:3,function(x,y){dhyper(y,3,5,3-x)}) -> probs
colnames(probs) = paste("B=",0:3,sep="")
rownames(probs) = paste("W=",0:2,sep="")
print(probs)
print(fractions(probs))
#
# P(R=r | W=w) = dhyper(r,5,3,3-w)
outer(0:2,0:3,function(x,y){dhyper(y,5,3,3-x)}) -> probs
colnames(probs) = paste("R=",0:3,sep="")
rownames(probs) = paste("W=",0:2,sep="")
print(probs)
print(fractions(probs))
