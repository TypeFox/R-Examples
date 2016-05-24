fac2list<-function(x,g){
#
# data are stored in x
# information about the level of the value in x is stored in g,
# which can be a matrix with up to 4 columns
#
# sort the data in x into groups based on values in g.
# store results in list mode.
#
#  Example: fac2list(m[,2],m[,4]) would sort the values
#  in column 2 of m according to the values in column 4 of m
#
g=as.data.frame(g)
L=ncol(g)
g=listm(g)
for(j in 1:L)g[[j]]=as.factor(g[[j]])
g=matl(g)
Lp1=L+1
if(L>4)stop("Can have at most 4 factors")
if(L==1){
res=selby(cbind(x,g),2,1)
group.id=res$grpn
res=res$x
}
if(L>1){
res=selby2(cbind(x,g),c(2:Lp1),1)
group.id=res$grpn
res=res$x
}
print("Group Levels:")
print(group.id)
res
}
