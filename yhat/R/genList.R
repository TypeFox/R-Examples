genList <-
function(ivlist, value){
########################################################################################

numlist = length(ivlist)
newlist<-ivlist
newlist=0
for (i in 1:numlist){
	newlist[i]=abs(ivlist[i])+abs(value)
	if (((ivlist[i]<0) && (value >= 0))|| ((ivlist[i]>=0) && (value <0)))newlist[i]=newlist[i]*-1
}
return(newlist)
}

