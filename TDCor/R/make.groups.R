make.groups <-
function(x1,x2,x3,x4,delaymax,np)

{



k=sort(unique(c(x1,x2,x3,x4)))

# Calculation of maximum number of groups and group limits


if (length(k)>1)
{
	u=abs(outer(k,k,"-"))
	b=shortest.path(1-sign(floor(u/delaymax)),1)$BS
	limgr=sign(apply(abs(as.matrix(b[,1:(length(k)-1)])-as.matrix(b[,2:length(k)])),2,sum))
	MaxNb=sum(limgr)+1
	GrpLim=c(1,seq(2,length(k))[limgr==1],length(k)+1)
}else
{
	MaxNb=1
	GrpLim=c(1,length(k)+1)
}

# Creating preliminary groups without repeats

groups=list()
length(groups)=MaxNb

for (i in 1:MaxNb)
{
	groups[[i]]=k[GrpLim[i]:(GrpLim[i+1]-1)]
}


# Breaking down groups in which one (or more) method(s) contribute(s) several times [which in practice should most certainly never happen...]

# and recomposing the repeats

for (i in 1:MaxNb)

{


grp=groups[[i]]

maxpart=max(length(na.omit(match(x1,grp))),length(na.omit(match(x2,grp))),length(na.omit(match(x3,grp))),length(na.omit(match(x4,grp))))



if (maxpart>1)

{

l1=grp[na.omit(match(x1,grp))]

l2=grp[na.omit(match(x2,grp))]

l3=grp[na.omit(match(x3,grp))]

l4=grp[na.omit(match(x4,grp))]

FirstSubGrp=TRUE


for (i1 in 1:max(1,length(l1)))

{for (i2 in 1:max(1,length(l2)))

{for (i3 in 1:max(1,length(l3)))

{for (i4 in 1:max(1,length(l4)))

{

ng=na.omit(c(l1[i1],l2[i2],l3[i3],l4[i4]))

if (max(ng)-min(ng)>=delaymax & length(ng)==4)
{
	ng=ng[substract(1:4,which.max(abs(ng-mean(ng))))]
}

if (max(ng)-min(ng)<=delaymax & length(ng)>=3)
{
	if (FirstSubGrp)
	{	groups[[i]]=ng
		FirstSubGrp=FALSE
	}else

	{

	groups[[length(groups)+1]]=ng

	}
}

}

}

}

if (FirstSubGrp)
{
	groups[[i]]=numeric()
}


}

}else

{


ng=c(grp[na.omit(match(x1,grp))],grp[na.omit(match(x2,grp))],grp[na.omit(match(x3,grp))],grp[na.omit(match(x4,grp))])

if (max(ng)-min(ng)>=delaymax & length(ng)==4)
{
	ng=ng[substract(1:4,which.max(abs(ng-mean(ng))))]
}

if (max(ng)-min(ng)<=delaymax & length(ng)>=3)
{
	groups[[i]]=ng
}else
{
	groups[[i]]=numeric()
}

}

}



# Keeping only the groups with at least three elements (i.e. three methods contribute) AND np OK
# AND not included in a bigger group



particip=rep(0,length(groups))
npOK=rep(1,length(groups))


for (i in 1:length(groups))

{

particip[i]=length(groups[[i]])



if (np[1]==0 & sum(sign(groups[[i]])>=0)<3)

{npOK[i]=0}



if (np[2]==0 & sum(sign(groups[[i]])<=0)<3)

{npOK[i]=0}


}

incl=rep(FALSE,length(groups))
iin=which(particip==3)
jin=which(particip==4)

for (i in iin)
{

	for (j in jin)
	{
		if (length(na.omit(match(groups[[i]],groups[[j]])))==3)
		{incl[i]=TRUE; break}
	}
}



output=groups[npOK*particip>=3&!incl]



return(output)

}
