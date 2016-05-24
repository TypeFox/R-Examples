multinom.stat <-
function(obs,pop=1,zloc)# obs are the categorical data
{
	fobs=factor(obs,levels=unique(obs))
	Ck=table(fobs);C=sum(Ck)
	Ckz=table(fobs[zloc]);CZ=sum(Ckz)
	xz=Ckz*log(Ckz/CZ)+(Ck-Ckz)*log(((Ck-Ckz)/(C-CZ)))-Ck*log(Ck/C)
	xz[is.nan(xz)]=0
	t.stat=sum(xz)
	p1=(Ckz/CZ);p10=(Ck-Ckz)/(C-CZ)	
	#equivalent to minimize vz
	return(c(t.stat,p10,p1))
}
