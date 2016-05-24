rmulti.one <-
function(size,p)
{
	re=as.numeric(rmultinom(1, size = size, prob = p))
	if (length(which(re>1))==0) return(re)
	else
	{
		floc=NULL #fix location
		repeat{
		w2=which(re>=1) #select the locations should be fixed
		floc=unique(c(floc,w2))
		s2=sum(re[floc])-length(floc)
		r2=c(1:length(p))[-floc]
		np=p[-floc]/sum(p[-floc])
		re2=rmultinom(1, size = s2, prob = np)
		re[floc]=1
		re[r2]=re2
		if(length(which(re>1))==0)break 
		}
		return(re)
	}
}
