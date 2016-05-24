eqsCov <-
function(eqs)
{
	if(length(eqs)<=1)
	{
		eqs		<- readLines(eqs, n=-1L)
	}
	options(warn=-1)
	loc		<- agrep("COVARIANCE MATRIX TO BE ANALYZED",eqs)
	loc2	<- grep("BENTLER-WEEKS STRUCTURAL REPRESENTATION",eqs)
	eqsCov	<- eqs[loc:(loc2-1)]
	loc		<- which(is.na(str_locate("",eqsCov)[,1])==FALSE)
	drop	<- NULL
	for(l in 1:(length(loc)-1))
	{
		if((loc[l]+1)==loc[l+1])
		{
			drop	<- c(drop,l)
		}
	}
	loc		<- loc[-drop]
	v		<- NULL
	for(l in 1:length(loc))
	{
		tmp		<- str_trim(str_split(eqsCov[loc[l]+1]," ")[[1]])
		tmp		<- tmp[which(nchar(tmp)>0)]
		if(length(which(is.na(tmp)))>0)
		{
			tmp		<- tmp[-which(is.na(tmp))]
		}
		v		<- c(v,tmp)
	}
	mat		<- matrix(0,length(v),length(v))
	colnames(mat)	<- v
	rownames(mat)	<- v
	x		<- NULL
	cnt		<- 1
	for(i in 1:(length(loc)-1))
	{
		tmp		<- eqsCov[loc[i]:loc[i+1]]
		lt		<- str_trim(str_split(tmp[grep(v[cnt],tmp)[1]]," ")[[1]])
		lt		<- lt[-which(lt=="")]
		cnt		<- cnt+length(lt)
		poi		<- grep("[.]",tmp)
		tm		<- matrix(0,length(poi),length(lt))
		colnames(tm)	<- lt
		for(j in 1:length(poi))
		{
			p		<- str_split(str_trim(tmp[poi[j]])," ")[[1]]
			p		<- p[-which(p=="")]
			p		<- as.vector(p[-which(is.na(as.numeric(p)))],"numeric")
			tm[j,]	<- c(p,rep(0,length(lt)-length(p)))
		}
		for(j in 1:length(lt))
		{
			x		<- which(lt[j]==colnames(mat))
			mat[,x]	<- c(rep(0,length(v)-length(tm[,j])),tm[,j])
		}
	}
	x		<- diag(mat)
	mat		<- mat+t(mat)
	diag(mat)	<- x
	return(mat)
}
