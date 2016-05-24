eqsDesc <-
function(eqs)
{
	options(warn=-1)
	if(length(eqs)<=1)
	{
		eqs		<- readLines(eqs, n=-1L)
	}
	loc		<- agrep("SAMPLE STATISTICS BASED ON COMPLETE CASES",eqs)
	if(length(loc)==0)
	{
		warning("Not an appropriate .out file.  No covariance calculated.")
	}
	loc2	<- agrep("MARDIA'S COEFFICIENT",eqs)[1]
	eqsM	<- eqs[loc:(loc2-2)]
	loc		<- grep("VARIABLE",eqsM)
	v		<- NULL
	for(i in 1:length(loc))
	{
		x	<- str_trim(str_split(eqsM[loc[i]]," ")[[1]])
		x	<- x[which(nchar(x)>0 & x!="VARIABLE")]
		v	<- c(v,x)
	}
	loc		<- grep("MEAN",eqsM)
	m		<- NULL
	for(i in 1:length(loc))
	{
		x	<- str_trim(str_split(eqsM[loc[i]]," ")[[1]])
		x	<- x[which(nchar(x)>0 & x!="MEAN")]
		m	<- c(as.numeric(m),as.numeric(x))
	}
	loc		<- grep("SKEWNESS",eqsM)
	s		<- NULL
	for(i in 1:length(loc))
	{
		x	<- str_trim(str_split(eqsM[loc[i]]," ")[[1]])
		x	<- x[which(nchar(x)>0 & x!="SKEWNESS")]
		s	<- c(as.numeric(s),as.numeric(x[-1]))
	}
	loc		<- grep("KURTOSIS",eqsM)
	k		<- NULL
	for(i in 1:length(loc))
	{
		x	<- str_trim(str_split(eqsM[loc[i]]," ")[[1]])
		x	<- x[which(nchar(x)>0 & x!="KURTOSIS")]
		k	<- c(as.numeric(k),as.numeric(x[-1]))
	}
	loc		<- grep("STANDARD DEV.",eqsM)
	d		<- NULL
	for(i in 1:length(loc))
	{
		x	<- str_trim(str_split(eqsM[loc[i]]," ")[[1]])
		x	<- x[which(nchar(x)>0 & x!="DEV.")]
		d	<- c(as.numeric(d),as.numeric(x[-1]))
	}
	mat		<- data.matrix(cbind(m,d,s,k))
	rownames(mat)	<- v
	colnames(mat)	<- list("mean","sd","skew","kurt")
	return(mat)
}
