cq2log <-
function(data,genecols,condcols,effic,noamp=38) {
ngenes=length(genecols)
nfix=length(condcols)

for (g in genecols) {
	if ((names(data)[g] %in% effic[,1])==FALSE) {
		print(paste("cq2counts error: missing efficiency for ",names(data)[g],sep=""))
		stop
	}
	e=effic[effic[,1]==names(data)[g],2]
	if (length(which(data[,g]==(-1))>0)) {
		warning(paste("Warning: ", length(which(data[,g]==(-1))), " no-amplification instances for gene ",names(data)[g],"; these will be converted to Cq=",noamp,". Consider using Poisson-lognormal model.",sep=""))
		data[which(data[,g]==(-1)),g]=noamp
	}
	ii=which(!is.na(data[,g]))
	data[ii,g]=data[ii,g]*log(e)*(-1)
}
# stacking the dataset for modeling, to make it all have one response variable: count
c=data.frame(round(data[,genecols],2))
cs=data.frame(stack(c))
names(cs)=c("count","gene")
cs=cbind(cs,rep(data[,condcols]))
for (c in c(3:(2+length(condcols)))){
	cs[,c]=as.factor(cs[,c])
}
nas=which(is.na(cs$count))
if (length(nas)>0) cs=cs[-which(is.na(cs$count)),]
return(cs)
}
