cq2genorm <-
function(data,genes,effic,noamp=38) {
ngenes=length(genes)
for (g in genes) {
	if ((g %in% effic[,1])==FALSE) {
		stop(paste("missing efficiency for ",g,sep=""))
	}
	e=effic[effic[,1]==g,2]
	if (length(which(data[,g]==(-1))>0)) {
		warning(paste("Warning: ", length(which(data[,g]==(-1))), " no-amplification instances for gene ",g,"; these will be converted to Cq=",noamp,sep=""))
		data[which(data[,g]==(-1)),g]=noamp
	}
	ii=which(!is.na(data[,g]))
	data[ii,g]=data[ii,g]*log(e)*(-1)
}

c=data[,genes]
for (g in names(c)){
	c[,g]=c[,g]-mean(c[,g],na.rm=T)
}

# removing rows with NAs
c.nona=data.frame()
for (r in c(1:length(c[,1]))){
	if (sum(is.na(c[r,]))) next
	c.nona=data.frame(rbind(c.nona,c[r,]))
}
return(exp(c.nona))
}
