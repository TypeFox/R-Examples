cq2counts <-
function(data,genecols,condcols,effic,Cq1=NULL) {
ngenes=length(genecols)
nfix=length(condcols)
if (is.null(Cq1)) { Cq1=37 } 
for (g in genecols) {
	if ((names(data)[g] %in% effic[,1])==FALSE) {
		print(paste("cq2counts error: missing efficiency for ",names(data)[g],sep=""))
		stop
	}
	e=effic[effic[,1]==names(data)[g],2]
		if (Cq1=="table"){
			cq1=effic[effic$gene==names(data)[g],3]
			if(is.na(cq1)) { cq1=37 }
		} else {
			if (Cq1=="formula") { 
				cq1=79-21.5*e
			} else { cq1=Cq1 }
		}
#	print(cq1)
	ii=which(!is.na(data[,g]))
	ii=which(data[,g]>0)
	data[ii,g]=e^(cq1-data[ii,g])
}

# stacking the dataset for modeling, to make it all have one response variable: count
c=data.frame(round(data[,genecols],0))
c[c<0]=0 # since -1 meant there was zero molecules at the start
cs=stack(c)
names(cs)=c("count","gene")
cs=cbind(cs,rep(data[,condcols]))
for (c in c(3:(2+length(condcols)))){
	cs[,c]=as.factor(cs[,c])
}
nas=which(is.na(cs$count))
if (length(nas)>0) cs=cs[-which(is.na(cs$count)),]
return(cs)
}
