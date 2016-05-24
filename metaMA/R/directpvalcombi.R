`directpvalcombi` <-
function(pvalonesided,nrep,BHth=0.05) 
{
listres=vector("list",2)
nbstudies=length(pvalonesided)
nbreptot=sum(nrep)
if (nbreptot <2) {stop("Error: the argument \"nrep\" must be a vector with at least two values higher than 1")} 
weight=sqrt(nrep/nbreptot)
fstatmeta=function(g){
vecptime=unlist(lapply(pvalonesided, FUN = function(x) x[g]))
vec = qnorm(1 - vecptime)
stattestg = sum(weight[1:length(pvalonesided)] * vec[1:length(pvalonesided)], na.rm = TRUE)
stattestg}
statpvalc=unlist(lapply(rep(1:length(as.vector(pvalonesided[[1]])), 1), function(x) fstatmeta(x)))
rpvalpvalc=2*(1-pnorm(abs(statpvalc)))
res=which(p.adjust(rpvalpvalc,method="BH")<=BHth)
listres[[1]]=res
listres[[2]]=statpvalc
names(listres)=c("DEindices","TestStatistic")
listres
}

