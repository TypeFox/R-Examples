pop.dist <-
function(DistFile=NA,distances=NA,HaploFile=NA,Haplos=NA,outType="O",logfile=TRUE, saveFile=TRUE,NameIniPopulations=NA,NameEndPopulations=NA,NameIniHaplotypes=NA,NameEndHaplotypes=NA)
{
if(is.na(HaploFile)==TRUE&is.na(Haplos[1])==TRUE) print("Error: Please, define either Haplotypes or Haplotype file")
if(is.na(HaploFile)==FALSE&is.na(Haplos[1])==FALSE) print("Error: Please, define either Haplotypes or Haplotype file")
if(is.na(HaploFile)==FALSE&is.na(Haplos[1])==TRUE)
Haplos<-read.table(file=HaploFile)

Haplos<-data.frame(Haplos)

if(is.na(NameEndPopulations))
{
NameIniPopulations<-1
NameEndPopulations<-max(nchar(row.names(Haplos)))
}
pops<-unique(substr(row.names(Haplos),NameIniPopulations,NameEndPopulations))

if(is.na(DistFile)==TRUE&is.na(distances[1])==TRUE) print("Error: Please, define either distances or distances file")
if(is.na(DistFile)==FALSE&is.na(distances[1])==FALSE) print("Error: Please, define either distances or distances file")
if(is.na(DistFile)==FALSE&is.na(distances[1])==TRUE)
distances<-read.table(DistFile)

uniques<-c() ##get haplotypes distances using only unique haplotypes
if(is.na(NameEndHaplotypes))
{
NameIniHaplotypes<-1
NameEndHaplotypes<-max(nchar(row.names(distances)))
}
for(i in 1:length(unique(substr(row.names(distances),NameIniHaplotypes,NameEndHaplotypes))))
uniques<-c(uniques,which(unique(substr(row.names(distances),NameIniHaplotypes,NameEndHaplotypes))[i]==substr(row.names(distances),NameIniHaplotypes,NameEndHaplotypes))[1])
distances<-distances[c(uniques),c(uniques)]


POPdist<-matrix(NA,nrow=length(pops),ncol=length(pops))
row.names(POPdist)<-pops
colnames(POPdist)<-pops

for (k in 1:nrow(POPdist))
POPdist[k,k]<-0

for(k in 1:(nrow(Haplos)-1))
for(l in (k+1):nrow(Haplos))
{
pop1<-Haplos[k,]
pop2<-Haplos[l,]

Hpop1<-rep(colnames(pop1)[which(pop1!=0)],pop1[which(pop1!=0)])
Hpop2<-rep(colnames(pop2)[which(pop2!=0)],pop2[which(pop2!=0)])

DIST1_2<-c()	
	for(i in 1:length(Hpop1))
	for(j in 1:length(Hpop2))
	DIST1_2<-c(DIST1_2,distances[which(substr(row.names(distances),NameIniHaplotypes,NameEndHaplotypes)==Hpop1[i]), which(substr(colnames(distances),NameIniHaplotypes,NameEndHaplotypes)==Hpop2[j])])
	if(length(DIST1_2)==0) DIST1_2<-0
	if(outType=="7"|outType=="O")
	POPdist[k,l]<-mean(DIST1_2)
	if(outType=="L"|outType=="O")
	POPdist[l,k]<-mean(DIST1_2)

	row.names(POPdist)<-pops
	colnames(POPdist)<-pops
}
if(logfile==T)
if(is.na(DistFile)==T) DistFile<-"from memory"
if(is.na(HaploFile)==T) HaploFile<-"from memory"
write.table(c(paste("Among haplotypes distance matrix used:",DistFile),paste("Haplotypes per population matix used:",HaploFile)),file="PopopulationDistances.r.txt.log",row.names=F,col.names=F, quote=F)

if(saveFile==T)
write.table(POPdist,file="PopulationDistances.r.txt",na="",row.names=F,col.names=F, quote=F)
as.data.frame(POPdist,nrow=nrow(POPdist))
}
