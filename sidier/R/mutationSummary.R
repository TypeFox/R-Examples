mutationSummary <-
function(align,addExtremes=FALSE,output="brief"){
ALIGN<-align
OutSeqs<-as.character(align)
align<-OutSeqs
colnames(align)<-paste("pos",c(1:ncol(align)),sep="")
if(addExtremes==TRUE)
{
align<-cbind(rep("X",nrow(align)),align,rep("X",nrow(align)))
ALIGN<-as.DNAbin(align)
}

if(output!="detailed"&output!="brief") stop("Define ouptut format!")

constants<-c()
for(i in 1:ncol(align))
if(length(unique(align[,i]))==1)
constants<-c(constants,i)
if(is.null(constants)==F) AlignConstants<-align[,c(constants)]
if(is.null(constants)) AlignConstants<-NULL

variable<-which(is.na(match(c(1:ncol(OutSeqs)),constants)))
if(length(variable)!=0) AlignVariable<-align[,variable]
if(length(variable)==0) AlignVariable<-NULL

singletons<-c()
if(is.null(AlignVariable)==FALSE)
{
AlignVariable<-as.matrix(AlignVariable)
for(i in 1:ncol(AlignVariable))
if(length(unique(AlignVariable[,i]))==2)
if(
length(which(AlignVariable[,i]==unique(AlignVariable[,i])[1]))==1|
length(which(AlignVariable[,i]==unique(AlignVariable[,i])[2]))==1)
singletons<-c(singletons,i)
}
if(is.null(singletons)==F)
{
AlignSingle<-AlignVariable[,c(singletons)]
AlignNoSingle<-AlignVariable[,-c(singletons)] #PI
}
if(is.null(singletons))
{
AlignSingle<-NULL
AlignNoSingle<-AlignVariable #PI
}
if(is.null(AlignNoSingle)==FALSE)
AlignNoSingle<-as.matrix(AlignNoSingle)
###

sinGap2<-c()
if(is.null(AlignVariable)==FALSE)
{
	for(i in 1:ncol(AlignVariable))
	{
	if(length(unique(AlignVariable[,i]))==2)
	if(
	unique(AlignVariable[,i])[1]!="-"&
	unique(AlignVariable[,i])[2]!="-")
	sinGap2<-c(sinGap2,i)
	if(is.null(sinGap2)==FALSE) AlignSinGap2<-AlignVariable[,sinGap2]
	if(is.null(sinGap2)) AlignSinGap2<-NULL
	}

if(is.null(AlignVariable))
AlignSinGap2<-NULL
}

more2<-c()
if(is.null(AlignVariable)==FALSE)
{
for(i in 1:ncol(AlignVariable))
if(length(unique(AlignVariable[,i]))>2
# & length(which(unique(AlignVariable[,i])=="-"))==0
)
more2<-c(more2,i)
if(is.null(more2)==F) AlignMore2<-AlignVariable[,more2]
if(is.null(more2)) AlignMore2<-NULL

if(is.null(more2)==FALSE | is.null(sinGap2)==FALSE) substitutions<-AlignVariable[,sort(c(more2,sinGap2))]
if(is.null(more2)==TRUE & is.null(sinGap2)==TRUE) substitutions<-NULL
}
if(is.null(AlignVariable))
substitutions<-NULL

	SubSingletons<-c()
	more2Single<-c() #singleton if A A - T A (excluding gaps)
if(is.null(AlignVariable)==FALSE & is.null(substitutions)==FALSE)
{
	substitutions<-as.matrix(substitutions)
	for(i in 1:ncol(substitutions))
	if(length(unique(substitutions[,i]))==2)
	if(
	length(which(substitutions[,i]==unique(substitutions[,i])[1]))==1|
	length(which(substitutions[,i]==unique(substitutions[,i])[2]))==1)
	SubSingletons<-c(SubSingletons,i)

	for(i in 1:ncol(substitutions))
	if(length(unique(substitutions[,i]))>2 & length(which(unique(substitutions[,i])=="-"))!=0)
		{
		more2SinGaps<-unique(substitutions[,i])[-which(unique(substitutions[,i])=="-")]		
			if(length(more2SinGaps)==2)
			if(
			length(which(substitutions[,i]==more2SinGaps[1]))==1|
			length(which(substitutions[,i]==more2SinGaps[2]))==1
			)
			more2Single<-c(more2Single,i)
		}
		more2Single<-unique(more2Single)

	if(is.null(SubSingletons)==FALSE|is.null(more2Single)==FALSE)
	{
	AlignSubSingle<-as.matrix(substitutions[,c(SubSingletons,more2Single)])
	AlignSubNoSingle<-as.matrix(substitutions[,-c(SubSingletons,more2Single)]) #PI
	}
	if(is.null(SubSingletons) & is.null(more2Single))
	{
	AlignSubSingle<-NULL
	AlignSubNoSingle<-substitutions #PI
	}
	AlignSubSingle<-as.matrix(AlignSubSingle)
	AlignSubNoSingle<-as.matrix(AlignSubNoSingle) #PI
}
if(is.null(AlignVariable) | is.null(substitutions))
{
AlignSubSingle<-NULL
AlignSubNoSingle<-NULL
}
#####

gapMAL<-c()
for (i in 1:ncol(align))
if(paste(align[,i],collapse="")==paste(rep("-",nrow(align)),collapse=""))
gapMAL<-c(gapMAL,i)
if (is.null(gapMAL)==F) 
align<-align[,-c(gapMAL)] #removes positions full of gaps
colnames(align)<-c(1:ncol(align))

seqNames<-row.names(align)

	recod<-align
	for(l2 in 1:nrow(align))
	for(l in 1:ncol(align))
		{
		if(align[l2,l]=="-")
		recod[l2,l]<-0
		if(align[l2,l]!="-")
		recod[l2,l]<-1
		}
	gapIni<-matrix(nrow=nrow(align),ncol=ncol(align))
	gapEnd<-matrix(nrow=nrow(align),ncol=ncol(align))

	if(
	length(which(recod[,ncol(recod)]==0))!=0|
	length(which(recod[,1]==0))!=0
	) stop(cat("\nSome sequences have gaps in extremes. Set 'addExtremes=TRUE' to analyse this dataset."))

	for (i2 in 1:(ncol(recod)-1))
	for (i3 in 1:nrow(recod))
	{
	if(paste(recod[i3,i2],recod[i3,(i2+1)])=="1 0")
	gapIni[i3,i2]<-i2
	if(paste(recod[i3,i2],recod[i3,(i2+1)])=="0 1")
	gapEnd[i3,i2]<-i2
	}

	GAPS<-list()
	for (i4 in 1:nrow(gapIni))
	GAPS[[i4]]<-rbind(which(is.na(gapIni[i4,])==F),which(is.na(gapEnd[i4,])==F))

GAPSjoin<-cbind(GAPS[[1]],GAPS[[2]])
for(i5 in 3:length(GAPS))
GAPSjoin<-cbind(GAPSjoin,GAPS[[i5]])
GAPSunique<-unique(GAPSjoin,MARGIN=2)
GAPSunique<-GAPSunique[,order(GAPSunique[1,],GAPSunique[2,])] #sort by size

GAPSunique<-as.matrix(GAPSunique)
if(length(GAPSunique)!=0)
{
AlGap<-as.matrix(align[,(GAPSunique[1,]+1)])
#GapSingle<-c()
#for (i in 1:ncol(AlGap))
#if(length(which(AlGap[,i]==unique(AlGap[,i])[1]))==1|
#length(which(AlGap[,i]==unique(AlGap[,i])[2]))==1)
#GapSingle<-c(GapSingle,i)
GapSingle<-c()
GAPS_PI<-c()
for(i in 1:ncol(GAPSunique))
 {
 if(length(unique(which(GAPSunique[,i]==GAPSjoin,arr.ind=T)[,2]))==1) GapSingle<-c(GapSingle,i)
 if(length(unique(which(GAPSunique[,i]==GAPSjoin,arr.ind=T)[,2]))>1) GAPS_PI<-c(GAPS_PI,i)
 }

}
if(length(GAPSunique)==0)
GapSingle<-NULL


if(is.null(GapSingle)==F) AlignGapSingle<-as.matrix(AlGap[,GapSingle])
if(is.null(GapSingle)) AlignGapSingle<-NULL
if(exists("AlGap"))
ifelse(is.null(GapSingle)==F,AlignGapNoSingle<-AlGap[,GAPS_PI],AlignGapNoSingle<-AlGap)
if(exists("AlGap")==FALSE)
AlignGapNoSingle<-NULL



####
if(addExtremes==TRUE) align<-align[,-c(1,ncol(align))]
ifelse(is.null(ncol(align)),o1<-0,o1<-ncol(align))
ifelse(is.null(length(constants)),o2<-0,o2<-length(constants))
if(is.null(length(constants))==FALSE & addExtremes==TRUE) o2<-(o2-2)
ifelse(is.null(length(variable)),o3<-0,o3<-length(variable))
ifelse(is.null(length(singletons)),o4<-0,o4<-length(singletons))
ifelse(is.null(ncol(AlignNoSingle)),o5<-0,o5<-ncol(AlignNoSingle))
OUT<-t(matrix(c(o1,o2,o3,o4,o5)))

ifelse(is.null(ncol(substitutions)),o6<-0,o6<-ncol(substitutions))
ifelse(is.null(ncol(AlignSubSingle)),o7<-0,o7<-ncol(AlignSubSingle))
ifelse(is.null(ncol(AlignSubNoSingle)),o8<-0,o8<-ncol(AlignSubNoSingle))
ifelse(is.null(ncol(GAPSunique)),o9<-0,o9<-ncol(GAPSunique))
ifelse(is.null(ncol(AlignGapSingle)),o10<-0,o10<-ncol(AlignGapSingle))
ifelse(is.null(ncol(AlignGapNoSingle)),o11<-0,o11<-ncol(AlignGapNoSingle))
OUT2<-t(matrix(c(o6,o7,o8,o9,o10,o11)))

colnames(OUT)<-c("Length","Constants","Variable","Singletons","Informative")
row.names(OUT)<-"Sites"
colnames(OUT2)<-c("Substitutions","Subst.Single", "Subst.Info","Gaps","Gaps.Single","Gaps.Info")
row.names(OUT2)<-"Events"
OUT3<-list(OUT,OUT2)
names(OUT3)<-c("Sites","Events")


if(output=="detailed")
{

ifelse(is.null(AlignConstants)==FALSE,
{
ACON<-AlignConstants
colnames(ACON)<-colnames(AlignConstants)
row.names(ACON)<-row.names(AlignConstants)
},
{
ACON<-"No constant sites"
}
)
OUT3[[3]]<-ACON
names(OUT3)[3]<-"Constants.Alignment"


ifelse(is.null(AlignVariable)==FALSE,
{
VAAR<-as.data.frame(AlignVariable)
colnames(VAAR)<-colnames(AlignVariable)
row.names(VAAR)<-row.names(AlignVariable)
},
{
VAR<-"No variables sites"
}
)
OUT3[[4]]<-VAAR
names(OUT3)[4]<-"Variables.Alignment"

if(is.null(singletons)==FALSE)
{
SING<-as.data.frame(AlignVariable[,singletons])
colnames(SING)<-colnames(AlignVariable[,singletons])
row.names(SING)<-row.names(AlignVariable[,singletons])
}
if(is.null(singletons)==TRUE)
SING<-"No singleton sites"
OUT3[[5]]<-SING
names(OUT3)[5]<-"Singletons.Alignment"


ifelse(is.null(singletons)==FALSE,
{
INF<-as.data.frame(AlignVariable[,-singletons])
colnames(INF)<-colnames(AlignVariable[,-singletons])
row.names(INF)<-row.names(AlignVariable[,-singletons])
}
,
{
	ifelse(is.null(AlignVariable)==FALSE,
		{
		INF<-as.data.frame(AlignVariable)
		colnames(INF)<-colnames(AlignVariable)
		row.names(INF)<-row.names(AlignVariable)
		}
		,
		{
		INF<-"No informative sites"
		})
}
)
OUT3[[6]]<-INF
names(OUT3)[6]<-"Informatives.Alignment"

if(is.null(c(more2,sinGap2))==FALSE)
{
SUS<-AlignVariable[,sort(c(more2,sinGap2))]
SUS<-as.data.frame(SUS)
colnames(SUS)<-colnames(AlignVariable[,sort(c(more2,sinGap2))])
row.names(SUS)<-row.names(AlignVariable[,sort(c(more2,sinGap2))])
}
if(is.null(c(more2,sinGap2))==TRUE)
SUS<-"No substitutions"
OUT3[[7]]<-SUS
names(OUT3)[7]<-"Substitutions"

if(is.null(c(more2,sinGap2))==FALSE)
{
SS<-SUS[,c(SubSingletons,more2Single)]
SS<-as.data.frame(SS)
colnames(SS)<-colnames(SUS[,c(SubSingletons,more2Single)])
row.names(SS)<-row.names(SUS[,c(SubSingletons,more2Single)])
}
if(is.null(c(more2,sinGap2))==TRUE)
SS<-"No singleton substitutions"

OUT3[[8]]<-SS
names(OUT3)[8]<-"Subst.Single"


ifelse(is.null(c(SubSingletons,more2Single))==FALSE,
{
SUI<-SUS[,-c(SubSingletons,more2Single)]
SUI<-as.data.frame(SUI)
colnames(SUI)<-colnames(SUS[,-c(SubSingletons,more2Single)])
row.names(SUI)<-row.names(SUS[,-c(SubSingletons,more2Single)])
},
{
SUI<-"No informative Substitutions"
}
)
OUT3[[9]]<-SUI
names(OUT3)[9]<-"Subst.Info"


GU<-SIC(align=ALIGN,saveFile=F,addExtremes=addExtremes)[[1]]
GU[1,]<-as.numeric(GU[1,])+1
GU<-as.data.frame(GU)
colnames(GU)<-colnames(GU)
row.names(GU)<-row.names(GU)
OUT3[[10]]<-GU
names(OUT3)[10]<-"Gaps"

ifelse(is.null(GapSingle)==FALSE,
{
GS<-as.data.frame(GU[,GapSingle])
colnames(GS)<-colnames(GU)[GapSingle]
row.names(GS)<-row.names(GU)
},
{
GS<-"No singleton gaps"
}
)
OUT3[[11]]<-GS
names(OUT3)[11]<-"Gaps.Single"

ifelse(is.null(GapSingle)==FALSE,
{
GNS<-as.data.frame(GU[,-GapSingle])
colnames(GNS)<-colnames(GU[,-GapSingle])
row.names(GNS)<-row.names(GU)
},
{
GNS<-as.data.frame(GU[,GapSingle])
colnames(GNS)<-colnames(GU[,GapSingle])
row.names(GNS)<-row.names(GU)
})

OUT3[[12]]<-GNS
names(OUT3)[12]<-"Gaps.Info"

}
OUT3
}
