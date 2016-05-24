`pvalcombination.paired` <-
function(logratios, moderated=c("limma","SMVar","t")[1],BHth=0.05) 
{
nbstudies=length(logratios)
if (!(moderated %in% c("limma","SMVar","t")))
	{
print("Wrong argument for moderated in pvalcombi->by default,
 limma moderated t-tests will be used")
moderated="limma"
	}
listgd=vector("list", (nbstudies+3))
if (moderated=="limma")
	{
for (i in 1:nbstudies)
		{
if(length(which(apply(logratios[[i]],1,FUN=function(x) sum(is.na(x)))[1:10]==dim(logratios[[i]])[2]))!=0)
{stop("Please delete genes with complete missing data in at least one of the studies.
 Only missing at random values are allowed in this package")}
design <- rep(1,dim(logratios[[i]])[2])
fit = lmFit(logratios[[i]], design)
fit2i <- eBayes(fit)
listgd[[i]]=which(p.adjust(fit2i$p.value,method="BH")<=BHth)
p1sidedLimma=pt(fit2i$t,df=(fit2i$df.prior+fit2i$df.residual))
assign(paste("p1sidedLimma",i,sep=""),p1sidedLimma)
		}
tempvec=paste("p1sidedLimma",1:nbstudies,sep="")
	}
if (moderated=="SMVar")
	{
for (i in 1:nbstudies)
		{
tempLR=logratios[[i]]
stati=as.data.frame(SMVar.paired(paste("gene",rep(1:dim(tempLR)[1],1),sep=""),tempLR,threshold=BHth)) 
listgd[[i]]=stati$GeneId[which(stati$AdjPValue<=BHth)]
p1sidedSMVartemp=as.vector(pt(stati$TestStat,stati$DegOfFreedom))
p1sidedSMVar=p1sidedSMVartemp[order(stati$GeneId)]
assign(paste("p1sidedSMVar",i,sep=""),p1sidedSMVar)
		}
tempvec=paste("p1sidedSMVar",1:nbstudies,sep="")
	}
if (moderated=="t")
	{
for (i in 1:nbstudies)
		{
sti=row.ttest.statp(logratios[[i]])
p1sidedsti=pt(sti,df=(dim(logratios[[i]])[2]-1))
assign(paste("p1sidedst",i,sep=""),p1sidedsti)
rpvalsti=2*(1-pt(abs(sti),df=(dim(logratios[[i]])[2]-1)))
listgd[[i]]=which(p.adjust(rpvalsti,method="BH")<=BHth)
		}
tempvec=paste("p1sidedst",1:nbstudies,sep="")
	}
lsinglep=lapply(tempvec,FUN=function(x) get(x,inherits=TRUE))
nrep=unlist(lapply(logratios,FUN=function(x) dim(x)[2]))
listgd[[(nbstudies+1)]]=unique(unlist(listgd[1:nbstudies]))
restempdirect=directpvalcombi(lsinglep,nrep,BHth)
listgd[[(nbstudies+2)]]=restempdirect$DEindices
listgd[[(nbstudies+3)]]=restempdirect$TestStatistic
names(listgd)=c(paste("study",1:nbstudies,sep=""),"AllIndStudies","Meta","TestStatistic")
restemp=IDDIRR(listgd$Meta,listgd$AllIndStudies)
print(restemp)
invisible(listgd)
}

