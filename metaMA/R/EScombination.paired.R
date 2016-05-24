`EScombination.paired` <-
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
ES=array(dim=c(dim(logratios[[1]])[1],4,nbstudies))
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
ES[,,i]=effectsize(fit2i$t,length(logratios[[i]]),(fit2i$df.prior+fit2i$df.residual))
		}
	}
if (moderated=="SMVar")
	{
for (i in 1:nbstudies)
		{
tempLR=logratios[[i]]
stati=as.data.frame(SMVar.paired(paste("gene",rep(1:dim(tempLR)[1],1),sep=""),tempLR,threshold=BHth)) 
listgd[[i]]=stati$GeneId[which(stati$AdjPValue<=BHth)]
ES[,,i]=effectsize(stati$TestStat[order(stati$GeneId)],length(logratios[[i]]),stati$DegOfFreedom[order(stati$GeneId)])
		}
	}
if (moderated=="t")
	{
for (i in 1:nbstudies)
		{
sti=row.ttest.statp(logratios[[i]])
rpvalsti=2*(1-pt(abs(sti),df=(dim(logratios[[i]])[2]-1)))
listgd[[i]]=which(p.adjust(rpvalsti,method="BH")<=BHth)
ES[,,i]=effectsize(sti,length(logratios[[i]]),(length(logratios[[i]])-1))
		}
	}
listgd[[(nbstudies+1)]]=unique(unlist(listgd[1:nbstudies]))
#only unbiased; for biased effect sizes, add ES[,1,],ES[,2,]
restempdirect=directEScombi(ES[,3,],ES[,4,],BHth)
listgd[[(nbstudies+2)]]=restempdirect$DEindices
listgd[[(nbstudies+3)]]=restempdirect$TestStatistic
names(listgd)=c(paste("study",1:nbstudies,sep=""),"AllIndStudies","Meta","TestStatistic")
restemp=IDDIRR(listgd$Meta,listgd$AllIndStudies)
print(restemp)
invisible(listgd)
}

