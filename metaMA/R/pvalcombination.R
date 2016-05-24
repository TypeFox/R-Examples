`pvalcombination` <-
function(esets, classes, moderated=c("limma","SMVar","t")[1],BHth=0.05) 
{
nbstudies=length(esets)
if (!(moderated %in% c("limma","SMVar","t")))
	{
print("Wrong argument for \"moderated\" in pvalcombi->by default,
 limma moderated t-tests will be used")
moderated="limma"
	}
if (nbstudies != length(classes)) 
        stop("Length of classes must be equal to length of esets.")
    for (i in 1:nbstudies) {
if(length(which(apply(esets[[i]],1,FUN=function(x) sum(is.na(x)))[1:10]==dim(esets[[i]])[2]))!=0)
{stop("Please delete genes with complete missing data in at least one of the studies.
 Only missing at random values are allowed in this package")}
        if (!is.factor(classes[[i]])) {
            classes[[i]] <- factor(classes[[i]])
        }
        if (nlevels(classes[[i]]) != 2) {
            stop("Error: Each list in the argument \"classes\" must contain exactly 2 levels.")
        }
        else {
            Ref <- levels(classes[[i]])[1]
            classes[[i]] <- sapply(classes[[i]], function(x) ifelse(x == 
                Ref, 0, 1))
        	}
	}
listgd=vector("list", (nbstudies+3))
if (moderated=="limma")
	{
for (i in 1:nbstudies)
		{
group <- as.factor(classes[[i]])
design <- model.matrix(~-1 + group)
fit = lmFit(esets[[i]], design)
contrast.matrix <- makeContrasts("group0 - group1", levels = design)
fit2i <- contrasts.fit(fit, contrast.matrix)
fit2i <- eBayes(fit2i)
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
tempC1=esets[[i]][,which(classes[[i]]==1)]
tempC2=esets[[i]][,which(classes[[i]]==0)]
stati=as.data.frame(SMVar.unpaired(paste("gene",rep(1:dim(tempC1)[1],1),sep=""), 
list(tempC1,tempC2),threshold=BHth)) 
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
sti=row.ttest.stat(esets[[i]][,which(classes[[i]]==1)],
esets[[i]][,which(classes[[i]]==0)])
p1sidedsti=pt(sti,df=(length(classes[[i]])-2))
assign(paste("p1sidedst",i,sep=""),p1sidedsti)
rpvalsti=2*(1-pt(abs(sti),df=(length(classes[[i]])-2)))
listgd[[i]]=which(p.adjust(rpvalsti,method="BH")<=BHth)
		}
tempvec=paste("p1sidedst",1:nbstudies,sep="")
	}
lsinglep=lapply(tempvec,FUN=function(x) get(x,inherits=TRUE))
nrep=unlist(lapply(classes,FUN=function(x)length(x)))
listgd[[(nbstudies+1)]]=unique(unlist(listgd[1:nbstudies]))
restempdirect=directpvalcombi(lsinglep,nrep,BHth)
listgd[[(nbstudies+2)]]=restempdirect$DEindices
listgd[[(nbstudies+3)]]=restempdirect$TestStatistic
names(listgd)=c(paste("study",1:nbstudies,sep=""),"AllIndStudies","Meta","TestStatistic")
restemp=IDDIRR(listgd$Meta,listgd$AllIndStudies)
print(restemp)
invisible(listgd)
}

