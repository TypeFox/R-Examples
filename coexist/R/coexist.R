#species coexistence modeling
#under asymmetric dispersal and fluctuating source-sink dynamics
#testing the proportion of coexistence scenarios
#driven by neutral process and niche process
#WRITTEN BY YOUHUA CHEN
#COPYRIGHT: GPL>=2.0
#2012/7/13

#TWO SPECIES MODELS

#initialization
spabundance<-function(island,abund=1000)
{
sabund<-vector()
length(sabund)<-island
sabund[1]=abund
sabund[2:island]=0
return(sabund)
}

#rate matrix for each species at each island
parsetting<-function(island,rate=1,scale=2,type="decrease")
{
parset<-vector()
length(parset)<-island
if(type=="decrease")
{
parset[1:as.integer(island/2)]=rate
parset[(as.integer(island/2)+1):island]=rate/scale
}
if(type=="increase")
{
parset[1:as.integer(island/2)]=rate/scale
parset[(as.integer(island/2)+1):island]=rate
}
if(type=="constant")
{
parset[1:island]=rate
}
if(type=="mosaiclow")
{
for(i in 1:island)
{
if(i%%2==0)
{
parset[i]=rate
}else
{
parset[i]=rate/scale
}
}#for i
}
if(type=="mosaichigh")
{
for(i in 1:island)
{
if(i%%2==0)
{
parset[i]=rate/scale
}else
{
parset[i]=rate
}
}#for i
}
return(parset)
}


#dispersal parameter matrix
dispvar<-function(island,rate=.5,scale=2,type="decrease")
{
dismat<-matrix(0,ncol=island,nrow=island)
if(type=="decrease")
{
for(i in 1:as.integer((island-1)/2))
{
dismat[i,i+1]=rate
}
for(i in (as.integer((island-1)/2)+1):(island-1))
{
dismat[i,i+1]=rate/scale
}
diag(dismat)<-1-rowSums(dismat)
}
if(type=="increase")#heterogeneous increasing dispersal
{
for(i in 1:as.integer((island-1)/2))
{
dismat[i,i+1]=rate/scale
}
for(i in (as.integer((island-1)/2)+1):(island-1))
{
dismat[i,i+1]=rate
}
diag(dismat)<-1-rowSums(dismat)
}
if(type=="constant")#homogenous dispersal
{
for(i in 1:(island-1))
{
dismat[i,i+1]=rate
}
diag(dismat)<-1-rowSums(dismat)
}
if(type=="mosaiclow")#homogenous dispersal
{
for(i in 1:(island-1))
{
if(i%%2==0)
{
dismat[i,i+1]=rate
}else
{
dismat[i,i+1]=rate/scale
}
}
diag(dismat)<-1-rowSums(dismat)
}
if(type=="mosaichigh")#homogenous dispersal
{
for(i in 1:(island-1))
{
if(i%%2==0)
{
dismat[i,i+1]=rate/scale
}else
{
dismat[i,i+1]=rate
}
}
diag(dismat)<-1-rowSums(dismat)
}
return(dismat)
}

#competition parameters matrix
compvar<-function(island,rate=.5,scale=2,type="constant")
{
comp<-matrix(0,ncol=island,nrow=2)
comp[1,]=parsetting(island,rate,scale,type)
comp[2,]=1-comp[1,]
return(comp)
}



##############################
#perform dispersal analysis
dispersal<-function(spvector,initp,dismat,allee=1)
{
spvector<-spvector%*%dismat
spvector[1,1]=initp #serve as a source
spvector[2,1]=initp #serve as a source
spvector[which(spvector<allee)]=0
return(spvector)
}


#perform competition analysis
competition<-function(spvector,resource,comp1,comp2,grow,allee=1)
{
islandnum<-dim(spvector)[2]
for(i in 1:islandnum)
{
if(spvector[1,i]!=0 & spvector[2,i]!=0)
{
s1<-spvector[1,i]
s2<-spvector[2,i]
spvector[1,i]=s1+(1-comp1[1,i]*s1/resource[1,i]-comp1[2,i]*s2/resource[1,i])*s1*grow[1,i]
spvector[2,i]=s2+(1-comp2[1,i]*s2/resource[2,i]-comp2[2,i]*s1/resource[2,i])*s2*grow[2,i]
}
if(spvector[1,i]!=0 & spvector[2,i]==0)
{
s1<-spvector[1,i]
spvector[1,i]=s1+(1-comp1[1,i]*s1/resource[1,i])*s1*grow[1,i]
}
if(spvector[2,i]!=0 & spvector[1,i]==0)
{
s2<-spvector[2,i]
spvector[2,i]=s2+(1-comp2[1,i]*s2/resource[2,i])*s2*grow[2,i]
}
#check impossible coexistence
if(spvector[1,i]<allee)
{
spvector[1,i]=0
}
if(spvector[2,i]<allee) #must be a full species body and non-negative
{
spvector[2,i]=0
}
}#end for
return(spvector)
}#end competition function


#perform species-specific dispersal and fluctuating source analysis
#dismat is now a list
flex.dispersal<-function(spvector,initp,dismat,allee=1,type="constant")
{
if(type=="constant")
{
spnum<-length(dismat)
for(i in 1:spnum)
{
spvector[i,]<-spvector[i,]%*%dismat[[i]]
spvector[i,1]=initp
}
}
if(type=="flexible")
{
spnum<-length(dismat)
for(i in 1:spnum)
{
spvector[i,]<-spvector[i,]%*%dismat[[i]]
spvector[i,1]=rnorm(1,mean=initp,sd=initp/10) #fluctuating source
}
}
if(type=="cochange")
{
spnum<-length(dismat)
newresource<-rnorm(1,mean=initp,sd=initp/10)
for(i in 1:spnum)
{
spvector[i,]<-spvector[i,]%*%dismat[[i]]
spvector[i,1]=newresource #fluctuating source
}
}
spvector[which(spvector<allee)]=0
spvector[which(spvector<0)]=0
return(spvector)
}


#perform flexible competition analysis allowing multiple species
flex.competition<-function(spvector,resource,grow,comp,allee=1)
{
spnum<-dim(spvector)[1]
islandnum<-dim(spvector)[2]
for(i in 1:islandnum)
{
s<-spvector[,i]
for(sp in 1:spnum)
{
spvector[sp,i]=s[sp]+(1-comp[sp,i]*s[sp]/resource[sp,i]-(1-comp[sp,i])*sum(s[-sp])/resource[sp,i])*s[sp]*grow[sp,i]
#check impossible coexistence
if(spvector[sp,i]<allee)#must be a full species body and non-negative
{
spvector[sp,i]=0
}
}
}#end for
return(spvector)
}#end competition function




#######################################
#batch analysis of the output scenarios

#batch anlaysis of coexistence summary tables
batch.coexistence<-function(out,island=10)
{
colist<-list()
scenarionum<-length(out)
length(colist)<-scenarionum
for(i in 1:scenarionum)
{
colist[[i]]<-sta.coexistence(out[[i]],island)
}
return(colist)
}

#batch anlaysis of multiple coexistence summary tables
#depend on sta.mcoexistence
batch.mcoexistence<-function(out,island=10,spnum=2)
{
colist<-list()
scenarionum<-length(out)
length(colist)<-scenarionum
for(i in 1:scenarionum)
{
colist[[i]]<-sta.mcoexistence(out[[i]],island=island,spnum=spnum)
}
return(colist)
}


#batch analysis of niche and neutral cases
#depend on function-batch.coexistence
batch.n2n<-function(colist,island)
{
resultlist<-list()
scenarionum=length(colist)
length(resultlist)<-scenarionum
for(i in 1:scenarionum)
{
resultlist[[i]]<-sta.fitness(colist[[i]],island)
}
return(resultlist)
}


#plot distribution of niche and neutral coexistence patterns
#depend on function-batch.n2n
plot_n2n<-function(resultlist,island=10,pagesetup=c(1,1),path=NULL)
{
prop<-list()
length(prop)<-length(resultlist)
scenarionum<-length(resultlist)
temp<-matrix(0,ncol=2,nrow=island-1)
for(i in 1:scenarionum)
{
temp[,1]=rev(resultlist[[i]][1:9,1]/resultlist[[i]][island,1])
temp[,2]=rev(resultlist[[i]][1:9,2]/resultlist[[i]][island,2])
prop[[i]]<-temp
}
if(length(path)!=0)
{
randnum<-runif(1)
pos<-unlist(gregexpr("/",path))
folder<-substr(path,1,pos[length(pos)]-1)
dir.create(folder,showWarnings=FALSE)
filename=paste(path,"00yh",randnum,".dat",sep="")
}else
{
randnum<-runif(1)
dir.create(folder,showWarnings=FALSE)
filename=paste(folder,"n2nbarplot",randnum,".dat",sep="")
}
pdf(filename)
par(mfrow=pagesetup)
for(i in 1:scenarionum)
{
nnplot<-barplot(prop[[i]],beside=T,axisnames=FALSE,col=1) #don't have x-axis
y<-max(prop[[i]])
text(nnplot[5,1],y-.05,"Neutrality",font=2,cex=2)
text(nnplot[5,2],y-.05,"Niche",font=2,cex=2)
axis(1,at=nnplot[,1],labels=c(1,2,3,4,5,6,7,8,9))
axis(1,at=nnplot[,2],labels=c(1,2,3,4,5,6,7,8,9))
#axis(1,at=(nnplot[dim(temp)[1],1]+nnplot[1,2])/2,labels=paste("Scenario",i,sep="-"),lty=0)
mtext(paste("Model",i,sep="-"))
}
dev.off()
}#end function


#batch analysis to explore coexistence density for concerned two parameters
#depend on function-batch.coexistence
batch.paircomp<-function(coexistlist,parnum,parameters)
{
scenarionum<-length(coexistlist)
pairlist<-list()
length(pairlist)<-scenarionum
for(i in 1:scenarionum)
{
pairlist[[i]]<-sta.paircomparison(coexistlist[[i]],parnum=parnum,parameters=parameters)
}
return(pairlist)
}

#multiple species case
#depend on function-the output of batch.mcoexistence
#coenum determines the coexisting number of species for each patch
batch.mpaircomp<-function(coexistlist,coenum,spnum,parameters)
{
if(is.list(coexistlist))
{
scenarionum<-length(coexistlist)
pairlist<-list()
length(pairlist)<-scenarionum
for(i in 1:scenarionum)
{
pairlist[[i]]<-sta.mpaircomparison(coexistlist[[i]],coenum,spnum=spnum,parameters=parameters)
}
return(pairlist)
}
}


#batch analysis to explore coexistence density for a varying parameter
#depend on function-batch.coexistence
batch.onepar<-function(coexistlist,parameters)
{
scenarionum<-length(coexistlist)
pairlist<-list()
length(pairlist)<-scenarionum
for(i in 1:scenarionum)
{
pairlist[[i]]<-sta.comparison(coexistlist[[i]],parameters=parameters)
}
return(pairlist)
}


#batch analysis to explore multiple species coexistence density for a varying parameter
#depend on function-batch.coexistence
batch.monepar<-function(coexistlist,coenum,island,spnum,parameters)
{
if(is.list(coexistlist))
{
scenarionum<-length(coexistlist)
pairlist<-list()
length(pairlist)<-scenarionum
for(i in 1:scenarionum)
{
pairlist[[i]]<-sta.mcomparison(coexistlist[[i]],coenum,island,spnum,parameters=parameters)
}
return(pairlist)
}
}

#batch analysis to plot matrix values for one parameter matrix
#depend on make.heatmap
batch.pdf.onepar<-function(parmatlist,pagesetup=c(2,2),path=NULL)
{
scenarionum<-length(parmatlist)
parnum<-length(parmatlist[[1]])
if(length(path)!=0)
{
randnum<-runif(1)
pos<-unlist(gregexpr("/",path))
folder<-substr(path,1,pos[length(pos)]-1)
dir.create(folder,showWarnings=FALSE)
filename=paste(path,"00yh",randnum,".pdf",sep="")
}else
{
randnum<-runif(1)
dir.create(folder,showWarnings=FALSE)
filename=paste(folder,"singleparameter",randnum,".pdf",sep="")
}
pdf(filename)
par(mfrow=pagesetup)
for(each in 1:parnum)
{
for(i in 1:scenarionum)
{
xname=paste("Model",i,sep="-")
title=names(parmatlist[[i]])[each]
t<-parmatlist[[i]][[each]]
t<-t[order(t[,1],decreasing=FALSE),]
t<-t[-1,-1]
make.heatmap(t,xname=xname,xlab=c(0.1,0.25,0.5,0.75,0.9),ylab=c(1:9),title=title)
}#i
}#each
dev.off()
}


#batch analysis to plot matrix values for pairwise parameter matrix
#depend on make.heatmap
batch.pdf.pairpar<-function(parmatlist,pagesetup=c(2,2),path=NULL)
{
scenarionum<-length(parmatlist)
parnum<-length(parmatlist[[1]][[1]])
if(length(path)!=0)
{
randnum<-runif(1)
pos<-unlist(gregexpr("/",path))
folder<-substr(path,1,pos[length(pos)]-1)
dir.create(folder,showWarnings=FALSE)
filename=paste(path,"00yh",randnum,".pdf",sep="")
}else
{
randnum<-runif(1)
dir.create(folder,showWarnings=FALSE)
filename=paste(folder,"pairwiseparameters",randnum,".pdf",sep="")
}

pdf(filename)
par(mfrow=pagesetup)
for(each in 1:parnum)
{
for(i in 1:scenarionum)
{
xname=paste("Model",i,sep="-")
title=names(parmatlist[[i]][[1]])[each]
t<-parmatlist[[i]][[1]][[each]]
t<-t(t)
make.heatmap(t,xname=xname,xlab=c(0.1,0.25,0.5,0.75,0.9),ylab=c(0.1,0.25,0.5,0.75,0.9),title=title)
}#i
}#each
dev.off()
}

##############################################
####IMPORTANT SIMULATION FUNCTIONS############
##############################################
#coarse parameter sampling based on parspace
sim.coarse<-function(island=10,scale=2,dispersalscale=51,allee=1,T=1000,prange,type,initp,path=NULL)
{
parnum=5
parlen<-length(prange)
outcome<-list()
length(outcome)<-parlen^parnum
outindex<-matrix(0,nrow=parlen^parnum,ncol=parnum)
colnames(outindex)<-c("r1","r2","dis","c11","c22")

habitat1<-parsetting(island,initp,scale,type[1])
habitat2<-parsetting(island,initp,scale,type[2])
resource<-rbind(habitat1,habitat2)

#begin simulation
count=0
outcomefile=filename.check(path)
for(i1 in 1:parlen)
{
for(i2 in 1:parlen)
{
for(i3 in 1:parlen)
{
for(i4 in 1:parlen)
{
for(i5 in 1:parlen)
{
grow1<-parsetting(island,rate=prange[i1],scale,type[3])
grow2<-parsetting(island,rate=prange[i2],scale,type[4])
grow<-rbind(grow1,grow2)
dismat<-dispvar(island,rate=prange[i3]/dispersalscale,scale,type[5])
comp1<-compvar(island,rate=prange[i4],scale,type[6])
comp2<-compvar(island,rate=prange[i5],scale,type[7])
spvector<-rbind(spabundance(island,1000),spabundance(island,1000))
count=count+1
outindex[count,]<-c(prange[i1],prange[i2],prange[i3],prange[i4],prange[i5])
for(j in 1:T)
{
spvector<-competition(spvector,resource,comp1,comp2,grow,allee)
spvector<-dispersal(spvector,dismat,allee)
}
outcome[[count]]<-spvector
write.table(outcome[[count]],file=outcomefile,sep="\t",append=TRUE)
outcome[[count]]<-list(abund=outcome[[count]],pars=outindex[count,])
}#for i1
}
}
}
}#for i5
#write.table(outindex,file=outindexfile,sep="\t")
return(outcome)
}#end function


#coarse parameter sampling based on parspace, can allow multiple species
flexsim<-function(scale=2,dispersalscale=51,allee=1,T=1000,prange,initp,spnum=2,island=10,sourcetype="constant",type,path=NULL)
{
#habitat, growth rate, dispersal rates and competition rates-four parameter set models
parnum=spnum*3
outcome<-list()
if(!is.list(prange))
{
length(outcome)<-length(prange)^parnum
}else
{
parlen<-vector()
length(parlen)<-length(prange)
for(i in 1:length(prange))
{
parlen[i]<-length(prange[[i]])
}
length(outcome)<-prod(parlen)
}
resource<-matrix(0,ncol=island,nrow=spnum)
grow<-resource
comp<-resource

for(i in 1:spnum)
{
resource[i,]<-parsetting(island,initp,scale,type[i])
}

#begin simulation
outcomefile=filename.check(path)

if(!is.list(prange))
{
parcomb<-comblist(prange,parnum)
}
if(is.list(prange))#list of parameters to save space
{
parcomb<-comblist2(prange)
}
comnum<-dim(parcomb)[1]
colnames(parcomb)<-c(paste("r",c(1:spnum),sep=""),paste("dis",c(1:spnum),sep=""),paste("com",c(1:spnum),sep=""))

#construct list for storing dispersal matrix
dismat<-list()
length(dismat)<-spnum

for(each in 1:comnum)
{
typenum<-spnum
for(i in 1:spnum)
{
grow[i,]<-parsetting(island,rate=parcomb[each,i],scale,type[typenum+i])
}
typenum<-spnum+typenum
for(i in 1:spnum)
{
dismat[[i]]<-dispvar(island,rate=parcomb[each,i+spnum]/dispersalscale,scale,type[typenum+i])
}
typenum<-spnum+typenum
for(i in 1:spnum)
{
comp[i,]<-parsetting(island,rate=parcomb[each,i+spnum*2],scale,type[typenum+i])
}#finish setup of parameters
spvector<-rbind(spabundance(island,1000),spabundance(island,1000))

for(j in 1:T)
{
spvector<-flex.competition(spvector,resource,grow,comp,allee)
spvector<-flex.dispersal(spvector,dismat,allee,type=sourcetype)
}

outcome[[each]]<-spvector
write.table(outcome[[each]],file=outcomefile,sep="\t",append=TRUE)
outcome[[each]]<-list(abund=outcome[[each]],pars=parcomb[each,])

}#end parcomb search

return(outcome)
}#end function





#faster multiple species simulation
#parcomb is matrix and store all possible combinations before simulation and apply to all simulations.
fast.flexsim<-function(scale=2,island,dispersalscale=51,allee=1,T=1000,initp,parcombination,spnum=2,sourcetype="constant",type,path=NULL)
{
#habitat, growth rate, dispersal rates and competition rates-four parameter set models
parnum=spnum*3
parcomb<-parcombination
outcome<-list()
length(outcome)<-dim(parcomb)[1]
resource<-matrix(0,ncol=island,nrow=spnum)
grow<-resource
comp<-resource
spvector<-resource

for(i in 1:spnum)
{
resource[i,]<-parsetting(island,initp,scale,type[i])
}

outcomefile<-filename.check(path)

comnum<-dim(parcomb)[1]
colnames(parcomb)<-c(paste("r",c(1:spnum),sep=""),paste("dis",c(1:spnum),sep=""),paste("com",c(1:spnum),sep=""))

#construct list for storing dispersal matrix
dismat<-list()
length(dismat)<-spnum

for(each in 1:comnum)
{
typenum<-spnum
for(i in 1:spnum)
{
grow[i,]<-parsetting(island,rate=parcomb[each,i],scale,type[typenum+i])
}
typenum<-spnum+typenum
for(i in 1:spnum)
{
dismat[[i]]<-dispvar(island,rate=parcomb[each,i+spnum]/dispersalscale,scale,type[typenum+i])
}
typenum<-spnum+typenum
for(i in 1:spnum)
{
comp[i,]<-parsetting(island,rate=parcomb[each,i+spnum*2],scale,type[typenum+i])
}#finish setup of parameters
for(i in 1:spnum)
{
spvector[i,]<-spabundance(island,1000)
}
#finish abundance setup
for(j in 1:T)
{
spvector<-flex.competition(spvector,resource,grow,comp,allee)
spvector<-flex.dispersal(spvector,dismat,allee,type=sourcetype)
}

outcome[[each]]<-spvector
write.table(outcome[[each]],file=outcomefile,sep="\t",append=TRUE)
outcome[[each]]<-list(abund=outcome[[each]],pars=parcomb[each,])

}#end parcomb search

return(outcome)
}#end function

#make the parameter space matrix
make.parcomb<-function(prange,parnum,path=NULL)
{
if(!is.list(prange))
{
parcomb<-comblist(prange,parnum)
}
if(is.list(prange))#list of parameters to save space
{
parcomb<-comblist2(prange)
}
if(length(path)!=0)
{
path=filename.check(path)
write.table(parcomb,path,sep="\t")
}
return(parcomb)
}#end function



###################################
####POST-STATISTIC SECTION#########
###################################
#coexistence analysis
sta.coexistence<-function(outcome,island=10)
{
sta<-matrix(0,ncol=3+length(outcome[[1]]$pars),nrow=length(outcome))
for(i in 1:length(outcome))
{
s12=0
s1=0
s2=0
for(j in 1:island)
{
if(outcome[[i]]$abund[1,j]!=0 & outcome[[i]]$abund[2,j]!=0)
{
s12=s12+1
}
if(outcome[[i]]$abund[1,j]!=0 & outcome[[i]]$abund[2,j]==0)
{
s1=s1+1
}
if(outcome[[i]]$abund[1,j]==0 & outcome[[i]]$abund[2,j]!=0)
{
s2=s2+1
}
}#end for j
sta[i,1]=s1
sta[i,2]=s2
sta[i,3]=s12-1 #not consider the first site
sta[i,4:dim(sta)[2]]=outcome[[i]]$pars
colnames(sta)<-c("s1win","s2win","coexist","r1","r2","disp","comp1","comp2")
}#end i
return(sta)
}#end function



# coexistence analysis for multiple species case
sta.mcoexistence<-function(outcome,island=10,spnum)
{
conum<-spnum*2 #s1,s2,s3,two,three..
sta<-matrix(0,ncol=conum+length(outcome[[1]]$pars),nrow=length(outcome))
for(i in 1:length(outcome))
{
for(j in 2:island)#not consider the first patch
{
num<-length(which(outcome[[i]]$abund[,j]!=0))
if(num==0)
{
sta[i,1]=sta[i,1]+1
}
if(num==1)
{
spindex<-which(outcome[[i]]$abund[,j]!=0)
sta[i,spindex+1]=sta[i,spindex+1]+1
}
if(num>1)
{
sta[i,spnum+num]=sta[i,spnum+num]+1
}
}
sta[i,(conum+1):dim(sta)[2]]=as.numeric(outcome[[i]]$pars) #important to make a transformation
}#end i
colnames(sta)<-c(paste("s",c(0:spnum),sep=""),paste("co",c(2:spnum),sep=""),names(outcome[[i]]$pars))
return(sta)
}#end function



#neutral versus niche cases
sta.fitness<-function(coexistence,island)
{
neutral<-coexistence[which(coexistence[,4]==coexistence[,5]),]
niche<-coexistence[which(coexistence[,4]!=coexistence[,5]),]
neutral.num<-dim(neutral)[1]
niche.num<-dim(niche)[1]
conum<-matrix(0,ncol=2,nrow=island)
colnames(conum)<-c("neutral","niche")
for(i in 1:(island-1))
{
conum[i,1]<-length(which(neutral[,3]==island-i))
conum[i,2]<-length(which(niche[,3]==island-i))
}
conum[island,1]=neutral.num
conum[island,2]=niche.num
return(conum)
}

#different parameter rate comparison
sta.comparison<-function(coexistence,parameters)
{
comparisonlist<-list()
length(comparisonlist)<-dim(coexistence)[2]-3
conum<-matrix(0,ncol=length(parameters)+1,nrow=island)

for(pars in 1:length(comparisonlist))
{
for(i in 1:(island-1))
{
conum[i,1]=island-i
for(j in 1:length(parameters))
{
conum[i,j+1]<-length(which(coexistence[,3]==island-i & coexistence[,pars+3]==parameters[j]))
}#end j
}#end i
for(j in 1:length(parameters))
{
conum[island,j+1]=length(which(coexistence[,pars+3]==parameters[j]))
}
comparisonlist[[pars]]<-conum
}
names(comparisonlist)<-colnames(coexistence[,4:dim(coexistence)[2]])
return(comparisonlist)
}


#different parameter rate comparison for multiple species cases
sta.mcomparison<-function(coexistence,coenum,island,spnum,parameters)
{
comparisonlist<-list()
length(comparisonlist)<-dim(coexistence)[2]-2*spnum
conum<-matrix(0,ncol=length(parameters)+1,nrow=island)
colnames(conum)<-c("coe.num",paste("=",parameters,sep=""))

for(pars in 1:length(comparisonlist))
{
for(i in 1:(island-1))
{
conum[i,1]=island-i #number of pacthes for which there are coenum species coexist
for(j in 1:length(parameters))
{
conum[i,j+1]<-length(which(coexistence[,spnum+coenum]==island-i & coexistence[,pars+2*spnum]==parameters[j]))
}#end j
}#end i
for(j in 1:length(parameters))
{
conum[island,j+1]=length(which(coexistence[,spnum+coenum]==parameters[j]))
}
comparisonlist[[pars]]<-conum
}
names(comparisonlist)<-colnames(coexistence[,(2*spnum+1):dim(coexistence)[2]])
return(comparisonlist)
}


#pairwise parameter comparison
sta.paircomparison<-function(coexistence,parnum,parameters)
{
comparisonlist<-list()
length(comparisonlist)<-parnum*(parnum-1)/2
varlist<-comparisonlist
namesvector<-vector()
length(namesvector)<-length(comparisonlist)
count=0
for(p1 in 1:(parnum-1))
{
for(p2 in (p1+1):parnum)
{
conum<-matrix(0,ncol=length(parameters),nrow=length(parameters))
varmat<-conum
count=count+1
for(i in 1:length(parameters))
{
for(j in 1:length(parameters))
{
temp<-coexistence[which(coexistence[,3+p1]==parameters[i] & coexistence[,3+p2]==parameters[j]),]
conum[i,j]<-mean(temp[,3])
varmat[i,j]<-var(temp[,3])
}#end j
}#end i
comparisonlist[[count]]<-conum
namesvector[count]<-paste(colnames(coexistence)[3+p1],colnames(coexistence)[3+p2],sep="-")
varlist[[count]]<-varmat
}
}
names(comparisonlist)<-namesvector
names(varlist)<-namesvector
return(list(mean=comparisonlist,var=varlist))
}#end function

#different parameter rate comparison
sta.comparison<-function(coexistence,parameters,island)
{
comparisonlist<-list()
length(comparisonlist)<-dim(coexistence)[2]-3
conum<-matrix(0,ncol=length(parameters)+1,nrow=island)

for(pars in 1:length(comparisonlist))
{
for(i in 1:(island-1))
{
conum[i,1]=island-i
for(j in 1:length(parameters))
{
conum[i,j+1]<-length(which(coexistence[,3]==island-i & coexistence[,pars+3]==parameters[j]))
}#end j
}#end i
for(j in 1:length(parameters))
{
conum[island,j+1]=length(which(coexistence[,pars+3]==parameters[j]))
}
comparisonlist[[pars]]<-conum
}
names(comparisonlist)<-colnames(coexistence[,4:dim(coexistence)[2]])
return(comparisonlist)
}


#pairwise parameter comparison for multiple species with multiple parameter space
sta.mpaircomparison<-function(coexistence,coenum,spnum,parameters)
{
if(coenum>spnum | coenum<=1)
{
cat("wrong number","\n")
coenum=spnum
}
parnum<-dim(coexistence)[2]-2*spnum
comparisonlist<-list()
length(comparisonlist)<-parnum*(parnum-1)/2
varlist<-comparisonlist
namesvector<-vector()
length(namesvector)<-length(comparisonlist)
count=0
for(p1 in 1:(parnum-1))
{
for(p2 in (p1+1):parnum)
{
if(!is.list(parameters))
{
conum<-matrix(0,ncol=length(parameters),nrow=length(parameters))
colnames(conum)=paste("=",parameters,sep="")
rownames(conum)=paste("=",parameters,sep="")
varmat<-conum
count=count+1
for(i in 1:length(parameters))
{
for(j in 1:length(parameters))
{
temp<-coexistence[which(coexistence[,2*spnum+p1]==parameters[i] & coexistence[,2*spnum+p2]==parameters[j]),]
if(length(temp)==0)
{
conum[i,j]=0
varmat[i,j]=0
}else
{
conum[i,j]<-mean(temp[,spnum+coenum])
varmat[i,j]<-var(temp[,spnum+coenum])
}
}#end j
}#end i
}else #is a list of with unequal parameter spaces
{
conum<-matrix(0,ncol=length(parameters[[p2]]),nrow=length(parameters[[p1]]))
colnames(conum)=paste("=",parameters[[p2]],sep="")
rownames(conum)=paste("=",parameters[[p1]],sep="")
varmat<-conum
count=count+1
for(i in 1:length(parameters[[p2]]))
{
for(j in 1:length(parameters[[p1]]))
{
temp<-coexistence[which(coexistence[,2*spnum+p1]==parameters[[p2]][i] & coexistence[,2*spnum+p2]==parameters[[p1]][j]),]
if(length(temp)==0)
{
conum[i,j]=0
varmat[i,j]=0
}else
{
conum[i,j]<-mean(temp[,spnum+coenum])
varmat[i,j]<-var(temp[,spnum+coenum])
}
}#end j
}#end i
}#list if
comparisonlist[[count]]<-conum
namesvector[count]<-paste(colnames(coexistence)[2*spnum+p1],colnames(coexistence)[2*spnum+p2],sep="-")
varlist[[count]]<-varmat
}
}
names(comparisonlist)<-namesvector
names(varlist)<-namesvector
return(list(mean=comparisonlist,var=varlist))
}#end function



#####################
#I/O FUNCTION SECTION
#specific for this package
#all dataset's name are coded
#as a xxxx00yh0.xxxx.dat format

read.patchdata<-function(path=NULL,spnum=2,islandnum=10)
{
if(length(path)!=0)
{
raw<-scan(path,what=character(),sep="\t")
fileline<-length(count.fields(path))
outlist<-list()
length(outlist)<-fileline/(spnum+1)

sp<-matrix(0,nrow=spnum,ncol=islandnum)

count=0
for(i in 1:length(raw))
{
if(raw[i]=="V1")
{
count=count+1
for(j in 1:spnum)
{
sp[j,]=as.numeric(raw[(i+10+j+(j-1)*islandnum):(i+10+j-1+j*islandnum)])
}
outlist[[count]]<-list(abund=sp)
}
}
}
return(outlist)
}#end function

#read data with index
read.data<-function(path=NULL,index=NULL,spnum=2,islandnum=10)
{
if(is.character(index))
{
indmat<-read.table(index,header=TRUE)
index=indmat
}

if(length(path)!=0)
{
raw<-scan(path,what=character(),sep="\t")
fileline<-length(count.fields(path))
outlist<-list()
length(outlist)<-fileline/(spnum+1)

sp<-matrix(0,nrow=spnum,ncol=islandnum)

count=0
for(i in 1:length(raw))
{
if(raw[i]=="V1")
{
count=count+1
for(j in 1:spnum)
{
sp[j,]=as.numeric(raw[(i+10+j+(j-1)*islandnum):(i+10+j-1+j*islandnum)])
}
#read pars
par<-index[count,] #save parameters
names(par)<-c(paste("r",1:spnum,sep=""),paste("dis",1:spnum,sep=""),paste("c",1:spnum,sep=""))
#save both abundance and pars data
outlist[[count]]<-list(abund=sp,pars=par)
}
}
}
return(outlist)
}#end function

#batch read different file data based on the order 1,2,3,4,5,6...
#it can directly read data from the raw file names without trucation
batch.read<-function(path,index=NULL,spnum=2,islandnum=10)
{
num<-length(path)
outlist<-list()
length(outlist)<-num
if(length(index)==0)
{
for(i in 1:num)
{
outlist[[i]]<-read.patchdata(path=path[i],spnum=spnum,islandnum=islandnum)
}
}
if(length(index)!=0)
{
for(i in 1:num)
{
outlist[[i]]<-read.data(path=path[i],index=index, spnum=spnum,islandnum=islandnum)
}
}
return(outlist)
}

#convert saved data sets' names in to a list
#based on the order of numbers
convert.filenames<-function(folder)
{
files<-list.files(path=folder,full.names=TRUE)
newf<-vector()
fnum<-length(files)
length(newf)<-fnum
dataorder<-rep(0,1,fnum)
for(i in 1:fnum)
{
pos.end<-unlist(gregexpr("00yh",files[i]))[1]-1
pos.start<-unlist(gregexpr(paste(folder,"/out",sep=""),files[i]))[1]+15
dataorder[i]<-as.numeric(substr(files[i],pos.start,pos.end))
}
for(i in 1:fnum)
{
newf[dataorder[i]]<-files[i]
}
newf<-newf[!is.na(newf)]
return(newf)
}


#open and create a new folder if not existed
filename.check<-function(path=NULL,return=TRUE)
{
randnum<-runif(1)
if(length(path)!=0)
{
pos<-unlist(gregexpr("/",path))
dot<-unlist(gregexpr(".",path,fixed=TRUE))
dd<-unlist(gregexpr(":",path,fixed=TRUE))
special<-unlist(gregexpr("00yh",path))
if(length(pos)>=2 & dot[1]!=-1 & dd[1]!=-1 & special[1]!=-1)
{
folder<-substr(path,1,pos[length(pos)]-1)
dir.create(folder,showWarnings=FALSE)
fname<-path
}
if(length(pos)>=2 & dot[1]!=-1 & dd[1]!=-1 & special[1]==-1)
{
folder<-substr(path,1,pos[length(pos)]-1)
dir.create(folder,showWarnings=FALSE)
fname<-paste(substr(path,1,dot[length(dot)]-1),"00yh",substr(path,dot[length(dot)],nchar(path)),sep="")
}
if(dot[1]==-1 & dd[1]!=-1 & special[1]==-1)
{
dir.create(path,showWarnings=FALSE)
fname<-paste(path,"/",randnum,"00yh",".dat",sep="")
}
if(dot[1]==-1 & dd[1]==-1 & special[1]==-1)
{
path<-paste("c://",path,sep="")
dir.create(path,showWarnings=FALSE)
fname<-paste(path,"/",randnum,"00yh",".dat",sep="")
}
}
if(length(path)==0)
{
if(length(folder)==0)
{folder="c://outcome"}
fname<-paste(folder,"/",randnum,"00yh",".dat",sep="")
dir.create(folder,showWarnings=FALSE)
}
if(return==TRUE)
{
return(fname)
}
}#end check

######################
#PLOT FUNCTION SECTION
######################
#make a heatmap based on matrix values
make.heatmap<-function(mat,type="gray",xname="x",yname="y",xlab=NULL,ylab=NULL,title=NULL)
{
xrange<-dim(mat)[2]+1
yrange<-dim(mat)[1]+1
if(xrange>=yrange)
{
maxrange<-xrange
xscale=1
yscale<-yrange/xrange
}
if(yrange>xrange)
{
maxrange<-yrange
xscale<-xrange/yrange
yscale=1
}
maxnum=ceiling(max(as.vector(mat)))
tmat<-t(mat)
plot((1:maxrange)*xscale,(1:maxrange)*yscale,type="n",xlab=xname,ylab=yname,axes=FALSE)
if(length(xlab)!=0)
{
axis(1,at=c(1:(xrange-1))+.5,labels=xlab)
}else
{
axis(1,at=c(1:(xrange-1))+.5,labels=c(1:(xrange-1)))
}
if(length(ylab)!=0)
{
axis(2,at=c(1:(yrange-1))+.5,labels=ylab)
}else
{
axis(2,at=c(1:(yrange-1))+.5,labels=c(1:(yrange-1)))
}
if(length(title)!=0)
{
mtext(title)
}
for(i in 1:(xrange-1))
{
for(j in 1:(yrange-1))
{
value=tmat[i,j]
x1=i
x2=i+1
y1=j
y2=j+1
if(type=="gray")
{
rect(x1,y1,x2,y2,col=gray((maxnum-value)/maxnum),border=gray((maxnum-value)/maxnum))
}
}#j
}#i
}#end function

################
#OTHER FUNCTIONS
################
#construct a full list of all combinations of parameters
#parvecotr is a parameter vector, vector
#parnum is the number of parameters, numeric
comblist<-function(parvector,parnum)
{
combnum<-length(parvector)^parnum
mat<-matrix(0,ncol=parnum,nrow=combnum)
for(i in 1:parnum)
{
leg<-length(parvector)^(parnum-i)
period=length(parvector)^i
repeated=length(parvector)^(i-1)
fullcircle=leg*length(parvector)
temp<-vector()
length(temp)<-fullcircle
for(k in 1:length(parvector))
{
temp[((k-1)*leg+1):(k*leg)]=parvector[k] 
}

for(ii in 1:repeated)
{
mat[((ii-1)*fullcircle+1):(ii*fullcircle),i]=temp
}
}
return(mat)
}#end function




#construct a full list of all combinations of parameters
#parvecotr is a parameter list, can vary on size
#parnum is the number of parameters, numeric
comblist2<-function(parlist)
{
combnum=1
parnum<-length(parlist)
lvector<-vector()
length(lvector)<-parnum
for(i in 1:parnum)
{
combnum<-length(parlist[[i]])*combnum
lvector[i]=length(parlist[[i]])
}
mat<-matrix(0,ncol=parnum,nrow=combnum)
for(i in 1:parnum)
{
leg<-prod(lvector[i:parnum])/lvector[i]
repeated=combnum/prod(lvector[i:parnum])
fullcircle=leg*lvector[i]
temp<-vector()
length(temp)<-fullcircle
for(k in 1:length(parlist[[i]]))
{
temp[((k-1)*leg+1):(k*leg)]=parlist[[i]][k] 
}

for(ii in 1:repeated)
{
mat[((ii-1)*fullcircle+1):(ii*fullcircle),i]=temp
}
}
return(mat)
}#end function