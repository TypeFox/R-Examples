ridit.raw <-
function(x,g,ref=NULL){
# x must be numeric vector but g and ref could be character
# ref= NULL means total of all groups used as reference group (length(ref)==0)
# ref= a single value from g means that group used as reference group (length(ref)==1)
# ref = a vector means that vector used as reference group (length(ref)>1)
x=as.numeric(x)
x=as.vector(x)
g=as.factor(g)			
levels=levels(g)			
levels(g)=1:length(levels)	
g=as.vector(g)
g=as.character(g)
code=is.numeric(ref)
ref=as.vector(ref)
ref=as.character(ref)
# checking reference group
if (length(ref)>1){
x=c(x,ref)
g=c(g,rep(".ref",length(ref)))
levels=c(".ref",levels)
}
# create contingency table
crosstab=t(as.matrix(table(x,g)))
rownames(crosstab)=levels		
# extract reference row from contingency table
refindex=NULL
if (length(ref)==1){
if(!code) refindex=which(levels==ref) ######
if(code && ref>=1 && ref<=nrow(crosstab)) refindex=as.numeric(ref) ######
} else if (length(ref)>1) refindex=which(levels==".ref")
# determine refrow using refindex variable
if(length(refindex)!=0) refrow=crosstab[refindex,] 
else refrow=apply(crosstab,2,sum)
# write message for announce reference row
if(length(refindex)==0) msg=paste("Reference: Total of all groups",sep="")
else msg=paste("Reference: Group = ",refindex,", Label = ", levels[refindex],sep="")
#######
# start calculate Ridit
nref=sum(refrow)
ridit=0.5*refrow[1]/nref
for(i in 2:length(refrow)){
iridit=(sum(refrow[1:i-1])+0.5*refrow[i])/nref
ridit=c(ridit,iridit)
}
# end calculate Ridit
# start calculate Mean Ridit
n=apply(crosstab,1,sum)
meanRidit=c()
for(i in 1:nrow(crosstab)){
itable=crosstab[i,]
meanRidit=c(meanRidit,sum(ridit*itable)/n[i])
}
# end calculate Mean Ridit
# start calculate rbar0 using formula 3.41 
n0=sum(n)
rbar0=sum(n*meanRidit)/n0
# end calculate rbar0 using formula 3.41 
# start calculate f adjustment factor using 3.37 formula
t=apply(crosstab,2,sum)
f=1-(sum(t*(t-1)*(t+1)))/(n0*(n0-1)*(n0+1))
# end calculate f adjustment factor using 3.37 formula
# start calculate chi-sq test statistic using 3.42 formula
teststatistic=(12*n0*sum(n*(meanRidit-rbar0)^2))/((n0+1)*f)
testdf=nrow(crosstab)-1
pvalue=pchisq(q=teststatistic,df=testdf,lower.tail = FALSE)
#  calculate chi-sq test statistic using 3.42 formula
# start making output summary
if(length(ref)==0) ref=NULL
names(meanRidit)=rownames(crosstab)
output=list(MeanRidit=meanRidit,TestStatistic=teststatistic,df=testdf,Sig=pvalue,x=x,g=g,ref=ref,crosstab=crosstab,msg=msg)
class(output) <- c("ridit", class(output))
output
}
