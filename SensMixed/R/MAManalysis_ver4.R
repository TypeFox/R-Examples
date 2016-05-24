# This is the MAManalysis function related to the contents of the two papers: 
#
# Brockhoff P.B., Schlich, P. and Skovgaard, I.M. (2013) Taking individual scaling 
# differences into account by analyzing prfile data with the Mixed Assessor Model. 
# Submitted to Food Quality and Preference. 
#
# C. Peltier, M. Visali, P. B. Brockhoff & Schlich, P. (2013). The mam-cap
# table: a new tool for monitoring panel performances. Food Quality and
# Preference. 
library(doBy)

MAManalysis=function(data,adjustedMAM=TRUE,alpha_conditionalMAM=0.2){
  # data:
  # data should be a data.frame of the following very specific form:  
  # The data MUST have the structure: ass prod rep att1 att2 ....attP
  # There MUST be at least TWO attribute in the data frame "data"
  # There MUST be at least 3 products
  # There MUST be at least 2 replicates
  # The data must be 100% complete and balanced (NO missing values)
  
  # adjustedMAM:  
  # The option adjustedMAM controls whether or not the negativity-adjustment 
  # described in Brockhoff et al (2013) is used. Default choice is TRUE meaning that 
  # the adjusted MAM is used by default. 
  # The option does not affect the result arrays number 1, 2 and 6.  
  
  # alpha_conditionalMAM:  
  # Controls the conditional MAM by only using the MAM if the P-value for the test
  # of scaling differences is below this value. If the adjustedMAM=TRUE the conditional
  # check of scaling differences is using the adjusted mean squares. The default level
  # is 0.20 as recommended in Brockhoff et al (2013). If this is set to 1 it corresponds
  # to the unconditional MAM, that is, applying the MAM always, no matter the degree of 
  # observed scaling differences
  # The option does not affect the result arrays number 1, 2.  
  #
  # Use e.g. the function by the R call: > result=MAManalysis(data) 
  #
  # And the function produces a list of 6 arrays    
  # And you have the following result-arrays: 
  # (refered to in R with double squared brackets - as shown)
  # For the explanations below: 
  # I= number of assessors, J= number of products, K=number of replicates
  # 
  #
  # result[[1]]:
  # A 3-way array with the individually decomposed ANOVA tables for each attribute, 
  # cf. Table 2 in Brockhoff et al (2013).
  # So result[[1]][,,1] is a 5xI matrix containing the individualized ANOVA table for attribute 1 
  # This result is not affected by the two options - giving always pure MAM decomposition 
  #
  # result[[2]]:
  # A 3-way array with individual performance tests for each attribute corresponding  
  # to Table 2 of Peltier et al (2013). (result[[2]][,,1] contains the results for attribute 1)
  # P-values are categorized using the usual R-symbols:
  # " "="p-value>=0.1", ."="p-value <0.1", "*"="p-value <0.05", "**"="p-value <0.01", "***"="p-value <0.001" 
  # The first four (double) row of these result matrices show the results 
  # corresponding exactly to the four rows of Table 2 in Peltier et al (2013)
  # The MAM-CAP table as such is NOT produced. 
  # Instead a descriptive statistic is given for each performance measure:
  #
  # PRODUCT DISCRIMINABILITY: The square root of the individual assessor product F
  # SCALING: The individual beta values (averaging to 1)
  # DISAGREEMENT: Individual disagreement statistic: sqrt(SS_dis(i)/(J-2))
  #             Will average to approximately sqrt(K*MS_DIS) 
  # REPEATABILITY: The individual error (within product) standard deviation
  #
  # In addition we provide two more statistics and hypothesis tests:
  #
  # LEVEL: The main effect of assessor (summing to zero) 
  #       We test whether the individual is different from the average
  # CORRELATION: The correlation between the individual product averages
  #       and the overall(consensus) product average. These will average to 
  #       something close to the so-called Cronbach's Alpha.  
  #       We test whether the correlation is different from zero, i.e. it 
  #       can also be seen as the significance test for negativity (if the 
  #       correlation is in fact negative)
  #
  # This result is not affected by the two options - giving always pure MAM 
  # decomposition based perfromance results 
  #
  #
  # result[[3]]:
  # A 3-way array with the MAM ANOVA table for each attribute, corresponding  
  # to Table 5 in Brockhoff et al (2013) and Table 1 in Peltier et al (2013).
  # The assessor effect is tested versus the (standard) assessor*product interaction MS
  # The product effect is tested versus the MS(disagreement) OR if the 
  # conditional MAM is used versus the assessor*product interaction MS, in case
  # a  scaling difference is not detected (with alpha_conditionalMAM)
  # And if the adjusted MAM is used the negativity adjusted MS(dis) will
  # be used (cf. Brockhoff et al, 2013)
  # result[[3]][,,1] contains the ANOVA table for attribute 1. 
  #
  #
  # result[[4]]:
  # A 3-way array with the MAM based post hoc analysis of product differences
  # for each attribute using the adjusted/conditional MAM as controlled by   
  # the two options. The product differences are shown together with the  
  # post hoc p-value. 
  # result[[4]][,,1] contains the post hoc analysis for attribute 1.
  #
  #
  # result[[5]]:
  # A matrix with the MAM based post hoc comparison of each product 
  # with the mean of the remaining products for each attribute 
  # using the adjusted/conditional MAM as controlled by   
  # the two options. The product differences are shown together with the  
  # post hoc p-values categorized the usual way. 
  #
  #
  # result[[6]]:
  # A 3-way array with the MAM based confidence limits for the product differences
  # for each attribute using the new method introduced in Brockhoff et al (2013)
  # The lower and upper 95% confidence limits are shown for each comparison
  # result[[6]][,,1] contains the Confidence limits for attribute 1.
  #
  #
  #  
attnames       =  labels(data)[[2]][-1:-3]                    ; natt  = length(attnames)                                                          # Save the attribute names and number.
ass            =  factor(data[,1]) ; assnames  = levels(ass)  ; nass  = length(assnames)                                                          # Define as factor, save level names (alphanum.) and #.
prod           =  factor(data[,2]) ; prodnames = levels(prod) ; nprod = length(prodnames)                                                         #                          -"-
rep            =  factor(data[,3]) ; repnames  = levels(rep)  ; nrep  = length(repnames)                                                          #                          -"-
nrow           =  dim(data)[1]                                ; ncol  = dim(data)[2]                                                             # Number of rows and selected columns in the data data.
attnames1      =  attnames2 = attnames ; 
asslevels=unique(ass)

ass_nr=rep(0,nrow)
for (i in 1:nass) ass_nr[ass==assnames[i]]=i

data[,1]=factor(ass_nr)

data[,2]=prod
data[,3]=rep
names(data)[1]="ass"
names(data)[2]="prod"
names(data)[3]="rep"

data=orderBy(~ass+prod,data)

if (natt>1) attnames1[seq(2,natt,2)] = '' ; 
attnames2[seq(1,natt,2)] = ''                                           
  # Splitting attnames in two to allow for longer names.

prodFs         =  prodPs = SEs = rep(0,natt)                                                                                                       # Initializing p-value, F-value and SE vectors.
prodco         =  matrix(rep(0,natt*nprod),nrow=natt)                                                                                              # Initialize a natt x nprod matrix for product coefficients. 
const          =  rep(1,nrow)                                                                                                                      # Vector with ones for the design matrix (maybe rubbish!).

indivTable=array(dim=c(7,nass,natt))
perfTable=array(dim=c(14,nass,natt))

IPTable=array(dim=c(6,nass+1,natt))
PIPTable=array(dim=c(6,nass,natt))

FinalTable=array(dim=c(12,nass+1,natt))
AOVtable=array(dim=c(5,5,natt))

prodmeans=array(dim=c(nprod,natt))

assmeans=array(dim=c(nass,natt))
proddiffs=array(dim=c(nprod*(nprod-1)/2,natt))

colnames(proddiffs) <- attnames ## added by alku

stddiffs=array(dim=c(nprod*(nprod-1)/2,natt))

colnames(stddiffs) <- attnames ## added by alku

#proddiffs=array(dim=c(nprod,nprod,natt))
prodeffects=array(dim=c(natt,nprod))

#print(c(nrow,nprod*nass*nrep))
mins=apply(data[,4:ncol],2,min)
maxs=apply(data[,4:ncol],2,max)
ranges=maxs-mins

nnegatives=rep(0,natt)
for (p in 4:ncol){ # running through attributes one by one
X=data[,p]
mu=mean(X)
data$X=X
xam=summaryBy(X~ass,data=data)
xpm=summaryBy(X~prod,data=data)
prodmeans[,p-3]=xpm$X.mean

rownames(prodmeans) <- xpm$prod ## added by alku
colnames(prodmeans) <- attnames ## added by alku

assmeans[,p-3]=xam$X.mean
xapm=summaryBy(X~ass+prod,data=data)
xapm$int=xapm$X.mean-rep(xam$X.mean,rep(nprod,nass))-rep(xpm$X.mean,nass)+mu
ssa=nprod*nrep*(xam$X.mean-mu)^2
#ssp=rep(nrep*sum((xpm$X.mean-mu)^2),nass)
ssp=rep(0,nass)
sse=rep(0,nass)
sssca=rep(0,nass)
ssdis=rep(0,nass)
beta=rep(0,nass)
xapm$xs=rep(xpm$X.mean,nass)
sses=X-ave(X,factor((data$ass):(data$prod)))

for (i in 1:nass){
ssp[i]=anova(lm(X~prod,data=subset(data,ass==i)))[1,2]
res=lm(int~xs,data=subset(xapm,ass==i))
aovres=anova(res)[,2]
sssca[i]=nrep*aovres[1]
ssdis[i]=nrep*aovres[2]
ssp[i]=ssp[i]-sssca[i]-ssdis[i]
beta[i]= res$coef[2]
sse[i]=sum(sses[data$ass==i]^2)
}

indivTable[1,,p-3]=ssa
indivTable[2,,p-3]=ssp
indivTable[3,,p-3]=sssca
indivTable[4,,p-3]=ssdis
indivTable[5,,p-3]=sse
indivTable[6,,p-3]=beta+1
indivTable[7,,p-3]=scale(assmeans[,p-3],scale=F)


k=0
names.diff <- NULL
for  (i in 1:(nprod-1)) for (j in (i+1):nprod){
k <- k+1
names.diff[k] <- paste(rownames(prodmeans)[i], "-", rownames(prodmeans)[j])
proddiffs[k,p-3] <- (prodmeans[i,p-3]-prodmeans[j,p-3])
}

rownames(stddiffs) <- rownames(proddiffs) <- names.diff


for  (i in 1:nprod) prodeffects[p-3,i]=prodmeans[i,p-3]-mean(prodmeans[-i,p-3])

}


aovs=apply(indivTable,c(1,3),sum)[1:5,]

msa=aovs[1,]/(nass-1)
msc=aovs[3,]/(nass-1)
msd=aovs[4,]/((nass-1)*(nprod-2))
msp=aovs[2,]/(nprod-1)
mse=aovs[5,]/(nass*nprod*(nrep-1))
msint=(aovs[3,]+aovs[4,])/((nprod-1)*(nass-1))

for (p in 4:ncol){ # running through attributes one by one
perfTable[1,,p-3]=nass*indivTable[1,,p-3]/(msint[p-3]*(nass-1))
perfTable[2,,p-3]=1-pf(perfTable[1,,p-3],1,(nass-1)*(nprod-1))
perfTable[3,,p-3]=scale(assmeans[,p-3],scale=F)
perfTable[4,,p-3]=((indivTable[2,,p-3]+indivTable[3,,p-3]+indivTable[4,,p-3])/(nprod-1))/(indivTable[5,,p-3]/(nprod*(nrep-1)))
perfTable[5,,p-3]=1-pf(perfTable[4,,p-3],nprod-1,nprod*(nrep-1))
perfTable[6,,p-3]=indivTable[3,,p-3]/(indivTable[4,,p-3]/(nprod-2))
perfTable[7,,p-3]=1-pf(perfTable[6,,p-3],1,(nprod-2))
perfTable[8,,p-3]=(indivTable[6,,p-3]^2)*(aovs[2,p-3]/(nass))/(indivTable[4,,p-3]/(nprod-2))
perfTable[9,,p-3]=1-pf(perfTable[8,,p-3],1,(nprod-2))
perfTable[10,,p-3]=indivTable[6,,p-3]*sqrt((aovs[2,p-3]/nass)/(indivTable[2,,p-3]+indivTable[3,,p-3]+indivTable[4,,p-3]))
perfTable[11,,p-3]=(indivTable[4,,p-3]/(nprod-2))/(indivTable[5,,p-3]/(nprod*(nrep-1)))
perfTable[12,,p-3]=1-pf(perfTable[11,,p-3],nprod-2,nprod*(nrep-1))
for (i in 1:nass) perfTable[13,i,p-3]=indivTable[5,i,p-3]/mean(indivTable[5,-i,p-3])
perfTable[14,,p-3]=1-pf(perfTable[13,,p-3],nprod*(nrep-1),(nass-1)*nprod*(nrep-1))

IPTable[1,1:nass,p-3]=scale(assmeans[,p-3],scale=F)
#IPTable[2,1:nass,p-3]=sqrt(((indivTable[2,,p-3]+indivTable[3,,p-3]+indivTable[4,,p-3])/(nrep*(nprod-1))))
IPTable[2,1:nass,p-3]=sqrt(perfTable[4,,p-3])
IPTable[3,1:nass,p-3]=indivTable[6,,p-3]
IPTable[4,1:nass,p-3]=indivTable[6,,p-3]*sqrt((aovs[2,p-3]/nass)/(indivTable[2,,p-3]+indivTable[3,,p-3]+indivTable[4,,p-3]))
IPTable[5,1:nass,p-3]=sqrt(indivTable[4,,p-3]/(nprod-2))
IPTable[6,1:nass,p-3]=sqrt(indivTable[5,,p-3]/(nprod*(nrep-1)))
IPTable[,nass+1,p-3]=apply(IPTable[,1:nass,p-3],1,mean)

PIPTable[1,,p-3]=perfTable[2,,p-3]
PIPTable[2,,p-3]=perfTable[5,,p-3]
PIPTable[3,,p-3]=perfTable[7,,p-3]
PIPTable[4,,p-3]=perfTable[9,,p-3]
PIPTable[5,,p-3]=perfTable[12,,p-3]
PIPTable[6,,p-3]=perfTable[14,,p-3]

IPtabT=matrix(as.character(round(IPTable[,,p-3],2)),nrow=6)
PIPtabT=IPtabT
for (i in 1:6){
PIPtabT[i,1:nass]=" "
PIPtabT[i,1:nass][PIPTable[i,,p-3]<=0.10]="."
PIPtabT[i,1:nass][PIPTable[i,,p-3]<=0.05]="*"
PIPtabT[i,1:nass][PIPTable[i,,p-3]<=0.01]="**"
PIPtabT[i,1:nass][PIPTable[i,,p-3]<=0.001]="***"
}

FinalTable[,,p-3]=
matrix(c(IPtabT[1,],PIPtabT[1,],IPtabT[2,],PIPtabT[2,],IPtabT[3,],PIPtabT[3,],IPtabT[4,],
PIPtabT[4,],IPtabT[5,],PIPtabT[5,],IPtabT[6,],PIPtabT[6,]),ncol=nass+1,byrow=T)

for (i in c(2,4,6,8,10,12)) FinalTable[i,nass+1,p-3]=""
}

Fmatr=matrix(rep(0,7*natt),nrow=7)
Pmatr=matrix(rep(0,7*natt),nrow=7)

if (adjustedMAM) {
  negatives=(indivTable[6,,]<0)
  nnegatives=apply(negatives,2,sum)
  for (p in 1:natt){
    aovs[3,p]=aovs[3,p]-sum(indivTable[3,,p][negatives[,p]])
    aovs[4,p]=aovs[4,p]+sum(indivTable[3,,p][negatives[,p]])
    msc[p]=(aovs[3,p])/(nass-1-nnegatives[p])
    msd[p]=(aovs[4,p])/((nass-1)*(nprod-2)+nnegatives[p])  
  }
}

Fmatr[1,]=msa/msint
Fmatr[2,]=msp/msd
Fmatr[3,]=msc/msd
Fmatr[4,]=msd/mse

Fmatr[5,]=msp/msint
Fmatr[6,]=msint/mse
Fmatr[7,]=((aovs[2,]+aovs[3,])/(nprod+nass-2))/msd

DFnums=c(nass-1,nprod-1,nass-1,(nprod-2)*(nass-1),nprod-1,(nprod-1)*(nass-1),(nprod+nass-2))

DFdens=c((nprod-1)*(nass-1),(nprod-2)*(nass-1),(nprod-2)*(nass-1),nass*nprod*(nrep-1),
         (nprod-1)*(nass-1),nass*nprod*(nrep-1),(nprod-2)*(nass-1))

for (i in 1:7) Pmatr[i,]=1-pf(Fmatr[i,],DFnums[i],DFdens[i])


AOVtable[,1,]=round(aovs,2)

if (adjustedMAM) {
  for (p in 1:natt) {
    DFnums=c(nass-1,nprod-1,nass-1-nnegatives[p],(nprod-2)*(nass-1)+nnegatives[p],nprod-1,(nprod-1)*(nass-1),(nprod+nass-2-nnegatives[p]))
    AOVtable[,3,p]=c(DFnums[1:4],nass*nprod*(nrep-1))
    DFdens=c((nprod-1)*(nass-1),(nprod-2)*(nass-1)+nnegatives[p],(nprod-2)*(nass-1)+nnegatives[p],nass*nprod*(nrep-1),
             (nprod-1)*(nass-1),nass*nprod*(nrep-1),(nprod-2)*(nass-1)+nnegatives[p])
    for (i in 1:7) {
      Pmatr[i,p]=1-pf(Fmatr[i,p],DFnums[i],DFdens[i])
    }
    AOVtable[,2,p]=round(aovs[,p]/AOVtable[,3,p],2) 
  }
}

if (!adjustedMAM) {
  dfs=c(nass-1,nprod-1,nass-1,(nprod-2)*(nass-1),nass*nprod*(nrep-1))
  AOVtable[,3,]=dfs
  AOVtable[,2,]=round(aovs/dfs,2)
}

AOVtable[1:4,4,]=round(Fmatr[1:4,],2)
AOVtable[1:4,5,]=round(Pmatr[1:4,],4)


# Here the MAM condition is applied, such that the product F may become the usual one again 
MAMcondition=(Pmatr[3,]<=alpha_conditionalMAM)
AOVtable[2,4,][!MAMcondition]=round(Fmatr[5,][!MAMcondition],2) 
AOVtable[2,5,][!MAMcondition]=round(Pmatr[5,][!MAMcondition],2) 

dimnames(FinalTable)[[1]]=c("Level","","Product"," ","Scaling", 
"  ","Correlation","   ","Disagreement","    ","Repeatability","      ")
#dimnames(FinalTable)[[2]]=c(paste(rep("Ass",nass),1:nass),"AVE")
dimnames(FinalTable)[[2]]=c(assnames,"AVE")
dimnames(FinalTable)[[3]]=attnames

dimnames(AOVtable)[[1]]=c("Assessor","Product","Scaling","Disagreement","Error")
dimnames(AOVtable)[[2]]=c("SS","MS","DF","F","Pval")
dimnames(AOVtable)[[3]]=attnames

# First the int-based p-values are found:
stddiffs[] <- matrix(rep(sqrt(2*msint/(nrep*nass)),nprod*(nprod-1)/2),
                   ncol=natt,byrow=T)

diffpvals=1-pt(abs(proddiffs)/matrix(rep(sqrt(2*msint/(nrep*nass)),nprod*(nprod-1)/2),
                                     ncol=natt,byrow=T),(nprod-1)*(nass-1))
effectspvals=1-pt(abs(prodeffects)/matrix(rep(sqrt(msint*(1/(nrep*nass)+1/(nrep*nass*(nprod-1)))),
rep(nprod,natt)),ncol=nprod,byrow=T),(nprod-1)*(nass-1))


for (p in 1:natt){
if (MAMcondition[p]) {
stddiffs[, p] <- matrix(rep(sqrt(2*msd/(nrep*nass)),nprod*(nprod-1)/2),
         ncol=natt,byrow=T)[, p]  
diffpvals[,p]=1-pt(abs(proddiffs[,p])/matrix(rep(sqrt(2*msd/(nrep*nass)),nprod*(nprod-1)/2),
                                       ncol=natt,byrow=T)[,p],(nprod-2)*(nass-1))
effectspvals[p,]=1-pt(abs(prodeffects[p,])/matrix(rep(sqrt(msd*(1/(nrep*nass)+1/(nrep*nass*(nprod-1)))),
                                              rep(nprod,natt)),ncol=nprod,byrow=T)[p,],(nprod-2)*(nass-1))
if (adjustedMAM) {
  dfmsds=matrix(rep(AOVtable[4,3,],nprod*(nprod-1)/2),ncol=natt,byrow=T)
  diffpvals[,p]=1-pt(abs(proddiffs[,p])/matrix(rep(sqrt(2*msd/(nrep*nass)),nprod*(nprod-1)/2),ncol=natt,byrow=T)[,p],dfmsds[p])
effectspvals[p,]=1-pt(abs(prodeffects[p,])/matrix(rep(sqrt(msd*(1/(nrep*nass)+1/(nrep*nass*(nprod-1)))),
                                              rep(nprod,natt)),ncol=nprod,byrow=T)[p,],dfmsds[p])
}
}
}
diffpvals=2*diffpvals
effectspvals=2*effectspvals

PtabT=matrix(as.character(round(prodeffects,2)),ncol=nprod)
PPtabT=effectspvals
for (i in 1:natt){
PPtabT[i,]=" "
PPtabT[i,][effectspvals[i,]<=0.10]="."
PPtabT[i,][effectspvals[i,]<=0.05]="*"
PPtabT[i,][effectspvals[i,]<=0.01]="**"
PPtabT[i,][effectspvals[i,]<=0.001]="***"
}

ProdTable=matrix(rep(" ",2*nprod*natt),ncol=nprod)

rownames(ProdTable)=c(rbind(attnames,rep(" ",natt)))
colnames(ProdTable)=prodnames

i1=(1:natt)*2-1
i2=(1:natt)*2
ProdTable[i1,]=PtabT
ProdTable[i2,]=PPtabT

i1=(1:nprod)*2-1
i2=(1:nprod)*2
DIFtable=array(dim=c(2*(nprod-1),nprod-1,natt))

## create difflsmeans table for each attribute
## similar to lmerTest, which also contains std errors
lsmeans.tab <- array(dim = c(nrow(proddiffs), 3 , natt)) ## added by alku


dimnames(lsmeans.tab)[[1]] <- rownames(proddiffs)
dimnames(lsmeans.tab)[[2]] <- c("Estimate", "Standard Error", "Pval")
dimnames(lsmeans.tab)[[3]] <- attnames

for(k in 1:nrow(proddiffs)){
  lsmeans.tab[k,1,] <-  round(proddiffs[k,], 3)
  lsmeans.tab[k,2,] <-  round(stddiffs[k,], 3)
  lsmeans.tab[k,3,] <-  round(diffpvals[k,], 4)
}


k=0
for  (i in 1:(nprod-1)) for (j in (i+1):nprod){
k=k+1
DIFtable[2*i-1,j-1,]=round(proddiffs[k,],2)
DIFtable[2*i,j-1,]=round(diffpvals[k,],4)
}

dimnames(DIFtable)[[1]]=c(rbind(prodnames[-nprod],rep(" ",nprod-1)))
dimnames(DIFtable)[[2]]=prodnames[-1]
dimnames(DIFtable)[[3]]=attnames

IndividualTable=indivTable[1:5,,]

for (i in 1:nass) IndividualTable[2,i,]=AOVtable[2,2,]/nass
dimnames(IndividualTable)[[1]]=c("Assessor","Product","Scaling","Disagreement","Error")
dimnames(IndividualTable)[[2]]=assnames
dimnames(IndividualTable)[[3]]=attnames

#Change row-orders in FinalTable:
FN=FinalTable
FinalTable=FN[c(3,4,5,6,9,10,11,12,1,2,7,8),,]


# Confidence limits:

# First the int-based:
  LLs=proddiffs-
  matrix(rep(qt(0.975,(nprod-1)*(nass-1))*sqrt(2*msint/(nrep*nass)),nprod*(nprod-1)/2),ncol=natt,byrow=T)
ULs=proddiffs+
  matrix(rep(qt(0.975,(nprod-1)*(nass-1))*sqrt(2*msint/(nrep*nass)),nprod*(nprod-1)/2),ncol=natt,byrow=T)

# Next the ms_dis based (to be used in the computation of the correct MAM ones) is found
# These are IN CASE of adjustedMAM=TRUE adjusted BACK to become the UNadjusted MAM ones
# The CI method in Brockhoff et al(2013) only covers this case

# The potential "unadjustment":
if (adjustedMAM) {
  negatives=(indivTable[6,,]<0)
  for (p in 1:natt){
    aovs[3,p]=aovs[3,p]+sum(indivTable[3,,p][negatives[,p]])
    aovs[4,p]=aovs[4,p]-sum(indivTable[3,,p][negatives[,p]])
    msc[p]=(aovs[3,p])/(nass-1)
    msd[p]=(aovs[4,p])/((nass-1)*(nprod-2))  
  }
}

LLs_dis=proddiffs-
  matrix(rep(qt(0.975,(nprod-2)*(nass-1))*sqrt(2*msd/(nrep*nass)),nprod*(nprod-1)/2),ncol=natt,byrow=T)
ULs_dis=proddiffs+
  matrix(rep(qt(0.975,(nprod-2)*(nass-1))*sqrt(2*msd/(nrep*nass)),nprod*(nprod-1)/2),ncol=natt,byrow=T)


for (j in 1:(ncol-3)) {
    if (MAMcondition[j]) {
      
    m=0
    if (Fmatr[7,j]>1)  for (i in 1:(nprod-1)) for (n in (i+1):nprod)   {
      m=m+1
      
      Uprob=function(x){
        prodmeansnow=prodmeans[,j]
        prodmeansnow[i]=prodmeans[i,j]+x/2
        prodmeansnow[n]=prodmeans[n,j]-x/2
        k=min(0.5*(proddiffs[m,j]+x)^2/(sum((prodmeansnow-mean(prodmeansnow))^2)),1)
        se=sqrt(2*(msc[j]*k+(1-k)*msd[j])/(nrep*nass))
        df=((k*msc[j]+(1-k)*msd[j])^2)/(((k*msc[j])^2)/(nass-1)+(((1-k)*msd[j])^2)/((nass-1)*(nprod-2)))
        1-pt(x/se,df)-0.025
      }
      
      Lprob=function(x){
        prodmeansnow=prodmeans[,j]
        prodmeansnow[i]=prodmeans[i,j]-x/2
        prodmeansnow[n]=prodmeans[n,j]+x/2
        k=min(0.5*(proddiffs[m,j]-x)^2/(sum((prodmeansnow-mean(prodmeansnow))^2)),1)
        se=sqrt(2*(msc[j]*k+(1-k)*msd[j])/(nrep*nass))
        df=((k*msc[j]+(1-k)*msd[j])^2)/(((k*msc[j])^2)/(nass-1)+(((1-k)*msd[j])^2)/((nass-1)*(nprod-2)))
        1-pt(x/se,df)-0.025
      }
      
      if (ULs_dis[m,j]<0) tryres <- try(uniroot(Uprob,c(0,abs(proddiffs[m,j]))),silent = TRUE)
      if (ULs_dis[m,j]<0){ if (class(tryres) == "try-error") ULs[m,j]=NA }
      if (ULs_dis[m,j]<0){ if (class(tryres) != "try-error") ULs[m,j]=proddiffs[m,j]+uniroot(Uprob,c(0,abs(proddiffs[m,j])))$root}
      if (ULs_dis[m,j]>=0){
        x=1
        while ((Uprob(x)> 0)& (x<ranges[j]+1)) x=x+1  
        if (x<ranges[j]+1) ULs[m,j]=proddiffs[m,j]+uniroot(Uprob,c(0,x))$root else ULs[m,j]=proddiffs[m,j]+ranges[j]+1
      }
      
      if (LLs_dis[m,j]>0) tryres <- try(uniroot(Lprob,c(0,abs(proddiffs[m,j]))),silent = TRUE)
      if (LLs_dis[m,j]>0){ if (class(tryres) == "try-error") LLs[m,j]=NA }
      if (LLs_dis[m,j]>0){ if (class(tryres) != "try-error") LLs[m,j]=proddiffs[m,j]-uniroot(Lprob,c(0,abs(proddiffs[m,j])))$root}
      if (LLs_dis[m,j]<=0){
        x=1
        while ((Lprob(x)> 0)& (x<ranges[j]+1)) x=x+1  
        if (x<ranges[j]+1) LLs[m,j]=proddiffs[m,j]-uniroot(Lprob,c(0,x))$root else LLs[m,j]=proddiffs[m,j]-(ranges[j]+1)
      } 
  } # for i and n
} # MAMcondition
} # for j

CLTable=DIFtable
k=0
for  (i in 1:(nprod-1)) for (j in (i+1):nprod){
  k=k+1
  CLTable[2*i-1,j-1,]=round(LLs[k,],2)
  CLTable[2*i,j-1,]=round(ULs[k,],4)
}


## create table for CI for the shiny app in SensMixed
ci.tab <- array(dim = c(nrow(proddiffs), 2 , natt)) ## added by alku


dimnames(ci.tab)[[1]] <- rownames(proddiffs)
dimnames(ci.tab)[[2]] <- c("Lower CI", "Upper CI")
dimnames(ci.tab)[[3]] <- attnames

for(k in 1:nrow(proddiffs)){
  ci.tab[k,1,] <-  round(LLs[k,], 4)
  ci.tab[k,2,] <-  round(ULs[k,], 4)
}


list(IndividualTable, FinalTable, AOVtable, DIFtable, lsmeans.tab, ProdTable,
     CLTable, ci.tab)

}

