#
# loading the data
#
# first, you can convert data from ped-files to internal format"
adpr <- paste(.Library,"/GenABEL/exdata/",sep="")
genofile <- paste(adpr,"pedin.18",sep="")
cat("location of ped-file:",genofile,"\n")
mapfile <- paste(adpr,"map.18",sep="")
cat("location of map-file:",mapfile,"\n")
phfile <- paste(adpr,"phenos.18",sep="")
cat("location of pheno-file:",phfile,"\n")
convert.snp.ped(pedfile=genofile,mapfile=mapfile,outfile="chr18.raw",bcast=50)
a<-readline("PRESS <Enter> TO CONTINUE...")
# These file can be read in using load.gwaa.data:"
srdta <- load.gwaa.data(geno="chr18.raw",pheno=phfile)
a<-readline("PRESS <Enter> TO CONTINUE...")
# Now, let us convert and read data from other human-readable format:"
genofile <- paste(adpr,"srgenos.dat",sep="")
cat("location of file with human-readable genotypes:",genofile,"\n")
a<-readline("PRESS <Enter> TO CONTINUE...")
convert.snp.text(infile=genofile,outfile="srgenos.raw",bcast=200)
a<-readline("PRESS <Enter> TO CONTINUE...")
phenofile <- paste(adpr,"srphenos.dat",sep="")
cat("location of file with phenotypes:",phenofile,"\n")
srdta <- load.gwaa.data(geno="srgenos.raw",pheno=phenofile)
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# manipulating phenotypic data
#
names(srdta@phdata)
summary(srdta@phdata)
a<-readline("PRESS <Enter> TO CONTINUE...")
hist(srdta@phdata$age)
a<-readline("PRESS <Enter> TO CONTINUE...")
summary(glm(srdta@phdata$qt2 ~ srdta@phdata$age + srdta@phdata$sex))
a<-readline("PRESS <Enter> TO CONTINUE...")
attach(srdta@phdata)
summary(glm(qt2 ~ age + sex))
a<-readline("PRESS <Enter> TO CONTINUE...")
summary(glm(bt ~ age + sex, family=binomial()))
detach(srdta@phdata)
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# phenotypic QC
#
check.trait("age",srdta)
a<-readline("PRESS <Enter> TO CONTINUE...")
check.trait("qt1",srdta)
a<-readline("PRESS <Enter> TO CONTINUE...")
alltra <- names(srdta@phdata)
alltra
check.trait(alltra,srdta@phdata[srdta@phdata$bt==0,],fdrate=0.2,graph=FALSE)
a<-readline("PRESS <Enter> TO CONTINUE...")
hist(srdta@phdata$qt2)
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta@phdata$qt2 <- replace(srdta@phdata$qt2,(srdta@phdata$qt2==888),NA)
check.trait("qt2",srdta)
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(srdta@phdata$age,srdta@phdata$qt2)
a<-readline("PRESS <Enter> TO CONTINUE...")
summary(glm(srdta@phdata$qt2 ~ srdta@phdata$age + srdta@phdata$sex))
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# Manipulating genetic data
#
srdta@gtdata@nids
srdta@gtdata@nsnps
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta@gtdata@idnames[1:12]
srdta@gtdata@male[1:12]
srdta@gtdata@male[c("p1","p2","p3","p4")]
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta@gtdata@snpnames[1:4]
srdta@gtdata@chromosome[1:4]
srdta@gtdata@map[1:4]
a<-readline("PRESS <Enter> TO CONTINUE...")
n4 <- c("rs18","rs655")
n4
srdta@gtdata@map[n4]
n4 <- c("rs18","rs65")
n4
srdta@gtdata@map[n4]
srdta@gtdata@chromosome[n4]
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta[1:12,1:4]
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta@phdata[1:12,]
srdta@gtdata[1:12,1:4]
a<-readline("PRESS <Enter> TO CONTINUE...")
as.numeric(srdta@gtdata[1:12,1:4])
a<-readline("PRESS <Enter> TO CONTINUE...")
as.numeric(srdta@gtdata[c("p1","p3","p4"),c("rs18","rs65")])
a<-readline("PRESS <Enter> TO CONTINUE...")
as.character(srdta@gtdata[c("p1","p3","p4"),c("rs18","rs65")])
a<-readline("PRESS <Enter> TO CONTINUE...")
as.genotype(srdta@gtdata[c("p1","p3","p4"),c("rs18","rs65")])
a<-readline("PRESS <Enter> TO CONTINUE...")
as.hsgeno(srdta@gtdata[c("p1","p3","p4"),c("rs18","rs65")])
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta@gtdata@snpnames[1:20]
srdta@gtdata@map[1:20]
snp.names(srdta,end=31000)
nams <- snp.names(srdta,end=31000)
nams
srdta@gtdata@map[nams]
srdta@gtdata@chromosome[nams]
a<-readline("PRESS <Enter> TO CONTINUE...")
snp.names(srdta,begin=5000,end=31000)
a<-readline("PRESS <Enter> TO CONTINUE...")
snp.names(srdta,end="rs143")
snp.names(srdta,begin="rs29",end="rs143")
a<-readline("PRESS <Enter> TO CONTINUE...")
snp.names(srdta,begin="rs29",end="rs143")
maprs114 <- srdta@gtdata@map["rs114"]
maprs114
snp.names(srdta,begin=maprs114-12000,end=maprs114+12000)
snp.names(srdta,begin=maprs114-12000,end="rs130")
a<-readline("PRESS <Enter> TO CONTINUE...")
snn <- snp.names(srdta,begin=maprs114-12000,end="rs130")
snn
ids <- srdta@gtdata@idnames[1:4]
ids
as.numeric(srdta@gtdata[ids,snn])
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# genetic data QC
#
# Not run: very long output
# summary(srdta@gtdata)
a<-readline("PRESS <Enter> TO CONTINUE...")
summary(srdta@gtdata[,1:10])
a<-readline("PRESS <Enter> TO CONTINUE...")
summary(srdta@gtdata[srdta@phdata$bt==0,1:10])
a<-readline("PRESS <Enter> TO CONTINUE...")
mc <- check.marker(srdta,ibs.mrk=200)
#
# generate summary per SNP and per ID
#
a<-readline("PRESS <Enter> TO CONTINUE...")
summary(mc)
#
# show HWE details for markers which failed HWE test
#
a<-readline("PRESS <Enter> TO CONTINUE...")
HWE.show(snps=mc$nohwe,data=srdta)
#
# plot mc
#
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(mc)
#
# other round of QC with user threshold values
#
a<-readline("PRESS <Enter> TO CONTINUE...")
mc <- check.marker(data=srdta,call=0.94,maf=0.01,p.level=0.01,ibs.exclude="both")
summary(mc)
plot(mc)
a<-readline("PRESS <Enter> TO CONTINUE...")
HWE.show(snps=mc$nohwe,data=srdta)
#
# use only "bt" controls for HWE checks
#
a<-readline("PRESS <Enter> TO CONTINUE...")
mc <- check.marker(data=srdta,call=0.94,maf=0.01,p.level=0.01,hweids=(srdta@phdata$bt==0),ibs.exclude="both")
summary(mc)
plot(mc)
a<-readline("PRESS <Enter> TO CONTINUE...")
mc <- check.marker(data=srdta,call=0.94,maf=0.01,p.level=0.01,mincon=.90,ibs.exclude="both")
summary(mc)
plot(mc)
a<-readline("PRESS <Enter> TO CONTINUE...")
mc1 <- snp.subset(mc,snp.names(srdta,begin=300000,end=500000,chrom="1"))
plot(mc1)
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# running analysis
#
mc <- check.marker(data=srdta,call=0.94,maf=(5/srdta@gtdata@nids),fdr=0.05,hweids=(srdta@phdata$bt==0),ibs.exclude="both")
summary(mc)
mymrk <- mc$snpok
#
# case-control analysis (not recommended, rather use qtscore)
#
a<-readline("PRESS <Enter> TO CONTINUE...")
res.fcc <- ccfast("bt",data=srdta,snps=mymrk)
#
# what is the minimal P-value obsreved?
#
min(res.fcc$P1df)
-log10(min(res.fcc$P1df))
#
# at which SNP the P-value <= 0.01 observed?
#
b<-which(res.fcc$P1df<=0.01)
b
res.fcc$snpnames[b]
res.fcc$P1df[b]
#
# plot the results
#
plot(res.fcc)
#
# estimate experimet-wise significance
#
a<-readline("PRESS <Enter> TO CONTINUE...")
res.boocc <- emp.ccfast("bt",data=srdta,snps=mymrk)
min(res.boocc$P1df)
plot(res.boocc)
#
# run simular analyses qith qtscore
#
a<-readline("PRESS <Enter> TO CONTINUE...")
res.score <- qtscore(bt,data=srdta,snps=mymrk,trait="binomial")
min(res.score$P1df)
-log10(min(res.score$P1df))
which(res.score$P1df<=0.01)
plot(res.score)
#
# overlay ccfast results over qtscore
#
a<-readline("PRESS <Enter> TO CONTINUE...")
add.plot(res.fcc,col="red")
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(res.fcc$P1df,res.score$P1df,pch=20)
a<-readline("PRESS <Enter> TO CONTINUE...")
res.score <- qtscore(bt ~ sex + age,data=srdta,snps=mymrk)
min(res.score$P1df)
-log10(min(res.score$P1df))
which(res.score$P1df<=0.01)
plot(res.score)
a<-readline("PRESS <Enter> TO CONTINUE...")
res.boosc <- emp.qtscore(srdta@phdata$bt,data=srdta,snps=mymrk)
min(res.boosc$P1df)
plot(res.boosc)
a<-readline("PRESS <Enter> TO CONTINUE...")

# analysis of a region

index <- which.min(res.score$P1df)
index
osnp <- res.score$snpname[index]
osnp
pososnp <- srdta@gtdata@map[osnp]
pososnp
pososnp <- res.score$map[osnp]
pososnp
# wrong -- it will get SNP names from all data, including these which did not pass QC!
reg <- snp.names(res.score,begin=pososnp-50000,end=pososnp+50000) 
reg
r1.fcc <- snp.subset(res.fcc,snps=reg)
min(r1.fcc$P1df)
r1.score <- snp.subset(res.score,snps=reg)
min(r1.score$P1df)
r1.1.glm <- scan.glm("bt ~ sex + age + CRSNP",family=binomial(),data=srdta,snps=reg,bcast=20)
min(r1.1.glm$P1df)
min(r1.1.glm$P2df)
ymax <- -log10(min(r1.1.glm$P1df))+.1
ymax
plot(r1.fcc,ylim=c(0,ymax))
add.plot(r1.score,df=1,col="red",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(r1.fcc,ylim=c(0,ymax))
points(r1.score$map,-log10(r1.score$P1df),col="red",type="l")
points(r1.1.glm$map,-log10(r1.1.glm$P1df),col="green",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(r1.fcc,ylim=c(0,ymax))
points(r1.score$map,-log10(r1.score$P1df),col="red",type="l")
add.plot(r1.1.glm,df=1,col="green",type="l")
add.plot(r1.1.glm,df=2,col="blue",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")

g2d <- scan.glm.2D("bt ~ sex + age + CRSNP",family=binomial(),data=srdta,snps=reg[10:30],bcast=20)
plot(g2d)
a<-readline("PRESS <Enter> TO CONTINUE...")

# haplotypic analysis for bt

library(haplo.stats)
x.adj <- matrix(c(srdta@phdata$sex,srdta@phdata$age),ncol=2)
srdta@phdata[1:5,]
x.adj[1:5,]
a<-readline("PRESS <Enter> TO CONTINUE...")
r1.2.hap <- scan.haplo("bt~sex+age+CRSNP",data=srdta,snps=reg,trait="binomial",n.slide=2)
min(r1.2.hap$P1df)
r1.3.hap <- scan.haplo("bt~sex+age+CRSNP",data=srdta,snps=reg,trait="binomial",n.slide=3)
min(r1.3.hap$P1df)
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(r1.1.glm)
points(r1.2.hap$map,-log10(r1.2.hap$P1df),col="red",type="l")
points(r1.3.hap$map,-log10(r1.3.hap$P1df),col="green",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")
sr1 <- snp.names(r1.1.glm,beg=2150000,end=2170000,chrom="1")
plot(snp.subset(r1.1.glm,snp=sr1))
points(r1.2.hap$map,-log10(r1.2.hap$P1df),col="red",type="l")
points(r1.3.hap$map,-log10(r1.3.hap$P1df),col="green",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")
sr1 <- snp.names(r1.1.glm,beg=2175000,end=2200000,chrom="1")
plot(snp.subset(r1.1.glm,snp=sr1))
points(r1.2.hap$map,-log10(r1.2.hap$P1df),col="red",type="l")
points(r1.3.hap$map,-log10(r1.3.hap$P1df),col="green",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")

h2d <- scan.haplo.2D("bt~sex+age+CRSNP",data=srdta,snps=reg[10:30],trait="binomial")
plot(h2d)
a<-readline("PRESS <Enter> TO CONTINUE...")


# analysis of trait qt2

attach(srdta@phdata)
check.trait("qt2",data=srdta)
summary(glm(qt2~age+sex))
a<-readline("PRESS <Enter> TO CONTINUE...")
r.qts <- qtscore(srdta@phdata$qt2,data=srdta,snps=mymrk)
min(r.qts$P1df)
plot(r.qts)
a<-readline("PRESS <Enter> TO CONTINUE...")
r.qts <- qtscore(srdta@phdata$qt2,strata=srdta@phdata$sex,data=srdta,snps=mymrk)
min(r.qts$P1df)
plot(r.qts)
a<-readline("PRESS <Enter> TO CONTINUE...")
r.aqts <- qtscore(qt2~age+sex,data=srdta,snps=mymrk)
min(r.aqts$P1df)
which.min(r.aqts$P1df)
plot(r.aqts)
add.plot(r.aqts,col="red",cex=2.)
a<-readline("PRESS <Enter> TO CONTINUE...")
r.booqts <- emp.qtscore(qt2~age+sex,data=srdta,snps=mymrk,times=200)
min(r.booqts$P1df)
which.min(r.booqts$P1df)
plot(r.booqts)
a<-readline("PRESS <Enter> TO CONTINUE...")

#investigation of the region
osnp <- r.booqts$snpname[which.min(r.booqts$P1df)]
osnp
pososnp <- srdta@gtdata@map[osnp]
pososnp
reg <- snp.names(r.booqts,beg=pososnp-20000,end=pososnp+20000,chr="1")
reg
a<-readline("PRESS <Enter> TO CONTINUE...")
rr.qts <- snp.subset(r.aqts,snps=reg)
min(rr.qts$P1df)
rr.glm <- scan.glm("qt2~age+sex+CRSNP",data=srdta,snps=reg)
min(rr.glm$P1df)
plot(rr.glm)
points(rr.qts$map,-log10(rr.qts$P1df),col="red")
a<-readline("PRESS <Enter> TO CONTINUE...")
r.2.hap <- scan.haplo("qt2",data=srdta,snps=reg,x.adj=x.adj,n.slide=2)
min(r.2.hap$P1df)
r.3.hap <- scan.haplo("qt2",data=srdta,snps=reg,x.adj=x.adj,n.slide=3)
min(r.3.hap$P1df)
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(rr.glm)
add.plot(rr.qts,col="cyan")
add.plot(r.2.hap,col="red",type="l")
add.plot(r.3.hap,col="green",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")


#
# Using gwaa data in other packages -- haplo.stats
#
osnp
index <- which(srdta@gtdata@snpnames==osnp)
index
r5 <- srdta@gtdata@snpnames[(index-2):(index+2)]
r5
haplo.score(qt2,as.hsgeno(srdta@gtdata[,r5]),x.adj=x.adj)
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# Using gwaa data in other packages -- genetics
#
if (!require(genetics)) stop("you need to install library 'genetics' to see this part of demo")
gdta <- as.genotype(srdta@gtdata[,reg])
ld <- LD(gdta)
plot(ld)
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# generating combined 2D - LD heatmaps
#

# figure out best SNP P-value"
qts <- qtscore(qt2~sex+age,data=srdta)
min(qts$P1df)

# get index of the best SNP"
b<-which.min(qts$P1df)
b
a<-readline("PRESS <Enter> TO CONTINUE...")

# get names of SNPs in +/- 25 kbp region"
reg1 <- snp.names(srdta,beg=srdta@gtdata@map[b]-25000,end=srdta@gtdata@map[b]+25000)
reg1
a<-readline("PRESS <Enter> TO CONTINUE...")

# 2D haplo analysis"
# 2D haplo analysis will be adjusted for sex and age
# CRSNP signifies current SNP
reg1.h2D <- scan.haplo.2D("qt2~sex+age+CRSNP",data=srdta,snps=reg1)
min(reg1.h2D$P1df,na.rm=T)
a<-readline("PRESS <Enter> TO CONTINUE...")


# LD analysis of the region"
# get regional data in genetics format"
reg1.gen <- as.genotype(srdta@gtdata[,reg1])
# get LD"
reg1.LD <- LD(reg1.gen)
a<-readline("PRESS <Enter> TO CONTINUE...")

# now let us plot this stuff"
# first, only 2D association"
# get map"
reg1.map <- srdta@gtdata@map[reg1]
reg1.map
# get -log10-P for association"
lgP <- -log10(reg1.h2D$P1df)
# get color scheme"
maxP <- ceiling(max(lgP,na.rm=T))
maxP
brk <- c(0:(maxP+1))
brk
a<-readline("PRESS <Enter> TO CONTINUE...")
# draw image"
image(x=reg1.map,y=reg1.map,z=lgP,breaks=brk,col=heat.colors(length(brk)-1))
a<-readline("PRESS <Enter> TO CONTINUE...")
# now let us draw LD"
# get D' and transpose"
Dp <- t(reg1.LD$"D'")
# draw it"
image(x=reg1.map,y=reg1.map,z=Dp,breaks=c(-.1,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.1),col=heat.colors(10))
a<-readline("PRESS <Enter> TO CONTINUE...")
# now let us put these together"
image(x=reg1.map,y=reg1.map,z=lgP,breaks=brk,col=heat.colors(length(brk)-1))
image(x=reg1.map,y=reg1.map,z=Dp,breaks=c(-.1,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.1),col=heat.colors(10),add=T)

