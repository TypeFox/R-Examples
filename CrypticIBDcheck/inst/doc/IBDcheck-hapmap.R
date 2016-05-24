### R code from vignette source 'IBDcheck-hapmap.Rnw'

###################################################
### code chunk number 1: IBDcheck-hapmap.Rnw:97-98
###################################################
library(CrypticIBDcheck)


###################################################
### code chunk number 2: IBDcheck-hapmap.Rnw:121-131 (eval = FALSE)
###################################################
## lwkdat<-vector(mode="list",length=22)
## names(lwkdat)<-paste("chr",1:22,sep="")
## for(i in 1:22) {
## uu<-paste("http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/",
## "latest_phaseIII_ncbi_b36/hapmap_format/polymorphic/genotypes_chr",
## i,
## "_LWK_phase3.2_nr.b36_fwd.txt.gz",
## sep="")
## lwkdat[[i]]<-read.HapMap.data(uu)
## }


###################################################
### code chunk number 3: IBDcheck-hapmap.Rnw:146-152 (eval = FALSE)
###################################################
## snp.data<-lwkdat[[1]]$snp.data
## snp.support<-lwkdat[[1]]$snp.support[,c("Chromosome","Position")]
## for(i in 2:22) {
##   snp.data<-cbind(snp.data,lwkdat[[i]]$snp.data)
##   snp.support<-rbind(snp.support,lwkdat[[i]]$snp.support[,c("Chromosome","Position")])
## }


###################################################
### code chunk number 4: IBDcheck-hapmap.Rnw:155-158 (eval = FALSE)
###################################################
## dd<-duplicated(snp.support)
## snp.support<-snp.support[!dd,]
## snp.data<-snp.data[,!dd]


###################################################
### code chunk number 5: IBDcheck-hapmap.Rnw:166-168 (eval = FALSE)
###################################################
## dat<-new.IBD(snp.data,Chromosome=snp.support$Chromosome,
## Position=snp.support$Position, popsam=rep(TRUE,nrow(snp.data)))


###################################################
### code chunk number 6: IBDcheck-hapmap.Rnw:180-181 (eval = FALSE)
###################################################
## system("plink --no-web --help")


###################################################
### code chunk number 7: IBDcheck-hapmap.Rnw:190-191
###################################################
source(file.path(system.file(package="CrypticIBDcheck"),"scripts","thin.R"))


###################################################
### code chunk number 8: IBDcheck-hapmap.Rnw:217-218 (eval = FALSE)
###################################################
## t.dat<-thin(dat,win=100,shift=25,r2thresh=0.005)


###################################################
### code chunk number 9: IBDcheck-hapmap.Rnw:243-247 (eval = FALSE)
###################################################
## ss<-sim.control(simulate=TRUE,fitLD=FALSE,
## rships=c("unrelated", "MZtwins", "parent-offspring", "full-sibs", 
## "half-sibs", "cousins"), nsim=rep(200,6))
## cibd<-IBDcheck(t.dat,simparams=ss)


###################################################
### code chunk number 10: IBDcheck-hapmap.Rnw:260-261 (eval = FALSE)
###################################################
## ibdpairs<-plot(cibd)


###################################################
### code chunk number 11: IBDcheck-hapmap.Rnw:437-444 (eval = FALSE)
###################################################
## uu<-paste("http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/",
## "latest_phaseIII_ncbi_b36/relationships_w_pops_121708.txt",sep="")
## hapmap.info<-read.table(uu,header=TRUE,as.is=TRUE)
## subject.support<-hapmap.info[hapmap.info$population=="CEU",]
## parent<-(subject.support$mom==0 | subject.support$dad==0)
## subject.support<-subject.support[parent,]
## rm(hapmap.info)


###################################################
### code chunk number 12: IBDcheck-hapmap.Rnw:457-459 (eval = FALSE)
###################################################
## id=rownames(snp.data)
## subject.support=subject.support[match(id,subject.support$IID),]


