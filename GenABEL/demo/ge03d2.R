"?" <- function(...) {
	a<-readline("press ENTER to continue")
}

###################
# Quality Control #
###################

?"LOAD THE DATA"
require(GenABEL.data)
data(ge03d2)

?"ATTACH PHENOTYPIC DATA"
attach(ge03d2@phdata)

?"DESCRIBE TRAIT DATA"
descriptives.trait(ge03d2)

?"DESCRIBE TRAIT DATA, COMPARE CASES AND CONTROLS"
descriptives.trait(ge03d2,by=dm2)

#"DESCRIBE MARKER DATA"
?"...IN ALL SUBJECTS"
descriptives.marker(ge03d2)
?"...IN CONTROLS"
descriptives.marker(ge03d2,ids=(dm2==0))
?"...IN CASES"
descriptives.marker(ge03d2,ids=(dm2==1))

#"GENERATE QQ PLOT FOR HWE P-VALUES IN CONTROLS"
?"...GET SUMMARY FOR SNPS IN CONTROLS"
s <- summary(ge03d2@gtdata[(dm2==0),])
?"...LOOK UP FIRST 10 SNPs"
s[1:10,]
?"...LOOK UP FIRST 10 SNPs, SORTED BY HWE-P-VALUE"
s[order(s$Pexact),][1:10,]
?"...GENERATE QQ PLOT"
estlambda(s[,"Pexact"])

#"GENERATE QQ PLOT FOR HWE P-VALUES IN CASES"
?"...GET SUMMARY FOR SNPS IN CASES"
s <- summary(ge03d2@gtdata[(dm2==1),])
?"...LOOK UP FIRST 10 SNPs, SORTED BY HWE-P-VALUE"
s[order(s$Pexact),][1:10,]
?"...GENERATE QQ PLOT"
estlambda(s[,"Pexact"])

#"QUALITY CONTROL, FIRST PASS"
?"...RUN CHECK.MARKER WITHOUT HWE CHECKS, KEEPING ALL POLYMORPHIC MARKERS"
qc1 <- check.marker(ge03d2,p.level=0,maf=0.0001)
?"...DETAILED SUMMARY OF ERRORS"
summary(qc1)
?"...GENERATE DATA SET 1, WHICH IS RELATIVELY CLEAN"
data1 <- ge03d2[qc1$idok,qc1$snpok]
?"...FIX SPORADIC X-ERRORS"
data1 <- Xfix(data1)
?"...BEFORE WORKING WITH NEW DATA, DETACH PREVIOUS PHENOTYPIC DATA"
detach(ge03d2@phdata)

#"CHECK DATA1 FOR STRONG GENETIC OUTLIERS"
#"...FIRST COMPUTE GENOMIC KINSHIP USING ALL AUTOSOMAL MARKERS"
?"...FROM SAMPLE OF 1000 RANDOM SNPS TO AVOID LD PROBLEMS (THIS PARTICULAR SET)"
data1.gkin <- ibs(data1[,sample(autosomal(data1),1000)],weight="freq")
?"...COMPUTE EUCLIDIAN DISTANCE AS .5-Kinship"
data1.dist <- as.dist(0.5-data1.gkin)
?"...PERFORM PRINCIPAL COMPONENT ANALYSIS OF Kinship (2 COMP)"
data1.mds <- cmdscale(data1.dist)
?"...PLOT THE RESULTS; YOU WILL SEE TWO DISTINCT CLUSTERS"
plot(data1.mds)

#"ANALYSE AND REMOVE OUTLIERS FROM THE DATA"
?"...IDENTIFY SMALLER CLUSTER"
km <- kmeans(data1.mds,centers=2,nstart=1000)
cl1 <- names(which(km$cluster==1))
cl2 <- names(which(km$cluster==2))
if (length(cl1) > length(cl2)) cl1 <- cl2;
cl1
?"...PAINT THE OUTLIERS IN RED"
points(data1.mds[cl1,],pch=19,col="red")
?"...WHO WERE THESE OUTLIERS? BLUE=CASES, GREEN=CONTROLS"
points(data1.mds[data1@phdata$dm2==1,],cex=1.5,col="blue")
points(data1.mds[data1@phdata$dm2==0,],cex=1.5,col="green")
?"...EXCLUDE THE OUTLIERS"
data2 <- data1[!(data1@gtdata@idnames %in% cl1),]
?"...CHECK THE DATA AGAIN, NOW CHECK HWE"
qc2 <- check.marker(data2,hweids=(data2@phdata$dm2==0),maf=0.0001)
?"...SUMMARY OF QC"
summary(qc2)
?"...GENERATE QC-ED DATA"
data3 <- data2[qc2$idok,qc2$snpok]
?"...THESE DATA SHOULD BE COMPLETELY CLEAN"
qc3 <- check.marker(data3,hweids=(data2@phdata$dm2==0),maf=0.0001)

?"DESCRIBE THE FINAL ANALYSIS PHENOTYPES"
descriptives.trait(data2,by=data2@phdata$dm2)

?"DESCRIBE THE FINAL ANALYSIS GENOTYPES"
descriptives.marker(data2)

################
# GWA analysis #
################
#"TRY ANALYSIS OF RAW DATA"
?"...RUN GWA WITH QTSCORE"
an0 <- qtscore(dm2,ge03d2)
?"...PLOT THE ANALYSIS RESULTS"
plot(an0)
?"...CHECK IF POPULATION IS HOMOGENIOUS (Lambda=1)"
lambda(an0)
?"...COMPARE TO CORRECTED P-value"
add.plot(an0,df="Pc1df",col="red")
?"...DESCRIPTIVES OF THE RESULTS, SORTED BY Pc1df"
descriptives.scan(an0,sort="Pc1df")

#"TRY ANALYSIS OF QC'ed DATA"
?"...RUN GWA WITH QTSCORE"
an3 <- qtscore(dm2,data3)
?"...CHECK IF POPULATION IS HOMOGENIOUS (Lambda=1)"
lambda(an3)
?"...PLOT THE ANALYSIS RESULTS FOR CORRECTED P-value"
plot(an3,df="Pc1df")

#"TRY ADJUSTMENT FOR SEX AND AGE"
?"...RUN GWA WITH QTSCORE"
an3.a <- qtscore(dm2~sex+age,data3)
?"...CHECK LAMBDA"
lambda(an3.a)
?"...PLOT THE ANALYSIS RESULTS FOR CORRECTED P-value"
plot(an3.a,df="Pc1df")

#"COMPARE RESULTS"
?"...QC'ED DATA, NO ADJUSTMENT = BLUE"
plot(an3,df="Pc1df",col="blue")
?"...QC'ED DATA, ADJUSTMENT = RED"
add.plot(an3.a,df="Pc1df",col="red")
?"...DATA BEFORE QC = GREEN"
add.plot(an0,df="Pc1df",col="green")

#"STRATIFIED ANALYSIS WITH OBESE CASES"
?"...FIRST ATTACH THE DATA"
attach(data3@phdata)
?"...RUN GWA WITH QTSCORE, EXCLUDING CASES WITH BMI<30"
an3.ob <- qtscore(dm2~sex+age,data3,ids=((bmi>=30 & dm2==1) | dm2==0))
?"...CHECK LAMBDA"
lambda(an3.ob)
?"...PLOT THE ANALYSIS RESULTS FOR CORRECTED P-value"
plot(an3.ob,df="Pc1df")
?"...RUN GWA WITH QTSCORE, EXCLUDING CASES WITH BMI>=30"
an3.nob <- qtscore(dm2~sex+age,data3,ids=((bmi<30 & dm2==1) | dm2==0))
?"...CHECK LAMBDA"
lambda(an3.nob)
?"...PLOT THE ANALYSIS RESULTS FOR CORRECTED P-value"
plot(an3.nob,df="Pc1df")

#"COMPARE RESULTS"
?"...OBESE CASES = RED"
plot(an3.ob,df="Pc1df",col="red")
?"...NON-OBESE CASES = BLUE"
add.plot(an3.nob,df="Pc1df",col="blue")
?"...ALL CASES = GREEN"
add.plot(an3.a,df="Pc1df",col="green")

#"COMPUTE GW SIGNIFICANCE FOR THREE BEST MODELS"
?"...QC'ED DATA, ADJUSTMENT"
an3.a.e <- qtscore(dm2~sex+age,data3,times=100)
?"...DESCRIPTIVES OF THE RESULTS"
descriptives.scan(an3.a.e,sort="Pc1df")
?"...EXCLUDING CASES WITH BMI<30"
an3.ob.e <- qtscore(dm2~sex+age,data3,ids=((bmi>=30 & dm2==1) | dm2==0),times=100)
?"...DESCRIPTIVES OF THE RESULTS"
descriptives.scan(an3.ob.e,sort="Pc1df")
?"...EXCLUDING CASES WITH BMI>=30"
an3.nob.e <- qtscore(dm2~sex+age,data3,ids=((bmi<30 & dm2==1) | dm2==0),times=100)
?"...DESCRIPTIVES OF THE RESULTS"
descriptives.scan(an3.nob.e,sort="Pc1df")

#"COMPARE RESULTS"
?"...OBESE CASES = RED"
plot(an3.ob.e,df="Pc1df",col="red")
?"...NON-OBESE CASES = BLUE"
add.plot(an3.nob.e,df="Pc1df",col="blue")
?"...ALL CASES = GREEN"
add.plot(an3.a.e,df="Pc1df",col="green")

?"DETACH THE DATA"
detach(data3@phdata)

###########################################
# Analysis with stratification correction #
###########################################

#"STRUCTURED ASSOCIATION (SA) ANALYSIS"
?"...CONSTRUCT STRATA"
strata <- 1*(data1@phdata$id %in% cl1)
?"...RUN ANALYSIS"
an4.sa <- qtscore(dm2~sex+age,data1,strata=strata)
?"...CHECK LAMBDA"
lambda(an4.sa)
?"...COMPARE RESULTS; SA = RED"
plot(an4.sa,ylim=c(0,7),col="red",df="Pc1df")
?"...ANALYSIS EXCLUSING OUTLIERS = BLUE"
add.plot(an3.a,col="blue",df="Pc1df")
?"...ANALYSIS WITH OUTLIERS = GREEN"
add.plot(an0,col="green",df="Pc1df")

#"EIGENVALUE ANALYSIS"
an4.eg <- egscore(dm2~sex+age,data1,kin=data1.gkin)
?"...CHECK LAMBDA"
lambda(an4.eg)
?"...COMPARE RESULTS; SA = RED, W/O OUTLIERS = BLUE, WITH OUTLIERS = GREEN"
plot(an4.sa,ylim=c(0,7),col="red",df="Pc1df")
add.plot(an3.a,col="blue",df="Pc1df")
add.plot(an0,col="green",df="Pc1df")
?"...EIGENVECTOR CORRECTION ANALYSIS = BLACK"
add.plot(an4.eg,col="black",cex=1.5,df="Pc1df")

#COMPUTE EMPIRICAL GW P-VALUES
?"...FOR STRUCTURED ASSOCIATION"
an4.sa.e <- qtscore(dm2~sex+age,data1,strata=strata,times=100)
?"...DESCRIPTIVES OF THE RESULTS"
descriptives.scan(an4.sa.e,sort="Pc1df")
?"...FOR EIGENVALUE ANALYSIS"
an4.eg.e <- egscore(dm2~sex+age,data1,kin=data1.gkin,times=100)
?"...DESCRIPTIVES OF THE RESULTS"
descriptives.scan(an4.eg.e,sort="Pc1df")

###############
# Replication #
###############

#"TRY REPLICATION IN OTHER DATA SET"
?"SELECT TOP 10 SNPs FROM ADJUSTED ANALYSIS"
top10 <- rownames(descriptives.scan(an4.sa,sort="Pc1df"))
top10
?"TRY TO REPLICATE IN THE SMALL DATA SET"
require(GenABEL.data)
data(ge03d2c)
confirm <- qtscore(dm2~sex+age,ge03d2c[,top10])
descriptives.scan(confirm)
?"IS EMPIRICAL EXPERIMENT-WISE P-VALUE ALSO OK?"
confirm.e <- emp.qtscore(dm2~sex+age,ge03d2c[,top10],times=10000,bcast=100)
descriptives.scan(confirm.e)
?"CONFIRMED SNPs"
csnps <- rownames(descriptives.scan(confirm))[c(1,3,4)]
csnps
?"CHECK WHAT ARE THE CONFIRMED SNPs ON NCBI"
show.ncbi(csnps)

#####################
# Regional analysis #
#####################

?"DO REGIONAL ANALYSIS"
?"SELECT SNPS IN 300 KB REGION AROUND SNP3"
snp <- csnps[1]
snppos <- data1@gtdata@map[which(data1@gtdata@snpnames==snp)]
reg <- snp.names(data1,begin=snppos-150000,end=snppos+150000)
reg
?"ONE SNP ANALYSES"
reg.qt <- qtscore(dm2~sex+age,data1[,reg],clam=lambda(an4.sa)$estimate)
print(min(reg.qt[,"P1df"]))
reg.glm <- mlreg(dm2~sex+age,data=data2[,reg],trait="binomial") #,clam=an4.sa$lam$est)
print(min(reg.glm[,"P1df"]))
?"HAPLOTYPE ANALYSIS IN SLIDING WINDOW"
reg.h2 <- scan.haplo(dm2~sex+age+CRSNP,data2[,reg],trait.type="binomial")
print(min(reg.h2[,"P1df"]))
reg.h3 <- scan.haplo(dm2~sex+age+CRSNP,data2[,reg],n.slide=3,trait.type="binomial")
print(min(reg.h3[,"P1df"]))
?"PLOT RESULTS"
minp <- min(reg.qt[,"P1df"],reg.glm[,"P1df"],reg.h2[,"P1df"],reg.h3[,"P1df"])
print(minp)
plot(reg.qt,ylim=c(0,ceiling(-log10(minp))))
add.plot(reg.glm,cex=2)
add.plot(reg.h2,col="green",typ="l")
add.plot(reg.h3,col="blue",typ="l")
?"DO 2D HAPLOTYPE ANALYSIS IN 50 KB REGION OF THE PEAK"
bpos <- annotation(reg.h2)$Position[which.min(reg.h2[,"P1df"])]
sreg <- snp.names(data2,begin=bpos-25000,end=bpos+25000)
sreg
sreg.h2D <- scan.haplo.2D(dm2~age+sex+CRSNP,data2[,sreg],trait.type="binomial")
print(min(sreg.h2D$P1df,na.rm=T))
plot(sreg.h2D)
?"DO LD ANALYSIS AND UPDATE THE PLOT"
if (!require(genetics)) stop("you need to install library 'genetics' to see this part of demo")
sreg.LD <- LD(as.genotype(gtdata(data2[,sreg])))
image(sreg.h2D$map,sreg.h2D$map,t(sreg.LD$"D'"),add=T)
