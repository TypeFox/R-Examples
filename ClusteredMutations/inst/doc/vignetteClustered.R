### R code from vignette source 'vignetteClustered.Rnw'

###################################################
### code chunk number 1: vignetteClustered.Rnw:45-47
###################################################
library(ClusteredMutations)
data(PD4107a)


###################################################
### code chunk number 2: vignetteClustered.Rnw:50-52
###################################################
data.showers<-showers(data=PD4107a, chr=Chr, position=Position, min=2, max=10)
head(data.showers, n=10)


###################################################
### code chunk number 3: vignetteClustered.Rnw:57-82
###################################################
extra <- factor(c(),levels=c("T>C","T>G","T>A","C>T","C>G","C>A"))
extra[PD4107a$Ref_base=="A" & PD4107a$Mutant_base=="G"]<-"T>C"
extra[PD4107a$Ref_base=="T" & PD4107a$Mutant_base=="C"]<-"T>C"
extra[PD4107a$Ref_base=="A" & PD4107a$Mutant_base=="C"]<-"T>G"
extra[PD4107a$Ref_base=="T" & PD4107a$Mutant_base=="G"]<-"T>G"
extra[PD4107a$Ref_base=="A" & PD4107a$Mutant_base=="T"]<-"T>A"
extra[PD4107a$Ref_base=="T" & PD4107a$Mutant_base=="A"]<-"T>A"
extra[PD4107a$Ref_base=="G" & PD4107a$Mutant_base=="A"]<-"C>T"
extra[PD4107a$Ref_base=="C" & PD4107a$Mutant_base=="T"]<-"C>T"
extra[PD4107a$Ref_base=="G" & PD4107a$Mutant_base=="C"]<-"C>G"
extra[PD4107a$Ref_base=="C" & PD4107a$Mutant_base=="G"]<-"C>G"
extra[PD4107a$Ref_base=="G" & PD4107a$Mutant_base=="T"]<-"C>A"
extra[PD4107a$Ref_base=="C" & PD4107a$Mutant_base=="A"]<-"C>A"
PD4107a$extra<-extra

rainfall<-imd(data=PD4107a,chr=Chr,position=Position,extra=extra)

plot(rainfall$number, rainfall$log10distance, col=c("yellow", "green", 
      "pink", "red", "black", "blue")[rainfall$extra], pch=20, 
      ylab="Intermutation distance (bp)", xlab="PD4107a", yaxt="n")
axis(2, at=c(0, 1, 2, 3, 4, 6), labels=c("1", "10", "100", "1000",
    "10000", "1000000"), las=2, cex.axis=0.6)
legend("topleft", legend = levels(rainfall$extra), col=c("yellow", 
      "green", "pink", "red", "black", "blue"), pch=20, horiz=TRUE, 
      text.font=4, bg='lightblue')


###################################################
### code chunk number 4: vignetteClustered.Rnw:93-122
###################################################
set.seed(42)
position<-c( c((runif(1001,min=1,max=10000001))),
c(c(10110001,10110011,10110021,10110031,10110041,10120000,
10120001,10120011,10120021,10120031,10130000,
10130001,10130011,10130021,10130031,10140000,
10140001,10140011,10140021,10140031,10150000,
10150001,10150011,10150021,10150031,10160000,
10160001,10160011,10160021,10160031,10170000,
10170001,10170011,10170021,10170031,10180000,
10180001,10180011,10180021,10180031,10190000,
10210001,10210011,10210021,10210031,10220000,
10220001,10220011,10220021,10220031,10230000,
10230001,10230011,10230021,10230031,10240000,
10240001,10240011,10240021,10240031,10250000,
10250001,10250011,10250021,10250031,10260000,
10260001,10260011,10260021,10260031,10270000,
10270001,10270011,10270021,10270031,10280000,
10280001,10280011,10280021,10280031,10290000) + 
round(runif(81,min=0,max=500))), 
c(round(runif(991,min=10296000,max=20200000))))

rainfall<-imd(position=position)
#Rainfall plot for PD4107a cancer sample;
plot(rainfall$number, rainfall$log10distance, pch=20, 
     ylab="Intermutation distance (bp)", xlab="Example", yaxt="n")
axis(2, at=c(0, 1, 2, 3, 4, 6), labels=c("1", "10", "100", "1000",
    "10000", "1000000"), las=2, cex.axis=0.6)
theta <- seq(0, 2 * pi, length = 200)
lines(x = 100 * cos(theta) + 1050, y = sin(theta) + 1.5, col="red")


###################################################
### code chunk number 5: vignetteClustered.Rnw:130-131
###################################################
showers(position=position)


###################################################
### code chunk number 6: vignetteClustered.Rnw:134-135
###################################################
showers(data=PD4107a,chr=Chr,position=Position)


###################################################
### code chunk number 7: vignetteClustered.Rnw:139-141
###################################################
example1<-c(1,101,201,299,301,306,307,317,318,320,418,518,528,628)
10**(dissmutmatrix(position=example1,upper=TRUE))


###################################################
### code chunk number 8: vignetteClustered.Rnw:155-160
###################################################
mut.matrix <- dissmutmatrix(data=PD4107a, chr=Chr,
position=Position, subset=6)
dissplot(mut.matrix, method=NA, options=list( col = c("black",
"navy", "blue", "cyan", "green", "yellow", "orange", "red",
"darkred", "darkred", "white")))


