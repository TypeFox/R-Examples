# Sex Difference in Myelodysplastic Syndrome Survival and Balance in Randomized
# Clinical Trials. Radivoyevitch T, Saunthararajah Y. J Clin Oncol. 2013 Nov 18.
# [Epub ahead of print] No abstract available. PMID:  24248691

# diffs in sex proportions between arms A and B, 2/22=9% vs 18/43=42%
(X <-  matrix(c(2, 20, 18, 25), nrow = 2))
fisher.test(X)



.seerHome="~/data/SEER" # Assume here that the SEER data is in /data/SEER. Change this if not. 
rm(list=ls()) # note: dot variables defined above persist through such cleanings
library(RSQLite) # install.packages("RSQLite")
m=dbDriver("SQLite")
con=dbConnect(m,dbname=file.path(.seerHome,"00/all.db"))
dbListTables(con)
dbListFields(con,"other")
d=dbGetQuery(con, 
             "SELECT * from other where histo3>9979 and histo3<9990")#MDS
head(d)
summary(d)
summary(as.factor(d$COD))
MFcnts=summary(as.factor(d$sex))
summary(as.factor(d$surv))
# nd=transform(d,dwd=((COD>=20010)&(COD<=37000)) ) # dead with disease
nd=transform(d,dwd=(COD>0) ) # not alive; 0=alive
# d[d$COD==0,c("COD","surv")]
summary(as.factor(nd$surv))
# nd=transform(nd,decade=cut(yrdx,br=c(1999,2003,2011),labels=c("2000-2003","2004-2010"))) 
# # nd=transform(nd,agedxF=cut(agedx,br=c(20,69,77,83,1010),labels=c("35-68","69-76","77-83","84-106"))) 
# nd=transform(nd,agedxF=cut(agedx,br=c(20,77,1010),labels=c("<77",">77"))) 
# head(nd[,c("yrdx","decade")],40)

library(survival)
graphics.off()
if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
if(length(grep("darwin",R.Version()$os))) windows <- function( ... ) quartz( ... )
windows(width=9,height=5)
par(mfrow=c(1,2),mar=c(4.5,4.1,1,1),cex=1.4,lwd=1.5)
rb=c("blue","red") #males blue, femals red 
bbr=c("black","blue","red") #males blue, femals red 
plot(survfit(Surv(surv,dwd)~sex,data = nd),  main="MDS", #: SEER 2000-2010",
     xlab="Months",ylab="Survival", xlim=c(0,60) , #ylim=c(0.6,1),
     col=rb)
print(S<-survdiff(Surv(surv,dwd)~sex,data = nd)) 
options(digits=10)
dput(S)
mtext("A",side=3,line=0,cex=1.5,adj=-.3,font=2)
# legend("topright",c("Males","Females"),text.col=rb,bty="n")
legend("bottomleft",expression(paste("P < ",10^-15)),bty="n")
legend("topright",c("Cases",paste("Males =",MFcnts[1]),paste("Females =",MFcnts[2])),text.col=bbr,bty="n")
# legend("topleft",c("MDS"),bty="n")

d=dbGetQuery(con, "SELECT * from lymyleuk where ICD9=2050")#AML
# nd=transform(d,dwd=((COD>=20010)&(COD<=37000)) ) # dead with disease
MFcnts=summary(as.factor(d$sex))
nd=transform(d,dwd=(COD>0) ) # not alive; 0=alive
plot(survfit(Surv(surv,dwd)~sex,data = nd),main="AML", #: SEER 2000-2010",
     xlab="Months",ylab="Survival", xlim=c(0,60), #ylim=c(0.6,1),
     col=rb)
print(S<-survdiff(Surv(surv,dwd)~sex,data = nd)) 
mtext("B",side=3,line=0,cex=1.5,adj=-.3,font=2)
# legend("topright",c("Males","Females"),text.col=rb,bty="n")
legend("bottomleft",c("P=0.051"),bty="n")
# legend("topleft",c("AML"),bty="n")
legend("topright",c("Cases",paste("Males =",MFcnts[1]),paste("Females =",MFcnts[2])),text.col=bbr,bty="n")
# legend("topright",c(paste("Males =",MFcnts[1]),paste("Females =",MFcnts[2])),text.col=rb,bty="n")
