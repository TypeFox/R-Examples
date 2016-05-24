.seerHome="~/data/SEER" # Assume here that the SEER data is in /data/SEER. Change this if not. 
rm(list=ls()) # note: dot variables defined above persist through such cleanings
library(RSQLite) # install.packages("RSQLite")
m=dbDriver("SQLite")
con=dbConnect(m,dbname=file.path(.seerHome,"00/all.db"))
da=dbGetQuery(con,"SELECT * from respir where ICD9>=1620 and ICD9<=1629 
and histo3=8140 and seqnum<2") # adenos
# da=transform(da,dwd=(COD>0),hist= "adenocarcinoma" ) # not alive; 0=alive
da=transform(da,dwd=((COD>=20010)&(COD<=37000)),hist= "adenocarcinoma" ) # dead with cancer, 0=alive
MFcnta=summary(as.factor(da$sex))
head(da)
ds=dbGetQuery(con,"SELECT * from respir where ICD9>=1620 and ICD9<=1629 
and histo3>=8070 and histo3<=8079 and seqnum<2") # squames
MFcnts=summary(as.factor(ds$sex))
# ds=transform(ds,dwd=(COD>0),hist= "squamous cell") # not alive; 0=alive
ds=transform(ds,dwd=((COD>=20010)&(COD<=37000)),hist= "squamous cell" ) # dead with cancer, 0=alive
d=rbind(da,ds)

sapply(d,class)
head(d)
library(survival)
# graphics.off()
if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
if(length(grep("darwin",R.Version()$os))) windows <- function( ... ) quartz( ... )
windows(width=9,height=5)
par(mfrow=c(1,2),mar=c(4.5,4.1,1,1))
rb=c("blue","red") #males blue, femals red 
plot(S<-survfit(Surv(surv,dwd)~hist,data = subset(d,sex==1)), # main="SEER 2000-2010",
     xlab="Months",ylab="Survival", xlim=c(0,60) , #ylim=c(0.6,1),
     col=rb, main="Adeno Vs Squame")
print(S) 
legend("topleft",c("Males"),bty="n")

plot(S<-survfit(Surv(surv,dwd)~hist,data = subset(d,sex==2)), ,
     xlab="Months",ylab="Survival", xlim=c(0,60), #ylim=c(0.6,1),
     col=rb, main="Adeno Vs Squame")
print(S) 
legend("topleft",c("Females"),bty="n")
