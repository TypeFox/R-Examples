data(esoph)
esophlong<-reshape(esoph,direction="long", varying=list(N=c("ncases","ncontrols")),v.names="N",timevar="case",times=1:0)

indiv<-esophlong[rep(1:nrow(esophlong),esophlong$N),]
names(indiv)[4]<-"cancer"

model1<-bigglm(cancer~agegp+tobgp*alcgp, data=indiv, chunksize=100)

library(RSQLite)
sqlite<-dbDriver("SQLite")
conn<-dbConnect(sqlite)
dbWriteTable(conn, "esophtable", indiv)

model2<-bigglm(cancer~agegp+tobgp*alcgp, data=conn, tablename="esophtable", chunksize=100)
