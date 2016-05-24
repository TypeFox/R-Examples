
library(survey)
library(RSQLite)

data(api)
apiclus1$api_stu<-apiclus1$api.stu
apiclus1$comp_imp<-apiclus1$comp.imp
dclus1<-svydesign(id=~dnum, weights=~pw, fpc=~fpc,data=apiclus1)
dbclus1<-svydesign(id=~dnum, weights=~pw, fpc=~fpc,
data="apiclus1",dbtype="SQLite", dbname=system.file("api.db",package="survey"))

m<-svymean(~api00+stype,dclus1)
m.db<-svymean(~api00+stype, dbclus1)
all.equal(coef(m),coef(m.db))
all.equal(vcov(m), vcov(m.db))

r<-svyratio(~api_stu, ~enroll, design=dclus1)
r.db<-svyratio(~api_stu, ~enroll, design=dbclus1)
all.equal(coef(r), coef(r.db))
all.equal(SE(r), SE(r.db))

b<-svyby(~api99+api00,~stype, design=dclus1, svymean, deff=TRUE)
b.db<-svyby(~api99+api00,~stype, design=dbclus1,svymean, deff=TRUE)
all.equal(coef(b), coef(b.db))
all.equal(SE(b), SE(b.db))
all.equal(deff(b), deff(b.db))

l<-svyglm(api00~api99+mobility, design=dclus1)
l.db<-svyglm(api00~api99+mobility, design=dbclus1)
all.equal(coef(l),coef(l.db))
all.equal(vcov(l), vcov(l.db))

dclus1<-update(dclus1, apidiff=api00-api99)
dclus1<-update(dclus1, apipct= apidiff/api99)
dbclus1<-update(dbclus1, apidiff=api00-api99)
dbclus1<-update(dbclus1, apipct= apidiff/api99)

u<-svymean(~api00+apidiff+apipct, dclus1)
u.db<-svymean(~api00+apidiff+apipct, dbclus1)
all.equal(u, u.db)

all.equal(nrow(dclus1),nrow(dbclus1))
all.equal(nrow(subset(dclus1,stype=="E")),
          nrow(subset(dbclus1,stype=="E")))

## replicate weights
rclus1<-as.svrepdesign(dclus1)
db_rclus1<-svrepdesign(weights=~pw, repweights="wt[1-9]+", type="JK1", scale=(1-15/757)*14/15,
data="apiclus1rep",dbtype="SQLite", dbname=system.file("api.db",package="survey"),combined.weights=FALSE)
m<-svymean(~api00+api99,rclus1)
m.db<-svymean(~api00+api99,db_rclus1)
all.equal(m,m.db)

summary(db_rclus1)

s<-svymean(~api00, subset(rclus1, comp_imp=="Yes"))
s.db<-svymean(~api00, subset(db_rclus1, comp_imp=="Yes"))
all.equal(s,s.db)
