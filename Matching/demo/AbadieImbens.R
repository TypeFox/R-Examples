# Replication of Guido Imbens lalonde_exper_04feb2.m file
# See http://elsa.berkeley.edu/~imbens/estimators.shtml
#
# Note that the implications of the 'exact' options differ between the
# two programs

data(lalonde)

X  <- lalonde$age
Z  <- X;             
V  <- lalonde$educ;
Y  <- lalonde$re78/1000;
T  <- lalonde$treat;
w.educ=exp((lalonde$educ-10.1)/2);

res  <- matrix(nrow=1,ncol=3)

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)

res[1,]  <- cbind(1,rr$est,rr$se)


X  <- cbind(lalonde$age, lalonde$educ, lalonde$re74, lalonde$re75)
rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)

res  <- rbind(res,cbind(2,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=3,BiasAdj=FALSE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(4,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATT",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(5,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATC",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(6,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=2,Var.calc=0,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(7,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=3,Var.calc=0,
             Weight.matrix=diag(4), sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(8,rr$est,rr$se))


rr  <- Match(Y=Y,Tr=T,X=X,Z=X,V=V,estimand="ATE",M=1,BiasAdj=TRUE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(9,rr$est,rr$se))

Z  <- cbind(lalonde$married, lalonde$age)
rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=TRUE,Weight=1,Var.calc=0,sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(10,rr$est,rr$se))

V  <- lalonde$age
rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,exact=TRUE,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(11,rr$est,rr$se))

V  <- cbind(lalonde$married, lalonde$u74)
rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,exact=TRUE,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(12,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,sample=FALSE);
summary(rr)
res  <- rbind(res,cbind(13,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=3,sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(14,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,
             weights=w.educ,sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(15,rr$est,rr$se))


V  <- lalonde$age
Z  <- cbind(lalonde$married, lalonde$age)
X  <- cbind(lalonde$age, lalonde$educ, lalonde$re74, lalonde$re75)
weight  <- w.educ
Weight.matrix  <- diag(4)

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,
             sample=FALSE, M=3, estimand="ATT", BiasAdj=TRUE, Weight=3, exact=TRUE,Var.calc=3,
             weights=w.educ, Weight.matrix=Weight.matrix);
summary(rr)
res  <- rbind(res,cbind(75,rr$est,rr$se))


V  <- lalonde$married;
Z  <- cbind(lalonde$age, lalonde$re75);
X  <- cbind(lalonde$age, lalonde$educ, lalonde$re74);

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,
             sample=TRUE, M=3, estimand="ATE", BiasAdj=TRUE, Weight=2, exact=TRUE,Var.calc=0,
             weights=w.educ);
summary(rr)
res  <- rbind(res,cbind(76,rr$est,rr$se))

cat("\nResults:\n")
print(res)
