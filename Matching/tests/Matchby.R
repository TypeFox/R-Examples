suppressMessages(library(Matching))

data(lalonde)

X  <- cbind(lalonde$black, lalonde$age, lalonde$educ)
Y  <- lalonde$re78
Tr  <- lalonde$treat

rr2  <- Matchby(Y=Y, Tr=Tr, X=X, M=1, exact=TRUE, by=X[,1],
                ties=TRUE, replace=TRUE, AI=TRUE)
summary(rr2)

rr  <- Match(Y=Y, Tr=Tr, X=X, M=1, exact=TRUE)
summary(rr, full=TRUE)

rr$est-rr2$est
rr$se-rr2$se
rr$se.standard-rr2$se.standard
