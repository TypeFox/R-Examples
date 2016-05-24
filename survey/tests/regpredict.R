library(survey)
data(api)
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)


## regression estimator of total, three ways
pop<-data.frame(enroll=sum(apipop$enroll, na.rm=TRUE))
npop <- sum(!is.na(apipop$enroll))

api.reg <- svyglm(api.stu~enroll, design=dstrat)
a <- predict(api.reg, newdata=pop, total=npop)
b <- svytotal(~api.stu, calibrate(dstrat, ~enroll, pop=c(npop, pop$enroll)))

all.equal(as.vector(coef(a)),as.vector(coef(b)))
all.equal(as.vector(SE(a)), as.vector(SE(b)))
if(!is.null(getOption("DEBUG"))){ ## uses 6194x6194 matrix
    d <- predict(api.reg, newdata=na.omit(apipop[,"enroll",drop=FALSE]))
    all.equal(as.vector(coef(a)), as.vector(sum(coef(d))))
    all.equal(as.vector(SE(a)), as.vector(sqrt(sum(vcov(d)))))
}

## classical ratio estimator, four ways.
api.reg2 <- svyglm(api.stu~enroll-1, design=dstrat,
                   family=quasi(link="identity", var="mu"))

a <- predict(api.reg2, newdata=pop, total=npop)
b <- svytotal(~api.stu,
              calibrate(dstrat, ~enroll-1, pop= pop$enroll, variance=2))
e <- predict(svyratio(~api.stu, ~enroll, dstrat),total=pop$enroll)

all.equal(as.vector(coef(a)),as.vector(coef(b)))
all.equal(as.vector(SE(a)), as.vector(SE(b)))
all.equal(as.vector(coef(a)),as.vector(e$total))
all.equal(as.vector(SE(a)), as.vector(e$se))
if(!is.null(getOption("DEBUG"))){## uses 6194x6194 matrix
    d <- predict(api.reg2, newdata=na.omit(apipop[,"enroll",drop=FALSE]))
    all.equal(as.vector(coef(a)), as.vector(sum(coef(d))))
    all.equal(as.vector(SE(a)), as.vector(sqrt(sum(vcov(d)))))
}
