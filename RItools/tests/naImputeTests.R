require("RItools")
data(nuclearplants)

testdata<-nuclearplants

testdata$cum.n.fac <- factor(ifelse(c(T,T,T,F),testdata$cum.n,NA))
testdata$pt.log <- as.logical(testdata$pt)
testdata$pt.log[c(3,4)] <- NA
testdata$pt.na <- testdata$pt
testdata$pt.na[3:4] <- NA
testdata$date.mod <- testdata$date
testdata$date.mod[1:2] <- NA
##Previous versions imputed with the mean, so for testing purposes, for now, we use the mean here.
RItools:::naImpute(cost~cap+date+cum.n+as.logical(pt),testdata,impfn=mean.default)

RItools:::naImpute(cost~cap+date.mod+cum.n+as.logical(pt),testdata,impfn=mean.default)

RItools:::naImpute(cost~cap+date+cum.n.fac+as.logical(pt),testdata,impfn=mean.default)

RItools:::naImpute(cost~cap+date+cum.n+pt.log,testdata,impfn=mean.default)

##Assess use of the median explicitly
dat1a<-RItools:::naImpute(pr~date.mod+cum.n.fac+pt.log+pt.na,testdata)

all(dat1a$date.mod[is.na(testdata$date.mod)]==median(testdata$date.mod,na.rm=TRUE))
all(dat1a$pt.log[is.na(testdata$pt.log)]==(median(dat1a$pt.log,na.rm=TRUE)>.5))
all(dat1a$pt.na[is.na(testdata$pt.na)]==(median(dat1a$pt.na,na.rm=TRUE)))


##Check that it works using mean.default() instead of median()
dat1b<-RItools:::naImpute(pr~date.mod+cum.n.fac+pt.log+pt.na,testdata,impfn=mean.default)

all(dat1b$date.mod[is.na(testdata$date.mod)]==mean(testdata$date.mod,na.rm=TRUE))
all(dat1b$pt.log[is.na(testdata$pt.log)]==(mean(dat1b$pt.log,na.rm=TRUE)>.5))
all(dat1b$pt.na[is.na(testdata$pt.na)]==(mean(dat1b$pt.na,na.rm=TRUE)))

##Check that the factor columns are appropriate.
test.mm<-model.matrix(pr~date.mod+cum.n.fac+pt.log+pt.na-1,model.frame(pr~date.mod+cum.n.fac+pt.log+pt.na,testdata,na.action=NULL))
dat2<-RItools:::naImpute(pr~date.mod+cum.n.fac1+cum.n.fac14+pt.logTRUE+pt.na,data.frame(test.mm),impfn=mean.default)
all(dat2$cum.n.fac1[is.na(testdata$cum.n.fac)]==(mean(dat2$cum.n.fac1,na.rm=TRUE)))
all(dat2$cum.n.fac14[is.na(testdata$cum.n.fac)]==(mean(dat2$cum.n.fac14,na.rm=TRUE)))
