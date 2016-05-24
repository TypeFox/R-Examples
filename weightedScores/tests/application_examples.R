library(weightedScores)
data(childvisit)
names(childvisit)
season1<-childvisit$quarter
season1[season1>1]<-0
xdat<-cbind(1,childvisit$sex,childvisit$age,childvisit$matsmst,season1)
ydat<-childvisit$hospvis
id<-childvisit$id
tvec<-childvisit$quarter

out.Poisson.exch<-wtsc.wrapper(xdat,ydat,id,tvec,margmodel="poisson",corstr="exch",iprint=T)

out.Poisson.ar<-wtsc.wrapper(xdat,ydat,id,tvec,margmodel="poisson",corstr="ar",iprint=T)

out.Poisson.unstr<-wtsc.wrapper(xdat,ydat,id,tvec,margmodel="poisson",corstr="unstr",iprint=T)

out.nb1.exch<-wtsc.wrapper(xdat,ydat,id,tvec,margmodel="nb1",corstr="exch",iprint=T)

out.nb1.ar<-wtsc.wrapper(xdat,ydat,id,tvec,margmodel="nb1",corstr="ar",iprint=T)

out.nb1.unstr<-wtsc.wrapper(xdat,ydat,id,tvec,margmodel="nb1",corstr="unstr",iprint=T)

out.nb2.exch<-wtsc.wrapper(xdat,ydat,id,tvec,margmodel="nb2",corstr="exch",iprint=T)

out.nb2.ar<-wtsc.wrapper(xdat,ydat,id,tvec,margmodel="nb2",corstr="ar",iprint=T)

out.nb2.unstr<-wtsc.wrapper(xdat,ydat,id,tvec,margmodel="nb2",corstr="unstr",iprint=T)



y2<-ydat
y2[ydat>1]<-1

out.Bernoulli.logit.exch<-wtsc.wrapper(xdat,y2,id,tvec,margmodel="bernoulli",corstr="exch",iprint=T)

out.Bernoulli.logit.ar<-wtsc.wrapper(xdat,y2,id,tvec,margmodel="bernoulli",corstr="ar",iprint=T)

out.Bernoulli.logit.unstr<-wtsc.wrapper(xdat,y2,id,tvec,margmodel="bernoulli",corstr="unstr",iprint=T)

out.Bernoulli.probit.exch<-wtsc.wrapper(xdat,y2,id,tvec,margmodel="bernoulli",corstr="exch",link="probit",iprint=T)

out.Bernoulli.probit.ar<-wtsc.wrapper(xdat,y2,id,tvec,margmodel="bernoulli",corstr="ar",link="probit",iprint=T)

out.Bernoulli.probit.unstr<-wtsc.wrapper(xdat,y2,id,tvec,margmodel="bernoulli",corstr="unstr",link="probit",iprint=T)





