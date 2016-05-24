# relInterIntra
# gives the reliability coefficients in the article by 
 # Eliasziw et. al. 1994; Phys. Therapy 74.8; 777-788.
# all references in this function are to this article.
# Arguments
# x: data frame representing a data structure as in Table 1 (p 779), 
#  with consecutive measurements for each rater in adjacent columns.
#  i.e. rater1measure1 rater1measure2 ... rater2measure1 rater2measure2 
# rho0inter: null hypothesis value of the interrater reliability coefficient
# rho0intra: null hypothesis value of the intrarater reliability coefficient
# conf.level: confidence level of the one-sided confidence intervals reported
# for the reliability coefficients
# output reformatted as an "irrlist" stucture - Jim Lemon 2009-05-27

relInterIntra<-function(x,nrater=1,raterLabels=NULL,
 rho0inter=0.6,rho0intra=0.8,conf.level=.95) {

 xdim<-dim(x)
 nsubj<-xdim[1]
 nmeas<-xdim[2]/nrater
 if(is.null(raterLabels)) raterLabels<-letters[1:nrater]
 Frame1<-data.frame(cbind(rep(1:nsubj,nrater*nmeas),
  rep(1:nrater,rep(nsubj*nmeas,nrater)),
  rep(rep(1:nmeas,rep(nsubj,nmeas)),nrater),
  matrix(as.matrix(x),ncol=1)))
 names(Frame1)<-c('Subject','Rater','Repetition','Result')
 Frame1$Subject<-factor(Frame1$Subject)
 Frame1$Rater<-factor(Frame1$Rater,labels=raterLabels)
 Frame1$Repetition<-factor(Frame1$Repetition)
 # this and following two commands:
 # aliases for compatibility with Eliasziw et. al. notation
 nn<-nsubj
 tt<-nrater
 mm<-nmeas
 aovFull<-aov(Result~Subject*Rater,data=Frame1)
 meanSquares<-summary(aovFull)[[1]][,3]
 for(raterAct in 1:tt) {
  raterActCat<-raterLabels[raterAct]
  aovAct<-aov(Result~Subject,data=Frame1[Frame1$Rater==raterActCat,])
  meanSquares<-c(meanSquares,summary(aovAct)[[1]][2,3])
 }
 names(meanSquares)<-
  c('MSS','MSR','MSSR','MSE',paste('MSE',levels(Frame1$Rater),sep=''))
 MSS<-meanSquares[1]
 MSR<-meanSquares[2]
 MSSR<-meanSquares[3]
 MSE<-meanSquares[4]
 # the same for random and fixed, see table 2 (p. 780) and 3 (p. 281)
 MSEpart<-meanSquares[-(1:4)] 
 sighat2Srandom<-(MSS-MSSR)/(mm*tt)
 sighat2Rrandom<-(MSR-MSSR)/(mm*nn)
 sighat2SRrandom<-(MSSR-MSE)/mm
 # the same for random and fixed, see table 2 (p. 780) and 3 (p. 281)
 sighat2e<-MSE 
 sighat2Sfixed<-(MSS-MSE)/(mm*tt)
 sighat2Rfixed<-(MSR-MSSR)/(mm*nn)
 sighat2SRfixed<-(MSSR-MSE)/mm
 # the same for random and fixed, see table 2 (p. 780) and 3 (p. 281)
 sighat2e.part<-MSEpart 
 rhohat.inter.random<-sighat2Srandom/
  (sighat2Srandom+sighat2Rrandom+sighat2SRrandom+sighat2e)
 rhohat.inter.fixed<-(sighat2Sfixed-sighat2SRfixed/tt)/
  (sighat2Sfixed+(tt-1)*sighat2SRfixed/tt+sighat2e)
 rhohat.intra.random<-(sighat2Srandom+sighat2Rrandom+sighat2SRrandom)/
  (sighat2Srandom+sighat2Rrandom+sighat2SRrandom+sighat2e)
 rhohat.intra.fixed<-(sighat2Sfixed+(tt-1)*sighat2SRfixed/tt)/
  (sighat2Sfixed+(tt-1)*sighat2SRfixed/tt+sighat2e)
 rhohat.intra.random.part<-(sighat2Srandom+sighat2Rrandom+sighat2SRrandom)/
  (sighat2Srandom+sighat2Rrandom+sighat2SRrandom+sighat2e.part)
 rhohat.intra.fixed.part<-(sighat2Sfixed+(tt-1)*sighat2SRfixed/tt)/
  (sighat2Sfixed+(tt-1)*sighat2SRfixed/tt+sighat2e.part)
 Finter<-(1-rho0inter)*MSS/((1+(tt-1)*rho0inter)*MSSR)
 Finter.p<-1-pf(Finter,df1=nn-1,df2=(nn-1)*(tt-1))
 alpha<-1-conf.level
 nu1<-(nn-1)*(tt-1)*
  (tt*rhohat.inter.random*(MSR-MSSR)+
  nn*(1+(tt-1)*rhohat.inter.random)*MSSR+
  nn*tt*(mm-1)*rhohat.inter.random*MSE)^2/
  ((nn-1)*(tt*rhohat.inter.random)^2*MSR^2+
  (nn*(1+(tt-1)*rhohat.inter.random)-tt*rhohat.inter.random)^2*MSSR^2+
  (nn-1)*(tt-1)*(nn*tt*(mm-1))*rhohat.inter.random^2*MSE^2)
 nu2<-(nn-1)*(tt-1)*
  (nn*(1+(tt-1)*rhohat.inter.fixed)*MSSR+
  nn*tt*(mm-1)*rhohat.inter.fixed*MSE)^2/
  ((nn*(1+(tt-1)*rhohat.inter.fixed))^2*MSSR^2+
  (nn-1)*(tt-1)*(nn*tt*(mm-1))*rhohat.inter.fixed^2*MSE^2)
 F1<-qf(1-alpha,df1=nn-1,df2=nu1)
 F2<-qf(1-alpha,df1=nn-1,df2=nu2)
 lowinter.random<-nn*(MSS-F1*MSSR)/
  (nn*MSS+F1*(tt*(MSR-MSSR)+nn*(tt-1)*MSSR+nn*tt*(mm-1)*MSE))
 lowinter.random<-min(c(lowinter.random,1))
 lowinter.fixed<-nn*(MSS-F2*MSSR)/
  (nn*MSS+F2*(nn*(tt-1)*MSSR+nn*tt*(mm-1)*MSE))
 lowinter.fixed<-min(c(lowinter.fixed,1))
 Fintra<-(1-rho0intra)*MSS/((1+(mm-1)*rho0intra)*MSE*tt)
 Fintra.p<-1-pf(Fintra,df1=nn-1,df2=nn*(mm-1))
 Fintra.part<-(1-rho0intra)*MSS/((1+(mm-1)*rho0intra)*MSEpart*tt)
 Fintra.part.p<-1-pf(Fintra.part,df1=nn-1,df2=nn*(mm-1))
 F3<-qf(1-alpha,df1=nn-1,df2=nn*(mm-1))
 lowintra<-(MSS/tt-F3*MSE)/(MSS/tt+F3*(mm-1)*MSE)
 lowintra<-min(c(lowintra,1))
 F4<-qf(1-alpha,df1=nn-1,df2=nn*(mm-1))
 lowintra.part<-(MSS/tt-F4*MSEpart)/(MSS/tt+F4*(mm-1)*MSEpart)
 for(raterAct in 1:tt)
  lowintra.part[raterAct]<-min(lowintra.part[raterAct],1)
 SEMintra<-sqrt(MSE)
 SEMintra.part<-sqrt(MSEpart)
 SEMinter.random<-sqrt(sighat2Rrandom+sighat2SRrandom+sighat2e)
 SEMinter.fixed<-sqrt(sighat2SRfixed+sighat2e)
 rels<-list(method="Inter/Intrarater reliability",
  subjects=nsubj,raters=nrater,irr.name="rhohat",
  value=list(rohat=c(rhohat.inter.random,rhohat.intra.random,
  rhohat.inter.fixed,rhohat.intra.fixed,
  rhohat.intra.random.part,rhohat.intra.fixed.part),
  Fs=c(Finter,Fintra,Fintra.part),
  pvalue=c(Finter.p,Fintra.p,Fintra.part.p),
  lowvalue=c(lowinter.random,lowinter.fixed,lowintra,lowintra.part),
  sem=c(SEMintra,SEMintra.part,SEMinter.random,SEMinter.fixed)),
  stat.name="nil",statistic=NULL)
  class(rels)<-"irrlist"
 names(rels$value$rohat)<-
  c('rhohat.inter.random','rhohat.intra.random',
  'rhohat.inter.fixed','rhohat.intra.fixed',
  paste('rhohat.intra.random.part',raterLabels,sep='.'),
  paste('rhohat.intra.fixed.part',raterLabels,sep='.'))
 names(rels$value$Fs)<-
  c('Finter','Fintra',paste('Fintra',raterLabels,sep='.'))
 names(rels$value$pvalue)<-
  c("pvalue.Finter","pvalue.Fintra",paste('pvalue.Fintra',raterLabels,sep='.'))
 names(rels$value$lowvalue)<-c('lowinter.random','lowinter.fixed','lowintra',
  paste('lowintra',raterLabels,sep='.'))
 names(rels$value$sem)<-
  c('SEMintra',paste('SEMintra.part',raterLabels,sep='.'),
  'SEMinter.random','SEMinter.fixed')
 return(rels)
}
