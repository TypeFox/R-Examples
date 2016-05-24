################################################################################
subsamp.within<-function (lev,ni=10) {
#subsample ni observations per levels of lev
    y <- 1:length(lev)
    nlev <- nlevels(factor(lev))
    nl <- as.integer(factor(lev))
    x <- list()
    for (i in 1:nlev) if (length(y[nl == i]) > ni) 
        x[[i]] <- sample(y[nl == i],size=ni)
    else x[[i]] <- y[nl == i]
    return(unlist(x))
}
##################################################################################
#' @title Calculates corrected Assignment Index
#' @description Calculates corrected Assignment Index as described in \href{http://onlinelibrary.wiley.com/doi/10.1046/j.1365-294X.2002.01496.x/abstract}{Goudet etal. (2002)}
#' @usage AIc(dat)
#' @param dat a data frane with nlocs+1 columns, 
#' @return aic  The corrected assignment index of each individual
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @references \href{http://onlinelibrary.wiley.com/doi/10.1046/j.1365-294X.2002.01496.x/abstract}{Goudet J, Perrin N, Waser P (2002)} Tests for sex-biased dispersal 
#'  using bi-parentally inherited genetic markers 11, 1103:1114
#' @export
##################################################################################
AIc<-function(dat){
dat[,1]<-as.integer(as.factor(dat[,1]))
popfreq<-pop.freq(dat) #freq per pop
indfreq<-pop.freq(cbind(1:dim(dat)[1],dat[,-1])) # freq per ind
id.al<-lapply(indfreq,function(y) apply(y,2,fun2<-function(z) which(z>0))) #which allele in which ind
nloc<-dim(dat)[2]-1 
npop<-length(table(dat[,1]))
nind<-dim(dat)[1]
ail<-matrix(numeric(nind*nloc),ncol=nloc)
for (il in 1:nloc){
  for (ii in 1:nind){
    dum<-as.vector(id.al[[il]][[ii]])
    if (length(dum)!=0){
      if (length(dum)==1)
        ail[ii,il]<-2*log(popfreq[[il]][[dum,dat[ii,1]]])
      else
        ail[ii,il]<-log(2*popfreq[[il]][[dum[1],dat[ii,1]]]*popfreq[[il]][[dum[2],dat[ii,1]]])
    }
    else ail[ii,il]<-NA
  }
}
#mean assignment per locus per pop
mean.ail<-apply(ail,2,function(x) {dum<-tapply(x,dat[,1],mean,na.rm=T);return(rep(dum,table(dat[,1])))})

nas<-which(is.na(ail))
ail[nas]<-mean.ail[nas] #set missing values to population mean
aic<-apply(ail-mean.ail,1,sum)
return(aic)
}
##################################################################################


#' @title Test for sex biased dispersal
#' @description Test whether one
#' sex disperses more than the other using the method described in 
#' \href{http://onlinelibrary.wiley.com/doi/10.1046/j.1365-294X.2002.01496.x/abstract}{Goudet etal. (2002)} 
#' @usage sexbias.test(dat,sex,nperm=NULL,test="mAIc",alternative="two.sided")
#' @param dat a data frame with n.locs+1 columns and n.inds rows
#' @param sex a vector containing the individual's sex
#' @param nperm the number of permutation to carry out
#' @param test one of "mAIc" (default), "vAIc","FIS" or "FST"
#' @param alternative one of "two.sided" (default),"less" or "greater"
#' @return call the function call
#' @return res the observation for each sex
#' @return statistic the observed statistic for the chosen test
#' @return p.value the p-value of the hypothesis
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @references \href{http://onlinelibrary.wiley.com/doi/10.1046/j.1365-294X.2002.01496.x/abstract}{Goudet J, Perrin N, Waser P (2002)} Tests for sex-biased dispersal 
#'  using bi-parentally inherited genetic markers 11, 1103:1114
#' @examples 
#'   data(crocrussula)
#'   sexbias.test(crocrussula$genot,crocrussula$sex)
#'   dat<-qn2.read.fstat(system.file("extdata","qn2_sex.dat",package="hierfstat"))
#'   sexbias.test(dat[[1]],sex=dat[[2]])
#'   \dontrun{
#'   sexbias.test(crocrussula$genot,crocrussula$sex,nperm=1000)
#'   sexbias.test(dat[[1]],sex=dat[[2]],nperm=100,test="FST",alternative="greater")
#'   }
#' @export
###################################################################################
sexbias.test<-function(dat,sex,nperm=NULL,test="mAIc",
                        alternative="two.sided"){
#calculate the corrected assignment index for all individuals
#dat is a genotype data frame
#sex is a vector of sex labels
#nperm is the number of permutations to be carried out in the permutation test
#2 tailed test on mean AIc for the time being
#TODO: 
#1. options for one tail tests
#2. include test on the variance and Fis / Fst (call wc)

  call<-match.call()
  if(!is.na(pmatch(test,"mAIc")))
    test<-"mAIc"
  TESTS<-c("mAIc","vAIc","FIS","FST")
  test<-pmatch(test,TESTS)
  if (is.na(test)) stop("Invalid test")
  if (test==-1) stop("Ambiguous test")
  
  if(!is.na(pmatch(alternative,"two.sided")))
    alternative<-"two.sided"
  ALTS<-c("two.sided","less","greater")
  alt<-pmatch(alternative,ALTS)
  if(is.na(alt)) stop("Invalid alternative")
  if(alt==-1) stop("Ambiguous alternative")
  

sexes<-names(table(sex))
sex1<-which(sex==sexes[1])
sex2<-which(sex==sexes[2])

if (test==3 | test==4){
  a<-wc(dat[sex1,])
  b<-wc(dat[sex2,])
  FIS<-c(a$FIS,b$FIS)
  names(FIS)<-sexes
  FST<-c(a$FST,b$FST)
  names(FST)<-sexes
}  

if (test==1 | test==2){
aic<-AIc(dat)
mAic1<-mean(aic[sex1])
mAic2<-mean(aic[sex2])

vAic1<-var(aic[sex1])
vAic2<-var(aic[sex2])

mAics<-data.frame(mAic1,mAic2)
names(mAics)<-sexes
vAics<-data.frame(vAic1,vAic2)
names(vAics)<-sexes
}

if (is.null(nperm)){
  if (test==1){
  tmAic<-t.test(aic[sex1],aic[sex2],alternative=alternative)
  res<-list(statistic=tmAic$statistic,p.value=tmAic$p.value)
  }
  if (test==2){
    vAic<-var.test(aic[sex1],aic[sex2],alternative=alternative)
    res<-list(statistic=vAic$statistic,p.value=vAic$p.value)
  }
  if (test==3 | test==4) stop("You must specify nperm for this test. Exiting")
}
else{
res<-list()
if (test==1){
tstats<-numeric(nperm)
tstats[1]<-t.test(aic[sex1],aic[sex2])$statistic

foo<-function(){
dum<-samp.within(dat[,1])
t.test(aic[dum][sex1],aic[dum][sex2])$statistic
}
tstats[-1]<-replicate(nperm-1,foo())
res$statistic<-tstats[1]
if (alt==1){res$p.value<-sum(abs(tstats)>=abs(tstats[1]))/nperm}
if (alt==2){res$p.value<-sum(tstats<=tstats[1])/nperm}
if (alt==3){res$p.value<-sum(tstats>=tstats[1])/nperm}
}

if (test==2){
  fstats<-numeric(nperm)
  fstats[1]<-var.test(aic[sex1],aic[sex2])$statistic
  
  foo<-function(){
    dum<-samp.within(dat[,1])
    return(var.test(aic[dum][sex1],aic[dum][sex2])$statistic)
  }
  fstats[-1]<-replicate(nperm-1,foo())
  res$statistic<-fstats[1]
  if (alt==1){res$p.value<-sum(abs(fstats)>=abs(fstats[1]))/nperm}
  if (alt==2){res$p.value<-sum(fstats<=fstats[1])/nperm}
  if (alt==3){res$p.value<-sum(fstats>=fstats[1])/nperm}
  
}
if (test==3){
  fstats<-numeric(nperm)
  fstats[1]<-FIS[1]-FIS[2]
  foo<-function(){
    dum<-samp.within(dat[,1])
    d1<-dum[sex1]
    d2<-dum[sex2]
    fis1<-wc(dat[d1,])$FIS
    fis2<-wc(dat[d2,])$FIS
    return(fis1-fis2)
  }
  fstats[-1]<-replicate(nperm-1,foo())
  res$statistic<-fstats[1]
  if (alt==1){res$p.value<-sum(abs(fstats)>=abs(fstats[1]))/nperm}
  if (alt==2){res$p.value<-sum(fstats<=fstats[1])/nperm}
  if (alt==3){res$p.value<-sum(fstats>=fstats[1])/nperm}
  
}
if (test==4){
  fstats<-numeric(nperm)
  fstats[1]<-FST[1]-FST[2]
  foo<-function(){
    dum<-samp.within(dat[,1])
    d1<-dum[sex1]
    d2<-dum[sex2]
    fst1<-wc(dat[d1,])$FST
    fst2<-wc(dat[d2,])$FST
    return(fst1-fst2)
  }
  fstats[-1]<-replicate(nperm-1,foo())
  res$statistic<-fstats[1]
  if (alt==1){res$p.value<-sum(abs(fstats)>=abs(fstats[1]))/nperm}
  if (alt==2){res$p.value<-sum(fstats<=fstats[1])/nperm}
  if (alt==3){res$p.value<-sum(fstats>=fstats[1])/nperm}
  
}
  
#if (plotit) {
#par(mfrow=c(3,1))
#yr<-c(floor(min(aic)),ceiling(max(aic)))
#hist(aic[sex=="M"],breaks=yr[1]:yr[2],xlab="AIc Males",main="")
#hist(aic[sex=="F"],breaks=yr[1]:yr[2],xlab="AIc Females",main="")
#hist(tstats,xlab="T-statistic",main="Null distribution of T-stats");abline(v=tstats[1],col="red",lwd=2)
#}
#2 tailed for the time being
}

return(list(call=match.call(),statistic=res$statistic,p.value=res$p.value))

}

 
