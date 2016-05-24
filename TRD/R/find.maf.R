find.maf<-function(sample){

dec=3

idx.case=grep(TRUE,(sample[,1]==1))
idx.ctrl=grep(TRUE,(sample[,1]==0))

case=sample[idx.case,]
ctrl=sample[idx.ctrl,]

# find MAF in whole sample
idx.m.0=length(grep(TRUE,sample[,2]==0))
idx.m.1=length(grep(TRUE,sample[,2]==1))
idx.m.2=length(grep(TRUE,sample[,2]==2))
idx.f.0=length(grep(TRUE,sample[,3]==0))
idx.f.1=length(grep(TRUE,sample[,3]==1))
idx.f.2=length(grep(TRUE,sample[,3]==2))
idx.b.0=length(grep(TRUE,sample[,4]==0))
idx.b.1=length(grep(TRUE,sample[,4]==1))
idx.b.2=length(grep(TRUE,sample[,4]==2))

# population size
psize=sum(idx.m.0,idx.m.1,idx.m.2,
          idx.f.0,idx.f.1,idx.f.2,
          idx.b.0,idx.b.1,idx.b.2)

# MAF (2 alleles per person)
p=((idx.m.1+idx.f.1+idx.b.1)+2*(idx.m.2+idx.f.2+idx.b.2))/(2*psize)

# Find MAF in child
psize.b=sum(idx.b.0,idx.b.1,idx.b.2)
p.b=(idx.b.1+2*idx.b.2)/2/psize.b

# Find MAF in parent
psize.p=sum(idx.m.0,idx.m.1,idx.m.2,idx.f.0,idx.f.1,idx.f.2)
p.p=(idx.m.1+2*idx.m.2+idx.f.1+2*idx.f.2)/2/psize.p

sample1=c(p.b,p.p,p)
names(sample1)=c('child','parents','trios')

if (length(idx.case)==0){
  cat('There are only control-trios in this sample.',fill=T)

  # find MAF in controls
  idx.m.0=length(grep(TRUE,ctrl[,2]==0))
  idx.m.1=length(grep(TRUE,ctrl[,2]==1))
  idx.m.2=length(grep(TRUE,ctrl[,2]==2))
  idx.f.0=length(grep(TRUE,ctrl[,3]==0))
  idx.f.1=length(grep(TRUE,ctrl[,3]==1))
  idx.f.2=length(grep(TRUE,ctrl[,3]==2))
  idx.b.0=length(grep(TRUE,ctrl[,4]==0))
  idx.b.1=length(grep(TRUE,ctrl[,4]==1))
  idx.b.2=length(grep(TRUE,ctrl[,4]==2))

  # population size
  psize.ctrl=sum(idx.m.0,idx.m.1,idx.m.2,
                 idx.f.0,idx.f.1,idx.f.2,
                 idx.b.0,idx.b.1,idx.b.2)

  # MAF (2 alleles per person)
  p.ctrl=((idx.m.1+idx.f.1+idx.b.1)+2*(idx.m.2+idx.f.2+idx.b.2))/(2*psize.ctrl)

  # Find MAF in child controls
  psize.ctrl.b=sum(idx.b.0,idx.b.1,idx.b.2)
  p.ctrl.b=(idx.b.1+2*idx.b.2)/2/psize.ctrl.b

  # Find MAF in parent controls
  psize.ctrl.p=sum(idx.m.0,idx.m.1,idx.m.2,idx.f.0,idx.f.1,idx.f.2)
  p.ctrl.p=(idx.m.1+2*idx.m.2+idx.f.1+2*idx.f.2)/2/psize.ctrl.p

  # Format data
  cases=c(NA,NA,NA)
  names(cases)=c('child','parents','trios')
  ctrls=c(p.ctrl.b,p.ctrl.p,p.ctrl)
  names(ctrls)=c('child','parents','trios')

} else if(length(idx.ctrl)==0){
  cat('There are only case-trios in this sample.',fill=T)

  # find MAF in cases
  idx.m.0=length(grep(TRUE,case[,2]==0))
  idx.m.1=length(grep(TRUE,case[,2]==1))
  idx.m.2=length(grep(TRUE,case[,2]==2))
  idx.f.0=length(grep(TRUE,case[,3]==0))
  idx.f.1=length(grep(TRUE,case[,3]==1))
  idx.f.2=length(grep(TRUE,case[,3]==2))
  idx.b.0=length(grep(TRUE,case[,4]==0))
  idx.b.1=length(grep(TRUE,case[,4]==1))
  idx.b.2=length(grep(TRUE,case[,4]==2))

  # population size
  psize.case=sum(idx.m.0,idx.m.1,idx.m.2,
                 idx.f.0,idx.f.1,idx.f.2,
                 idx.b.0,idx.b.1,idx.b.2)

  # MAF (2 alleles per person)
  p.case=((idx.m.1+idx.f.1+idx.b.1)+2*(idx.m.2+idx.f.2+idx.b.2))/(2*psize.case)

  # Find MAF in child cases
  psize.case.b=sum(idx.b.0,idx.b.1,idx.b.2)
  p.case.b=(idx.b.1+2*idx.b.2)/2/psize.case.b

  # Find MAF in parent cases
  psize.case.p=sum(idx.m.0,idx.m.1,idx.m.2,idx.f.0,idx.f.1,idx.f.2)
  p.case.p=(idx.m.1+2*idx.m.2+idx.f.1+2*idx.f.2)/2/psize.case.p

  # Format data
  cases=c(p.case.b,p.case.p,p.case)
  names(cases)=c('child','parents','trios')
  ctrls=c(NA,NA,NA)
  names(ctrls)=c('child','parents','trios')

} else {

  # find MAF in cases
  idx.m.0=length(grep(TRUE,case[,2]==0))
  idx.m.1=length(grep(TRUE,case[,2]==1))
  idx.m.2=length(grep(TRUE,case[,2]==2))
  idx.f.0=length(grep(TRUE,case[,3]==0))
  idx.f.1=length(grep(TRUE,case[,3]==1))
  idx.f.2=length(grep(TRUE,case[,3]==2))
  idx.b.0=length(grep(TRUE,case[,4]==0))
  idx.b.1=length(grep(TRUE,case[,4]==1))
  idx.b.2=length(grep(TRUE,case[,4]==2))

  # population size
  psize.case=sum(idx.m.0,idx.m.1,idx.m.2,
                 idx.f.0,idx.f.1,idx.f.2,
                 idx.b.0,idx.b.1,idx.b.2)

  # MAF (2 alleles per person)
  p.case=((idx.m.1+idx.f.1+idx.b.1)+2*(idx.m.2+idx.f.2+idx.b.2))/(2*psize.case)

  # Find MAF in child cases
  psize.case.b=sum(idx.b.0,idx.b.1,idx.b.2)
  p.case.b=(idx.b.1+2*idx.b.2)/2/psize.case.b

  # Find MAF in parent cases
  psize.case.p=sum(idx.m.0,idx.m.1,idx.m.2,idx.f.0,idx.f.1,idx.f.2)
  p.case.p=(idx.m.1+2*idx.m.2+idx.f.1+2*idx.f.2)/2/psize.case.p

  # find MAF in controls
  idx.m.0=length(grep(TRUE,ctrl[,2]==0))
  idx.m.1=length(grep(TRUE,ctrl[,2]==1))
  idx.m.2=length(grep(TRUE,ctrl[,2]==2))
  idx.f.0=length(grep(TRUE,ctrl[,3]==0))
  idx.f.1=length(grep(TRUE,ctrl[,3]==1))
  idx.f.2=length(grep(TRUE,ctrl[,3]==2))
  idx.b.0=length(grep(TRUE,ctrl[,4]==0))
  idx.b.1=length(grep(TRUE,ctrl[,4]==1))
  idx.b.2=length(grep(TRUE,ctrl[,4]==2))

  # population size
  psize.ctrl=sum(idx.m.0,idx.m.1,idx.m.2,
                 idx.f.0,idx.f.1,idx.f.2,
                 idx.b.0,idx.b.1,idx.b.2)

  # MAF (2 alleles per person)
  p.ctrl=((idx.m.1+idx.f.1+idx.b.1)+2*(idx.m.2+idx.f.2+idx.b.2))/(2*psize.ctrl)

  # Find MAF in child controls
  psize.ctrl.b=sum(idx.b.0,idx.b.1,idx.b.2)
  p.ctrl.b=(idx.b.1+2*idx.b.2)/2/psize.ctrl.b

  # Find MAF in parent controls
  psize.ctrl.p=sum(idx.m.0,idx.m.1,idx.m.2,idx.f.0,idx.f.1,idx.f.2)
  p.ctrl.p=(idx.m.1+2*idx.m.2+idx.f.1+2*idx.f.2)/2/psize.ctrl.p

  # Format data
  cases=c(p.case.b,p.case.p,p.case)
  names(cases)=c('child','parents','trios')
  ctrls=c(p.ctrl.b,p.ctrl.p,p.ctrl)
  names(ctrls)=c('child','parents','trios')
}

# Output

stat=list(signif(cases,dec),signif(ctrls,dec),signif(sample1,dec))
names(stat)=c('cases','ctrls','sample')
stat
}

