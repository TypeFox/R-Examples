gcount<-function(sample){
  idx.case=grep(TRUE,(sample[,1]==1))
  idx.ctrl=grep(TRUE,(sample[,1]==0))

  case=sample[idx.case,]
  ctrl=sample[idx.ctrl,]

  sample.g2=length(grep(TRUE, sample[,4]==2))
  sample.g1=length(grep(TRUE, sample[,4]==1))
  sample.g0=length(grep(TRUE, sample[,4]==0))
  sample.m2=length(grep(TRUE, sample[,2]==2))
  sample.m1=length(grep(TRUE, sample[,2]==1))
  sample.m0=length(grep(TRUE, sample[,2]==0))
  sample.f2=length(grep(TRUE, sample[,3]==2))
  sample.f1=length(grep(TRUE, sample[,3]==1))
  sample.f0=length(grep(TRUE, sample[,3]==0))

  sample1=data.frame(matrix(c(sample.g0,sample.g1,sample.g2,
                              sample.m0,sample.m1,sample.m2,
                              sample.f0,sample.f1,sample.f2),
                            byrow=TRUE,ncol=3))
  names(sample1)=c('g0','g1','g2')
  rownames(sample1)=c('child','mother','father')


  if(length(idx.case)==0){
    cat('There are only control-trios in this sample.',fill=T)

    ctrl.g2=length(grep(TRUE, ctrl[,4]==2))
    ctrl.g1=length(grep(TRUE, ctrl[,4]==1))
    ctrl.g0=length(grep(TRUE, ctrl[,4]==0))
    ctrl.m2=length(grep(TRUE, ctrl[,2]==2))
    ctrl.m1=length(grep(TRUE, ctrl[,2]==1))
    ctrl.m0=length(grep(TRUE, ctrl[,2]==0))
    ctrl.f2=length(grep(TRUE, ctrl[,3]==2))
    ctrl.f1=length(grep(TRUE, ctrl[,3]==1))
    ctrl.f0=length(grep(TRUE, ctrl[,3]==0))

    cases=data.frame(matrix(rep(NA,9),
                            byrow=TRUE,ncol=3))
    names(cases)=c('g0','g1','g2')
    rownames(cases)=c('child','mother','father')


    ctrls=data.frame(matrix(c(ctrl.g0,ctrl.g1,ctrl.g2,
                              ctrl.m0,ctrl.m1,ctrl.m2,
                              ctrl.f0,ctrl.f1,ctrl.f2),
                            byrow=TRUE,ncol=3))
    names(ctrls)=c('g0','g1','g2')
    rownames(ctrls)=c('child','mother','father')


  } else if (length(idx.ctrl)==0){
    cat('There are only case-trios in this sample.',fill=T)

    case.g2=length(grep(TRUE, case[,4]==2))
    case.g1=length(grep(TRUE, case[,4]==1))
    case.g0=length(grep(TRUE, case[,4]==0))
    case.m2=length(grep(TRUE, case[,2]==2))
    case.m1=length(grep(TRUE, case[,2]==1))
    case.m0=length(grep(TRUE, case[,2]==0))
    case.f2=length(grep(TRUE, case[,3]==2))
    case.f1=length(grep(TRUE, case[,3]==1))
    case.f0=length(grep(TRUE, case[,3]==0))

    cases=data.frame(matrix(c(case.g0,case.g1,case.g2,
                              case.m0,case.m1,case.m2,
                              case.f0,case.f1,case.f2),
                            byrow=TRUE,ncol=3))
    names(cases)=c('g0','g1','g2')
    rownames(cases)=c('child','mother','father')


    ctrls=data.frame(matrix(rep(NA,9),
                            byrow=TRUE,ncol=3))
    names(ctrls)=c('g0','g1','g2')
    rownames(ctrls)=c('child','mother','father')


  } else {
    case.g2=length(grep(TRUE, case[,4]==2))
    case.g1=length(grep(TRUE, case[,4]==1))
    case.g0=length(grep(TRUE, case[,4]==0))
    case.m2=length(grep(TRUE, case[,2]==2))
    case.m1=length(grep(TRUE, case[,2]==1))
    case.m0=length(grep(TRUE, case[,2]==0))
    case.f2=length(grep(TRUE, case[,3]==2))
    case.f1=length(grep(TRUE, case[,3]==1))
    case.f0=length(grep(TRUE, case[,3]==0))

    ctrl.g2=length(grep(TRUE, ctrl[,4]==2))
    ctrl.g1=length(grep(TRUE, ctrl[,4]==1))
    ctrl.g0=length(grep(TRUE, ctrl[,4]==0))
    ctrl.m2=length(grep(TRUE, ctrl[,2]==2))
    ctrl.m1=length(grep(TRUE, ctrl[,2]==1))
    ctrl.m0=length(grep(TRUE, ctrl[,2]==0))
    ctrl.f2=length(grep(TRUE, ctrl[,3]==2))
    ctrl.f1=length(grep(TRUE, ctrl[,3]==1))
    ctrl.f0=length(grep(TRUE, ctrl[,3]==0))

    cases=data.frame(matrix(c(case.g0,case.g1,case.g2,
                              case.m0,case.m1,case.m2,
                              case.f0,case.f1,case.f2),
                            byrow=TRUE,ncol=3))
    names(cases)=c('g0','g1','g2')
    rownames(cases)=c('child','mother','father')


    ctrls=data.frame(matrix(c(ctrl.g0,ctrl.g1,ctrl.g2,
                              ctrl.m0,ctrl.m1,ctrl.m2,
                              ctrl.f0,ctrl.f1,ctrl.f2),
                            byrow=TRUE,ncol=3))
    names(ctrls)=c('g0','g1','g2')
    rownames(ctrls)=c('child','mother','father')
  }

  stat= list(cases,ctrls,sample1)
  names(stat)=c('cases','ctrls','sample')
  stat

}

