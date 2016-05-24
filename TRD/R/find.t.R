find.t<-function(sample){

  dec=3
  idx.case=grep(TRUE,(sample[,1]==1))
  idx.ctrl=grep(TRUE,(sample[,1]==0))

  case=sample[idx.case,]
  ctrl=sample[idx.ctrl,]

  if (dim(sample)[2]==5){

  if (length(idx.case)==0){
    cat('There are only control-trios in this sample.',fill=T)
    b2.m=tdt.cnt(ctrl)[1]
    c2.m=tdt.cnt(ctrl)[2]
    r2.m=b2.m/(b2.m+c2.m)

    b2.f=tdt.cnt(ctrl)[3]
    c2.f=tdt.cnt(ctrl)[4]
    r2.f=b2.f/(b2.f+c2.f)

    b2=tdt.cnt(ctrl)[5]
    c2=tdt.cnt(ctrl)[6]
    r2=b2/(b2+c2)

    r.out=c(NA,r2.m,NA,r2.f,NA,r2)

  } else if  (length(idx.ctrl)==0){
    cat('There are only case-trios in this sample.',fill=T)
    b1.m=tdt.cnt(case)[1]
    c1.m=tdt.cnt(case)[2]
    r1.m=b1.m/(b1.m+c1.m)

    b1.f=tdt.cnt(case)[3]
    c1.f=tdt.cnt(case)[4]
    r1.f=b1.f/(b1.f+c1.f)

    b1=tdt.cnt(case)[5]
    c1=tdt.cnt(case)[6]
    r1=b1/(b1+c1)

    r.out=c(r1.m,NA,r1.f,NA,r1,NA)

  } else {
    b1.m=tdt.cnt(case)[1]
    c1.m=tdt.cnt(case)[2]
    r1.m=b1.m/(b1.m+c1.m)

    b1.f=tdt.cnt(case)[3]
    c1.f=tdt.cnt(case)[4]
    r1.f=b1.f/(b1.f+c1.f)

    b2.m=tdt.cnt(ctrl)[1]
    c2.m=tdt.cnt(ctrl)[2]
    r2.m=b2.m/(b2.m+c2.m)

    b2.f=tdt.cnt(ctrl)[3]
    c2.f=tdt.cnt(ctrl)[4]
    r2.f=b2.f/(b2.f+c2.f)

    b1=tdt.cnt(case)[5]
    c1=tdt.cnt(case)[6]
    r1=b1/(b1+c1)

    b2=tdt.cnt(ctrl)[5]
    c2=tdt.cnt(ctrl)[6]
    r2=b2/(b2+c2)
    r.out=c(r1.m,r2.m,r1.f,r2.f,r1,r2)

  }} else if (dim(sample)[2]==4){
    cat('No parent-of-origin information is entered.',fill=T)
    if (length(idx.case)==0){
      cat('There are only control-trios in this data.',fill=T)
      b2=tdt.cnt(ctrl)[5]
      c2=tdt.cnt(ctrl)[6]
      r2=b2/(b2+c2)

      r.out=c(NA,NA,NA,NA,NA,r2)

    } else if  (length(idx.ctrl)==0){
      cat('There are only case-trios in this data.',fill=T)
      b1=tdt.cnt(case)[5]
      c1=tdt.cnt(case)[6]
      r1=b1/(b1+c1)

      r.out=c(NA,NA,NA,NA,r1,NA)

    } else {
      b1=tdt.cnt(case)[5]
      c1=tdt.cnt(case)[6]
      r1=b1/(b1+c1)

      b2=tdt.cnt(ctrl)[5]
      c2=tdt.cnt(ctrl)[6]
      r2=b2/(b2+c2)
      r.out=c(NA,NA,NA,NA,r1,r2)

    }}

  r.out=signif(r.out,dec)
  tcase=c(r.out[1],r.out[3],r.out[5])
  tctrl=c(r.out[2],r.out[4],r.out[6])
  names(tcase)=c('mother','father','parents')
  names(tctrl)=c('mother','father','parents')
  out=list(tcase,tctrl)
  names(out)=c('cases','ctrls')
  out
}
